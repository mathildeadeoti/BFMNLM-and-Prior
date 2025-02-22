# Set directory and load functions
setwd("./")
####Data importation  and packages load

library(reshape2) 
library(tidyverse)
library(nlme)
library(mvtnorm) 

library(ggplot2)
library(ggeasy)
library(forcats)
library(ggpubr)

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE) # Allows to automatically save a bare version of a compiled Stan program to the hard disk so that it does not need to be recompiled.
if (file.exists(".RData")) file.remove(".RData")

##############################################################
######################################################################################*#######
######################################################################################*
############## Case of WA francofone data ##################################################
###################bayesian fit #########################
#######################################
#reading cumulative reported cases data  Y;
#run on of the first two next lines
row_data<-readxl::read_excel("Reported_case_COVID19_WA_data.xlsx");attach(row_data);head(row_data[,1:6])# direct importation
row_data<-read.table("clipboard", header=TRUE, sep = "\t");attach(row_data);head(row_data[,1:6])# copy and use this line to import
countries <- c("Benin","Burkina_Faso","Ivory_Coast","Guinea","Mali","Niger","Senegal")# Considered countries
row_data$Country.Region 

# melt dataset  
dati <- melt(row_data[-c(1,3,4)],id.vars = c("Country.Region")) ;head(dati)
dati<-droplevels(dati) 
dati$date<- as.Date.character(sub("X","",dati$variable),format="%m.%d.%y");head(dati) 
max(dati$date) 

dati <- filter(dati,dati$date<= "2020-10-01")  # data up to 2020-10-01
dati$Country.Region=as.factor(dati$Country.Region)


# creating new reported cases variable 
for (country in levels(dati$Country.Region)) dati$newcases[dati$Country.Region==country] =  
  c(0,diff(dati$value[dati$Country.Region==country])) ;head(dati)

## Plot to view how clustered are the data
windows()
ggplot(dati,aes(x=date,y=newcases,colour=Country.Region))+geom_line()+ 
  ylab("reported cases") + xlab("days since first reported case") 

############################################### 
#using data since first reported case for each country 
############################################### 
dat <- filter(dati,dati$value>=1) ;head(dat)
scalecte <- min(tapply(dat$newcases,dat$Country.Region,sd)) #constant used for numerical stability 
dat$newcasesT <- dat$newcases/scalecte #transformed variable 
dat <- dat[order(dat$Country.Region),] ;head(dat)  # Ordering countries

# creating variable day (different for each country) 
for (country in levels(dat$Country.Region)) dat$day[dat$Country.Region==country] = 1:sum(dat$Country.Region==country) ;head(dat)

## Plot to visualize countries separately at their different stade of epidemic curve
windows()
Fig1=ggplot(dat,aes(x=date,y=newcases))+geom_line()+facet_wrap(~Country.Region,scales = "free_y")+ 
  ylab("Cases") + xlab("Days since first reported cases") 

##################################################################################### 
# fitting 
## creating generalized logistic function and derivative as used in modeling process
derivlog3 <- function(x,loga,logb,logC) exp(logC)*exp(loga)* 
  exp(-exp(logC)*x)/((exp(logb)+exp(-exp(logC)*x))^2) #function without random effects, used for initial values 

logist3<- function(x1,beta1,beta2,beta3,beta4,b1,b2) {  
  exp(beta4+beta1+beta3+b2-exp(beta3+b2)*x1)/ 
    ((exp(beta2+b1)+exp(-exp(beta3+b2)*x1))^(exp(beta4)+1)) 
} 
der3logist <- deriv( ~   exp(beta4+beta1+beta3+b2-exp(beta3+b2)*x1)/ 
                       ((exp(beta2+b1)+exp(-exp(beta3+b2)*x1))^(exp(beta4)+1)), 
                     c("beta1","beta2","beta3","beta4","b1","b2"), function(x1,beta1,beta2,beta3,beta4,b1,b2){}) 

##################################################################################### 
# fitting nlme for getting initial values for Maximum likelihood approach of Schumasher et al. 2021 and Bayesian fitting
#The estimates from this Maximum likelihood approach of Schumasher et al. 2021 would serve as reference for prior setting
#This fitting would help to plot density and qqplot of residual
lognlme3.1=nlme(newcasesT~derivlog3(day,loga,logb,logC),data=dat, 
                fixed=loga+logb+logC~1,random=logb+logC~1|Country.Region, 
                start=c(loga=5,logb=-4,logC=-3), 
                control = list(msMaxIter=200,maxIter=200,returnObject=T))  
fixef(lognlme3.1) 
beta1=as.numeric(lognlme3.1$coefficients$fixed) 
sigmae=lognlme3.1$sigma^2 
D1=var(ranef(lognlme3.1)) 

###################################################################
############# RESIDUAL PLOT  ################
library(ggpubr)
resd=lognlme3.1$residuals;
windows() 
Fig2a=ggdensity(resd[,2], fill = "lightgray", xlab ="Residual" )
Fig2b=ggqqplot(resd[,2], xlab = "Quantiles of Standard Normal",ylab = "Quantiles of Residual"); shapiro.test(resd[,2])

#########################################
#Necessary function Maximum likelihood approach of Schumasher et al. 2021 and Bayesian fitting
source("smsn-nlmm.R") 
source("smn-nlmm.R") 
source("BFMNLM_Fitfun.R") 
#############################  fit NLMM based on Maximum likelihood approach of Schumasher et al. 2021
fitN<- EM.sim.NL(y=dat$newcasesT,x=matrix(dat$day),ind=dat$Country.Region, 
                 nlfder=der3logist,nlf=logist3,beta1=c(5,-4,-3,log(1)), 
                 sigmae=sigmae, 
                 D1=D1,distr="sn",nu=NULL,lb=NULL,lu=NULL, 
                 precisao=1e-4,max.iter=500,showiter=T,showerroriter=T) 
fitN$theta 
############################################################################## 
 
#### creating variable CLUSTER as stan code use numeric (different for each country) 
dat$cluster =as.numeric(factor(dat$Country.Region));head(dat)
datab <- list(N = length(dat$newcasesT), Y = dat$newcasesT,
              K_phi1=1, K_phi2=1,K_phi3=1, K_phi4=1,
              X_phi1=matrix(1,length(dat$newcasesT)) , X_phi2=matrix(1,length(dat$newcasesT)),X_phi3=matrix(1,length(dat$newcasesT)), X_phi4=matrix(1,length(dat$newcasesT)), C_1=dat$day,
              N_3=length(levels(dat$Country.Region)), N_2=length(levels(dat$Country.Region)), M_3=1, M_2=1,
              J_3=dat$cluster, J_2=dat$cluster,
              Z_3_phi3_1=rep(1,length(dat$newcasesT)), Z_2_phi2_1=rep(1,length(dat$newcasesT)),
              prior_only=0)

beta=c(2.477393, -3.480682, -3.318145) # nlme fit estimates to be used as initial values
##############################" stan model ##########################
##initial value
n_initfun1      <- function(chain_id) {list(b_phi1 = as.array( beta[1]), b_phi2 = as.array(beta[2]), b_phi3 = as.array(beta[3]),b_phi4 = as.array(log(1)),   sig2 = 1.5192, sd_2 = as.array(1), sd_3 = as.array(1), sd_32 = 0) }
st_initfun1     <- function(chain_id) {list(b_phi1 = as.array( beta[1]), b_phi2 = as.array(beta[2]), b_phi3 = as.array(beta[3]),b_phi4 = as.array(log(1)),   sig2 = 1.5192, sd_2 = as.array(1), sd_3 = as.array(1), sd_32 = 0, U_0 =rep(1,length(dat$newcasesT)),U_01 =rep(1,length(levels(dat$Country.Region))),nu = 3, delta =0, delta_3 = as.array(0), delta_2 = as.array(0), nu_b = as.array(3))}
snpk1_initfun  <- function(chain_id) {list(b_phi1 = as.array( beta[1]), b_phi2 = as.array(beta[2]), b_phi3 = as.array(beta[3]),b_phi4 = as.array(log(1)),   sig2 = 1.5192, sd_2 = as.array(1), sd_3 = as.array(1), sd_32 = 0,  phi_y1=0,phib1=0, phib2=0) }
##Fitting stan model
n_fiti <- stan_model(file ="./n_stan.stan")         #standard NLMM
sst_fiti <- stan_model(file ="./st_stan.stan")      #smsn-NLMM (flexible multilevel nonlinear model)
snp1_fiti <- stan_model(file ="./snpk1_stan.stan")  #snp-NLMM
### 
n_fit <- sampling(n_fiti, init = n_initfun1, data = datab,seed = 9891, chains = 3, iter = 400,warmup = 200, thin = 1, control = list(adapt_delta = .999, max_treedepth = 15));print(summary(n_fit, digits=2, pars = c("b_phi1", "b_phi2","b_phi3", "b_phi4", "sig2","sd_3","sd_2", "sd_32"))$summary)
sst_fit <- sampling(sst_fiti, init = st_initfun1, data = datab,seed = 9891, chains = 3, iter = 400,warmup = 200, thin = 1, control = list(adapt_delta = .999, max_treedepth = 15)); print(summary(sst_fit, digits=2, pars = c("b_phi1", "b_phi2","b_phi3", "b_phi4", "sig2","sd_3","sd_2", "sd_32","delta","delta_3","delta_2","nu","nu_b"))$summary)
snp1_fit <- sampling(snp1_fiti, init = snpk1_initfun, data = datab,seed = 9891, chains = 3, iter = 400,warmup = 200, thin = 1, control = list(adapt_delta = .999, max_treedepth = 15)); print(summary(snp1_fit, digits=2, pars = c("b_phi1", "b_phi2","b_phi3", "b_phi4", "sig2","sd_3","sd_2", "sd_32","phi_y1","phib1","phib2"))$summary)
save(n_fit, file = "./estimresult/resultn_wacase.RData");save(sst_fit, file = "./estimresult/resultsst_wacase.RData");save(snp1_fit, file = "./estimresult/resultsnp1_wacase.RData");
###
#load the ran flexible model to continue directly if one do not want to wait until the end of running. This is to continue with the code
load("./resultsst_wacase.RData");


##get fitting output 
######### Get nonlinear parameters
##st model
fixparamst=summary(sst_fit, digits=2, pars = c("b_phi1", "b_phi2","b_phi3", "b_phi4"))$summary[,1]; fixparamst
randparamst=matrix(summary(sst_fit, digits=2, pars = c("r_2_phi2_1", "r_3_phi3_1"))$summary[,1],nrow = length(levels(dat$Country.Region))); row.names(randparamst) <- levels(dat$Country.Region); randparamst
indpar_stfit <- exp(matrix(fixparamst,nrow=length(countries),ncol=4,byrow = T) + cbind(0,randparamst,0));indpar_stfit

#epidemic parameters for dynamic exploration 
tot_estim_cases_st= indpar_stfit[,1]/indpar_stfit[,2]^indpar_stfit[,4]*scalecte;tot_estim_cases_st# ;tot_estim_casesinf= indparinf[,1]/indparinf[,2]^indparinf[,4]*scalecte;tot_estim_casessup= indparsup[,1]/indparsup[,2]^indparsup[,4]*scalecte  
peak_time_st = -log(indpar_stfit[,2]/indpar_stfit[,4])/indpar_stfit[,3];peak_time_st# ;peak_timeinf = -log(indparinf[,2]/indparinf[,4])/indparinf[,3];peak_timesup = -log(indparsup[,2]/indparsup[,4])/indparsup[,3]
G=7#################Generation interval: 
reprodst=exp(indpar_stfit[,3]*G);reprodst

## fitted plot
fitsampst<- mufit(x=matrix(dat$day),ind=dat$Country.Region, 
                  nlf=logist3,beta1=fixparamst, 
                  D1=var(ranef(lognlme3.1)),#q x q matrix with initial values for D
                  bi=randparamst) 
fitlogfitst<- data.frame(value=dat$newcases,day=dat$day, 
                         date = dat$date,Country.Region=dat$Country.Region, 
                         fitted=fitsampst*scalecte)
windows()
ggplot(fitlogfitst,aes(x=day,y=value))+geom_line()+geom_point()+ylab("deaths")+ 
geom_line(aes(y=fitted,color=I("blue"))) + facet_wrap(~Country.Region,scales = "free_y") 

##snp1 case
fixparamsnp1=summary(snp1_fit, digits=2, pars = c("b_phi1", "b_phi2","b_phi3", "b_phi4"))$summary[,1]; fixparamsnp1
randparamsnp1=matrix(summary(snp1_fit, digits=2, pars = c("r_2_phi2_1", "r_3_phi3_1"))$summary[,1],nrow = length(levels(dat$Country.Region))); row.names(randparamsnp1) <- levels(dat$Country.Region); randparamsnp1
indpar_snp1fit <- exp(matrix(fixparamsnp1,nrow=length(countries),ncol=4,byrow = T) + cbind(0,randparamsnp1,0));indpar_snp1fit

#epidemic parameters for dynamic exploration 
tot_estim_cases_snp1= indpar_snp1fit[,1]/indpar_snp1fit[,2]^indpar_snp1fit[,4]*scalecte ;tot_estim_cases_snp1# ;tot_estim_casesinf= indparinf[,1]/indparinf[,2]^indparinf[,4]*scalecte;tot_estim_casessup= indparsup[,1]/indparsup[,2]^indparsup[,4]*scalecte  
peak_time_snp1 = -log(indpar_snp1fit[,2]/indpar_snp1fit[,4])/indpar_snp1fit[,3];peak_time_snp1# ;peak_timeinf = -log(indparinf[,2]/indparinf[,4])/indparinf[,3];peak_timesup = -log(indparsup[,2]/indparsup[,4])/indparsup[,3]
G=7#################Generation interval: 
reprodsnp1=exp(indpar_snp1fit[,3]*G);reprodsnp1

## fitted plot
fitsampsnp1<- mufit(x=matrix(dat$day),ind=dat$Country.Region, 
                    nlf=logist3,beta1=fixparamsnp1, 
                    D1=var(ranef(lognlme3.1)),#q x q matrix with initial values for D
                    bi=randparamsnp1) 
fitlogfitsnp1<- data.frame(value=dat$newcases,day=dat$day, 
                           date = dat$date,Country.Region=dat$Country.Region,
                           fitted=fitsampsnp1*scalecte)
windows()
ggplot(fitlogfitsnp1,aes(x=day,y=value))+geom_line()+geom_point()+ylab("deaths")+
  geom_line(aes(y=fitted,color=I("blue"))) + facet_wrap(~Country.Region,scales = "free_y") 
windows()
ggplot(fitlogfitsnp1,aes(x=day,y=fitted))+geom_line()+geom_point()+ylab("deaths") + facet_wrap(~Country.Region,scales = "free_y")

##normal model
fixparamn=summary(n_fit, digits=2, pars = c("b_phi1", "b_phi2","b_phi3", "b_phi4"))$summary[,1]; fixparamn
randparamn=matrix(summary(n_fit, digits=2, pars = c("r_2_phi2_1","r_3_phi3_1" ))$summary[,1],nrow = length(levels(dat$Country.Region))); row.names(randparamn) <- levels(dat$Country.Region); randparamn
indpar_nfit <- exp(matrix(fixparamn,nrow=length(countries),ncol=4,byrow = T) + cbind(0,randparamn,0));indpar_nfit

#epidemic parameters for dynamic exploration 
tot_estim_cases_n= indpar_nfit[,1]/indpar_nfit[,2]^indpar_nfit[,4]*scalecte ;tot_estim_cases_n;#tot_estim_casesinf= indparinf[,1]/indparinf[,2]^indparinf[,4]*scalecte;tot_estim_casessup= indparsup[,1]/indparsup[,2]^indparsup[,4]*scalecte  
peak_time_n = -log(indpar_nfit[,2]/indpar_nfit[,4])/indpar_nfit[,3] ;peak_time_n;#peak_timeinf = -log(indparinf[,2]/indparinf[,4])/indparinf[,3];peak_timesup = -log(indparsup[,2]/indparsup[,4])/indparsup[,3]
G=7#################Generation interval: 
reprodn=exp(indpar_nfit[,3]*G);reprodn

## fitted plot
fitsampn<- mufit(x=matrix(dat$day),ind=dat$Country.Region, 
                 nlf=logist3,beta1=fixparamn, 
                 D1=var(ranef(lognlme3.1)),#q x q matrix with initial values for D
                 bi=randparamn) 
fitlogfitn<- data.frame(value=dat$newcases,day=dat$day, 
                        date = dat$date,Country.Region=dat$Country.Region, 
                        fitted=fitsampn*scalecte)
windows()
ggplot(fitlogfitn,aes(x=day,y=value))+geom_line()+geom_point()+ylab("deaths")+ 
  geom_line(aes(y=fitted,color=I("blue"))) + facet_wrap(~Country.Region,scales = "free_y") 

