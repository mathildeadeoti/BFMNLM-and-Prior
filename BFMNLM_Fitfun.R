#################################################################################################### 
## Bayesian Nonlinear Mixed Effects Model ## 
#################################################################################################### 
## last updated in 2025-02-20 

#Function to get fitted values of Y for plotting
mufit <- function(x,#covariates matrix (N x p1) 
                  ind,#factor whose levels indicates the subjects or groups (N x 1)  
                  nlf,#nonlinear function with args: x1,...,xp1, beta1,...,betap, b1,...,bq                     
                  beta1,#p x 1 vector with initial values for beta 
                  D1,#q x q matrix with initial values for D 
                  bi # estimated random effects matrix(ncol=q1,nrow=m) with m<-n_distinct(ind) 
){   
  p1<-ncol(x) 
  p <- length(beta1) 
  q1 <- dim(D1)[2] 
  fittedval<- numeric(length(ind)) 
  ind_levels <- levels(ind) 
  for (i in seq_along(ind_levels)) { 
    seqi <- which(ind==ind_levels[i]) 
    xfiti <- matrix(x[seqi,],ncol=ncol(x)) 
    ub1 <- matrix(bi[i,],ncol=1) 
    for (j in seq_along(seqi))  
    { 
      farg <- as.list(c(xfiti[j,],beta1,ub1)) 
      names(farg) = c(paste("x",1:p1,sep=""),paste("beta",1:p,sep=""),paste("b",1:q1,sep="")) 
      formals(nlf) <- farg 
      fittedval[seqi[j]] <- nlf() 
    } 
  } 
  obj.out <- fittedval 
  obj.out
}