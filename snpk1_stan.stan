// generated with brms 2.13.5
functions {
}
data {
  int<lower=1> N;  // number of observations
  vector[N] Y;  // response variable
  int<lower=1> K_phi1;  // number of population-level effects
  matrix[N, K_phi1] X_phi1;  // population-level design matrix
  int<lower=1> K_phi2;  // number of population-level effects
  matrix[N, K_phi2] X_phi2;  // population-level design matrix
int<lower=1> K_phi3;  // number of population-level effects
  matrix[N, K_phi1] X_phi3;  // population-level design matrix
  int<lower=1> K_phi4;  // number of population-level effects
  matrix[N, K_phi2] X_phi4;  // population-level design matrix

 // covariate vectors for non-linear functions
  vector[N] C_1;
  // data for group-level effects of ID 2
  int<lower=1> N_2;  // number of grouping levels
  int<lower=1> M_2;  // number of coefficients per level
    int<lower=1> J_2[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_2_phi2_1;
 // data for group-level effects of ID 3
  int<lower=1> N_3;  // number of grouping levels
  int<lower=1> M_3;  // number of coefficients per level
  int<lower=1> J_3[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_3_phi3_1;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
 real halfpi;
  halfpi=((exp(lgamma(.5)))^2)/2;
  }
parameters {
 vector[K_phi1] b_phi1;  // population-level effects
 vector[K_phi2] b_phi2;  // population-level effects
 vector[K_phi3] b_phi3;  // population-level effects
 vector[K_phi4] b_phi4;  // population-level effects

  real<lower=0,upper=5> sig2;  // residual SD
   real<lower=-halfpi,upper=halfpi> phib1;real<lower=-halfpi,upper=halfpi> phib2;real<lower=-halfpi,upper=halfpi> phi_y1;// snp shape parameter 
               vector<lower=0,upper=5>[M_2] sd_2;  // group-level standard deviations
                vector[N_2] z_2[M_2];  // standardized group-level effects
                vector<lower=0,upper=5>[M_3] sd_3;  // group-level standard deviations
               vector[N_3] z_3[M_3];  // standardized group-level effects
                real  sd_32;  // random effect correlation
}
transformed parameters {
  vector[N_2] r_2_phi2_1;  // actual group-level effects
  vector[N_3] r_3_phi3_1;  // actual group-level effects
  vector[N_3] pkb;vector[3] pb;real Eb3; real Eb2;  //for bi : polynom , vector of SNP parameters,location shifts used to zero means
  vector[2] py;real Ey;  //for y :  polynom , vector of SNP parameters,location shifts used to zero means
  
  pb[1]= cos(phib1);pb[2]=sin(phib1)*cos(phib2);pb[3]=sin(phib1)*sin(phib2);
  pkb = fabs(pb[1]+pb[2]*z_3[1]+pb[3]*z_2[1]);
 // Eb1=(-sd_3[1])*((pb'*A1)*pb);Eb2=(-sd_2[1])*((pb'*A2)*pb);
 Eb2=(-sd_2[1])*(2*pb[1]*pb[3]);
 Eb3=(-sd_3[1])*(2*pb[1]*pb[2]);
  
      py[1]= cos(phi_y1);py[2]=sin(phi_y1);

    // Ey=(-sqrt(sigma))*((py'*A)*py);
     Ey=(sqrt(sig2))*(2*py[1]*py[2]);

  r_2_phi2_1 = ((sd_2[1] * (z_2[1])+Eb2));
  r_3_phi3_1 = ((sd_32 * (z_2[1]) + sd_3[1] * (z_3[1])+Eb3));
}
model {
  // initialize linear predictor term
  vector[N] nlp_phi1 = X_phi1 * b_phi1;
  // initialize linear predictor term
  vector[N] nlp_phi2 = X_phi2 * b_phi2;
 // initialize linear predictor term
  vector[N] nlp_phi3 = X_phi3 * b_phi3;
  // initialize linear predictor term
  vector[N] nlp_phi4 = X_phi4 * b_phi4;

  // initialize non-linear predictor term
       vector[N] mu;  
       vector [N] z_y;
       vector [N] pky;
  for (n in 1:N) {
    // add more terms to the linear predictor
    nlp_phi2[n] += r_2_phi2_1[J_2[n]] * Z_2_phi2_1[n];
  }
        for (n in 1:N) {
    // add more terms to the linear predictor
    nlp_phi3[n] += r_3_phi3_1[J_3[n]] * Z_3_phi3_1[n];
  }
  for (n in 1:N) {
     //compute non-linear predictor values     
    mu[n] = exp(nlp_phi4[n] + nlp_phi1[n] + nlp_phi3[n] - exp(nlp_phi3[n]) * C_1[n]) / ((exp(nlp_phi2[n]) + exp( - exp(nlp_phi3[n]) * C_1[n])) ^ (exp(nlp_phi4[n]) + 1)); // mean function
    
    z_y[n] = (Y[n] - mu[n] + Ey) / sqrt(sig2);
    pky[n] = fabs(py[1]+py[2]*z_y[n]);
  }
  
  // priors including all constants
  target += 2 * sum(log(pkb));
  target += - 0.5 * log(sig2) * N + 2 * sum(log(pky));
  target += normal_lpdf(b_phi1 | 2.466, 0.302);
  target += normal_lpdf(b_phi2 | -7.784, 0.695);
  target += normal_lpdf(b_phi3 | -2.985, 0.249);
  target += normal_lpdf(b_phi4 | -0.518, 0.034);
   target += uniform_lpdf(phi_y1 |-halfpi,halfpi);
   target += uniform_lpdf(phib1 | -halfpi,halfpi);
  target += uniform_lpdf(phib2 | -halfpi,halfpi);
  target += student_t_lpdf(sig2 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
    target += student_t_lpdf(sd_2 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
   target += std_normal_lpdf(z_2[1]);
  target += student_t_lpdf(sd_3 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
    target += std_normal_lpdf(z_3[1]);
   target += normal_lpdf(sd_32 | 0, 10^3);
 // likelihood including all constants
  if (!prior_only) {
    target += normal_lpdf(z_y | 0, 1);
  }
}
generated quantities {
 }
