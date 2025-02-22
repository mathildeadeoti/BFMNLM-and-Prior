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
  real c=sqrt(2) * exp(-lgamma(.5));
     }
parameters {
  vector[K_phi1] b_phi1;  // population-level effects
  vector[K_phi2] b_phi2;  // population-level effects
  vector[K_phi3] b_phi3;  // population-level effects
  vector[K_phi4] b_phi4;  // population-level effects
  real<lower=0,upper=5> sig2;  // residual scale SD
  real delta;  // residual skew
   real nu;  // parameter of response mixture variable
  vector[N] T_0;  // Half normal
  vector[N] U_0;  // nu function for response
  vector<lower=0,upper=5>[M_2] sd_2;  // group-level standard deviations
  vector[N_2] z_2[M_2];  // standardized group-level effects
    vector<lower=0,upper=5>[M_3] sd_3;  // group-level standard deviations
  vector[N_3] z_3[M_3];  // standardized group-level effects
  vector[M_3] nu_b;  // parameter of random effects mixture variable 
  vector[M_2] delta_2;  // group-level
vector[M_3] delta_3;  // group-level skew
  vector[N_3] T_01[M_3];  // Half normal
  vector[N_3] U_01;  // nu function for random effect
  real sd_32;  // random effect correlation
}
transformed parameters {
  vector[N_2] r_2_phi2_1;  // actual group-level effects
  vector[N_3] r_3_phi3_1;  // actual group-level effects
  real<lower=0> ctaue1;real<lower=0> cU_0e1;
  ctaue1 = c *sqrt(nu_b[1]/2)*exp(lgamma((nu_b[1]-1)/2)-lgamma((nu_b[1])/2)) ;
  cU_0e1 = c *sqrt(nu/2)*exp(lgamma((nu-1)/2)-lgamma((nu)/2)) ;
 r_2_phi2_1 = ((sd_2[1] * (z_2[1]) + (fabs(T_01[1]) * delta_2[1]))/sqrt(U_01[1]) - ctaue1 * delta_2[1]);
 r_3_phi3_1 = ((sd_32 * (z_2[1]) + sd_3[1] * (z_3[1]) + (fabs(T_01[1]) * delta_3[1]))/sqrt(U_01[1]) - ctaue1 * delta_3[1]);
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
   // actualized variance covariance of random effect
  for (n in 1:N) {
    // add more terms to the linear predictor
    nlp_phi2[n] += r_2_phi2_1[J_2[n]] * Z_2_phi2_1[n];
  }
       for (n in 1:N) {
    // add more terms to the linear predictor
    nlp_phi3[n] += r_3_phi3_1[J_3[n]] * Z_3_phi3_1[n];
  }
  for (n in 1:N) {
    // compute non-linear predictor values
     mu[n] = exp(nlp_phi4[n] + nlp_phi1[n] + nlp_phi3[n] - exp(nlp_phi3[n]) * C_1[n]) / ((exp(nlp_phi2[n]) + exp( - exp(nlp_phi3[n]) * C_1[n])) ^ (exp(nlp_phi4[n]) + 1)) + delta * (fabs(T_0[n])/sqrt(U_0[n]) - cU_0e1);
  }
  // priors including all constants
  target += normal_lpdf(b_phi1 | 2.466, 0.302);
  target += normal_lpdf(b_phi2 | -7.784, 0.695);
  target += normal_lpdf(b_phi3 | -2.985, 0.249);
  target += normal_lpdf(b_phi4 | -0.518, 0.034);
  target += student_t_lpdf(sig2 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += student_t_lpdf(sd_2 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
    target += std_normal_lpdf(z_2[1]);
  target += student_t_lpdf(sd_3 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(z_3[1]);
 target += normal_lpdf(sd_32 | 0, 10^3);
  target += std_normal_lpdf(T_0);
  target += gamma_lpdf(nu | 2,0.1);
  target += gamma_lpdf(nu_b | 2,0.1);
   target += gamma_lpdf(U_0| nu,  nu);
  target += normal_lpdf(delta | 0, 10^3);
  target += normal_lpdf(delta_3 | 0, 10^3);
  target += normal_lpdf(delta_2 | 0, 10^3);
  target += std_normal_lpdf(T_01[1]);
 target += gamma_lpdf(U_01[1]| nu_b,  nu_b);
  // likelihood including all constants
  if (!prior_only) {
    target += normal_lpdf(Y | mu, sig2/sqrt(U_0[1]));
  }
}
generated quantities {
}
