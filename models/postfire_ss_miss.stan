data {
  int<lower = 0> N_obs; //n observations
  int<lower = 0> N_mis; //n missing observations
  int<lower = 1, upper = N_obs + N_mis> ii_obs[N_obs]; //indexes of obs
  int<lower = 1, upper = N_obs + N_mis> ii_mis[N_mis]; //indexes of missing
  real y_obs[N_obs]; //ndvi data
  
  int<lower=0> T; //n time steps per pixel
  int<lower=0> J; // n pixels
  
  vector[J] envg1; //env covariate 1
  vector[J] envg2; //env covariate 2
  
  vector[J] z0; // Initial state value
  
  int fire[J,T]; //fire occurence matrix
}

transformed data {
  int<lower = 0> N = N_obs + N_mis; //n obs total (T*J)
}

parameters {
  
  real y_mis[N_mis]; //missing data parameters
  
  //env regressions
  vector[J] alpha;
  vector[J] gamma;
  vector[J] lambda;
  real alpha_mu;
  real gamma_b1;
  real gamma_b2;
  real lambda_b1;
  real lambda_b2;
  real<lower=0> tau_sq;
  real<lower=0> gamma_tau_sq;
  real<lower=0> lambda_tau_sq;
  real<lower=0> alpha_tau_sq;
  
  real<lower=0> sdp; // Standard deviation of the process equation
  real<lower=0> sdo; // Standard deviation of the observation equation
  
  real z[J, T]; // State time series
}

transformed parameters {
    
  real y[N]; //missing and observed
  y[ii_obs] = y_obs; //insert obs
  y[ii_mis] = y_mis; //insert missing
  
  vector[J] gamma_mu;
  vector[J] lambda_mu;
  real tau = sqrt(tau_sq);
  real gamma_tau = sqrt(gamma_tau_sq);
  real lambda_tau = sqrt(lambda_tau_sq);
  real alpha_tau = sqrt(alpha_tau_sq);
  
  for (j in 1:J){
    gamma_mu[j] = (envg1[j]*gamma_b1) + (envg2[j]*gamma_b2);
    lambda_mu[j] = (envg1[j]*lambda_b1) + (envg2[j]*lambda_b2);
  }
}

model {

 //TODO
 //pixel-wise process and obs sd
 //observation error on fire
 //posterior predictive
 //vectorize state update i.e no looping through T
  
  //process eq error
  sdp ~  normal(0, 1);
  //obs eq error
  sdo ~  normal(0, 1);
  //env regressions error
  gamma_tau ~ inv_gamma(0.01, 0.01);
  lambda_tau ~ inv_gamma(0.01, 0.01);
  alpha_tau ~ inv_gamma(0.01, 0.01);
  
  //env regression global effect
  alpha_mu ~ normal(0.15,3);
  gamma_b1 ~ normal(0,3);
  gamma_b2 ~ normal(0,3);
  lambda_b1 ~ normal(0,3);
  lambda_b2 ~ normal(0,3);
  
  //env regression pixels effect
  alpha ~ normal(alpha_mu, alpha_tau);
  gamma ~ normal(gamma_mu,gamma_tau); 
  lambda ~ normal(lambda_mu,lambda_tau); 

  //loop through pixels/groups
   for (j in 1:J) {
    int pos=1
    
    // Distribution for the first state
    z[j,1] ∼ normal(z0[j], sdp);
    
    // Distributions for all other states
    // logistic growth
    // fire swtiches state to alpha
    for(t in 2:T){
      z[j,t] ~ normal((z[j,t-1,] * lambda[j] * (1 - (z[j,t-1] / gamma[j]))*(1-fire[j,t-1])+(fire[j,t-1]*alpha[j])),sdp)
    }
    
    // Distributions for the observations
    
    //vector sampling
    //segment(y, pos, s[j]) ~ normal(z[j,~], sdo);
    
    //pixels sampling
    for(t in 1:TT){
      y[pos] ∼ normal(z[j,t], sdo);
      pos=pos+1
    }
  }
  
}