data {
  int<lower=0> TT; //n time steps per pixel
  int<lower=0> J; // n pixels
  int<lower=0> M; //N observations to simulate
  int<lower=0> N; //N observations total (obs + miss)
  
  
  int<lower = 0> N_obs; //n observations
  int<lower = 0> N_mis; //n missing observations
  int<lower = 1, upper = N_obs + N_mis> ii_obs[N_obs]; //indexes of obs
  int<lower = 1, upper = N_obs + N_mis> ii_mis[N_mis]; //indexes of missing
  real<lower = 0, upper = 1> y_obs[N_obs]; //ndvi data
  
  vector[J] envg1; //env covariate 1
  vector[J] envg2; //env covariate 2
  
  vector[J] z0; // Initial state value
  
  int fire[N]; //observed fire occurence matrix
  int firenew[M]; //predict fire occurence matrix
  
  real<lower = 0> doy_rad[J, TT]; // doy in radian for each pixel
  vector[M] doy_rad_new; //dates for predictions
    
  int<lower=0> eg_id; //pid for example pixel
}

parameters {
  real y_mis[N_mis]; //missing data parameters
  
  //seasonality
  //A*sin(x+phi) = A*cos(phi)*sin(x) + A*sin(phi)*cos(x)
  real<lower = -1, upper = 1> phase1[J]; // A*cos*phi
  real<lower = -1, upper = 1> phase2[J]; // A*sin*phi
  real phase_mu1;
  real phase_mu2;
  real<lower=0> phase_tau1;
  real<lower=0> phase_tau2;
  
  //env regressions
  vector<lower = 0, upper = 1> [J] alpha;
  vector<lower = 0, upper = 1> [J] gamma;
  vector<lower = 0> [J] lambda;
  
  //non regression params
  real alpha_mu;
  //regression params
  real gamma_b1;
  real gamma_b2;
  real lambda_b1;
  real lambda_b2;

  //variance priors
 real<lower=0> gamma_tau;
 real<lower=0> lambda_tau;
 real<lower=0> alpha_tau;
  
  real<lower=0> sdp; // Standard deviation of the process equation
  real<lower=0> sdo; // Standard deviation of the observation equation
  
  real<lower = 0, upper = 1> z[J, TT]; // State time series
  vector<lower = 0, upper = 1> [N] firefit; //fire occurence vector
}

transformed parameters {
  real y[N]; //missing and observed
  y[ii_obs] = y_obs; //insert obs
  y[ii_mis] = y_mis; //insert missing
  
  //env regression
  vector[J] gamma_mu;
  vector[J] lambda_mu;
  
  for (j in 1:J){
    gamma_mu[j] = (envg1[j]*gamma_b1) + (envg2[j]*gamma_b2);
    lambda_mu[j] = (envg1[j]*lambda_b1) + (envg2[j]*lambda_b2);
  }
}

model {
  //TODO: change normal priors on variance to exp(2) or exp(10)
  
  //process eq error
  sdp ~  exponential(2);
  //obs eq error
  //sdo ~  exp(2); weak prior
  sdo ~  normal(0.1, 0.001); //strng prior
  
  //fire measure prior
  firefit ~ normal(0, 1);
  
  //env regressions error
  gamma_tau ~  exponential(2);
  lambda_tau ~  exponential(2);
  alpha_tau ~  exponential(2);
  phase_tau1 ~  exponential(2);
  phase_tau2 ~  exponential(2);
  
  //env regression global effect
  alpha_mu ~ normal(0,1);
  phase_mu1 ~ normal(0,0.5);
  phase_mu2 ~ normal(0,0.5);

  gamma_b1 ~ normal(0,1);
  gamma_b2 ~ normal(0,1);
  lambda_b1 ~ normal(0,1);
  lambda_b2 ~ normal(0,1);
  
  //env regression pixels effect
  phase1 ~ normal(phase_mu1, phase_tau1);
  phase2 ~ normal(phase_mu2, phase_tau2);
  alpha ~ lognormal(alpha_mu, alpha_tau);
  gamma ~ lognormal(gamma_mu,gamma_tau); 
  lambda ~ lognormal(lambda_mu,lambda_tau); 
  
  int pos;
  
  //fire measuremet error model
   fire ~ bernoulli_logit(firefit);    // fire obsmeasurement model
   
  //loop through pixels/groups
   pos = 0;
   for (j in 1:J) {
     //loc for start of pixel ts

    // Distribution for the first state
    z[j,1] ~ normal(z0[j], sdp);
    pos = pos + 1;
    
    // Distributions for all other states
    // logistic growth with seasonality
    // fire swtiches state to alpha
    
    for(t in 2:TT){
      z[j,t] ~ normal(((z[j,t-1] + 
      (z[j,t-1] * lambda[j] * (1 - (z[j,t-1] / gamma[j]))))*(1-inv_logit(firefit[pos]))) +
      (inv_logit(firefit[pos])*alpha[j]) + 
      (phase1[j]*(sin(doy_rad[j,t-1]))+phase2[j]*(cos(doy_rad[j,t-1]))),sdp);
      
      pos = pos + 1;
    }
    // Distributions for the observations
    //vector sampling
    segment(y, pos-TT+1, TT) ~ normal(z[j], sdo);
  }
}

generated quantities {
  vector<lower = 0, upper = 1>[M] znew;
  znew[1] = z[eg_id,TT];
  
  //generate a sereis of predictions for a pixel
  for (m in 2:M){
        znew[m] = fmax(normal_rng(((znew[m-1] + (znew[m-1] * lambda[eg_id] * (1 - (znew[m-1] / gamma[eg_id]))))*(1-firenew[m-1])) +
        (firenew[m-1]*alpha[eg_id]) + 
        (phase1[eg_id]*(sin(doy_rad[eg_id,m-1]))+phase2[eg_id]*(cos(doy_rad[eg_id,m-1]))),sdp),0);
  
  }
}
