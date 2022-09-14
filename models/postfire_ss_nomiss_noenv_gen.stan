data {
  int<lower=0> N; //n obs total (T*J)
  int<lower=0> TT; //n time steps per pixel
  int<lower=0> J; // n pixels
  int<lower=0> M; //N observations to simulate
  //int<lower=1,upper=J> pid[N];
  vector<lower = 0, upper = 1> [N] y; // obs
//  vector[J] envg1; //env covariate 1
//  vector[J] envg2; //env covariate 2
  
  vector[J] z0; // Initial state value
  
  int fire[J,TT]; //fire occurence matrix
  int firenew[M]; //predict fire occurence matrix
}
parameters {
  //env regressions
  vector<lower = 0, upper = 1> [J] alpha;
  vector<lower = 0, upper = 1> [J] gamma;
  vector<lower = 0> [J] lambda;
  real alpha_mu;
  real gamma_mu;
  real lambda_mu;
//  real gamma_b1;
//  real gamma_b2;
//  real lambda_b1;
//  real lambda_b2;
//  real<lower=0> tau_sq;
//  real<lower=0> gamma_tau_sq;
//  real<lower=0> lambda_tau_sq;
//  real<lower=0> alpha_tau_sq;
 real<lower=0> gamma_tau;
 real<lower=0> lambda_tau;
 real<lower=0> alpha_tau;
  
  real<lower=0> sdp; // Standard deviation of the process equation
  real<lower=0> sdo; // Standard deviation of the observation equation
  
  real<lower = 0, upper = 1> z[J, TT]; // State time series
}

transformed parameters {
//  vector[J] gamma_mu;
//  vector[J] lambda_mu;
//  real tau = sqrt(tau_sq);
//  real gamma_tau = sqrt(gamma_tau_sq);
//  real lambda_tau = sqrt(lambda_tau_sq);
//  real alpha_tau = sqrt(alpha_tau_sq);
  
//  for (j in 1:J){
//    gamma_mu[j] = (envg1[j]*gamma_b1) + (envg2[j]*gamma_b2);
//    lambda_mu[j] = (envg1[j]*lambda_b1) + (envg2[j]*lambda_b2);
//  }
}

model {

 //TODO
 //pixel-wise process and obs sd
 //observation error on fire
 //posterior predictive
 //vectorize state update i.e no looping through T
  
  //process eq error
  sdp ~  normal(0, 0.5);
  //obs eq error
  sdo ~  normal(0, 0.5);
  //env regressions error
  //gamma_tau ~ inv_gamma(0.01, 0.01);
  //lambda_tau ~ inv_gamma(0.01, 0.01);
  //alpha_tau ~ inv_gamma(0.01, 0.01);
  gamma_tau ~  normal(0, 0.5);
  lambda_tau ~  normal(0, 0.5);
  alpha_tau ~  normal(0, 0.5);
  
  //env regression global effect
  alpha_mu ~ normal(0.15,0.5);
  lambda_mu ~ normal(0.15,0.5);
  gamma_mu ~ normal(0.3,0.5);
//  gamma_b1 ~ normal(0,3);
//  gamma_b2 ~ normal(0,3);
//  lambda_b1 ~ normal(0,3);
//  lambda_b2 ~ normal(0,3);
  
  //env regression pixels effect
  alpha ~ normal(alpha_mu, alpha_tau);
  gamma ~ normal(gamma_mu,gamma_tau); 
  lambda ~ normal(lambda_mu,lambda_tau); 
  int pos;
  pos = 1;
  //loop through pixels/groups
   for (j in 1:J) {
    
    // Distribution for the first state
    z[j,1] ~ normal(z0[j], sdp);
    
    // Distributions for all other states - logistic growth
    // fire swtiches state to alpha
    for(t in 2:TT){
      z[j,t] ~ normal(((z[j,t-1] + (z[j,t-1] * lambda[j] * (1 - (z[j,t-1] / gamma[j]))))*(1-fire[j,t-1]))+(fire[j,t-1]*alpha[j]),sdp);
    }
    // Distributions for the observations
    //vector sampling
    //segment(y, pos, s[j]) ~ normal(z[~,j], sigma);
    
    //pixels sampling
    for(t in 1:TT){
      y[pos] ~ normal(z[j,t], sdo);
      pos=pos+1;
    }
  }
}

generated quantities {
  vector<lower = 0, upper = 1>[M] znew;
  znew[1] = z[16,TT];
  
  for (m in 2:M){
        znew[m] = normal_rng(((znew[m-1] + (znew[m-1] * lambda[16] * (1 - (znew[m-1] / gamma[16]))))*(1-firenew[m-1]))+(firenew[m-1]*alpha[16]),sdp);
  
  }
}
