
// Written by Zachary Susswein for the Pandemic Prevention Institute
data {
  int<lower=1> N;        // total number of observations (days)
  int<lower = 2> K;      // number of countries
  int<lower=2> ncat;     // number of multinomial categories (variants)
  int Y[K, N, ncat];     // response array, count of country x day x variant
  int trials[K, N];      // total number of trials (seqs) in the country x day
  int prior_only;        // should the likelihood be ignored?
  int num_days_nowcast;  // number of days to project
}

transformed data{
  
 ///////////////////////////
 // Initialize variables
 ///////////////////////////

  // The number of cases (i.e., variants) being compared to the base case
  // base case should be most prevalent variant, for numerical stability
  int<lower = 1> ncomp;
  
  // Params to z-score time
  // puts most predictors on unit scale (or pretty close), which helps the 
  // sampler and simplifies many priors
  vector[N+num_days_nowcast] t_scale;
  vector[N+num_days_nowcast] t;
  real mean_t;
  real sd_t;
  
 ///////////////////////////
 // Evaluate
 ///////////////////////////

  // hacky, but only way to set up vector of timepoints within Stan
  for(i in 1:N+num_days_nowcast){

      t[i] = i;
  }
  
  // Center and scale time  based on the number of _observed_ days that the 
  // model is being fit to, not the number of observed days + the number of days 
  // being nowcast (i.e., the length of t). This should only make a small 
  // difference, but prevents weirdness from occurring if we eventually start 
  // trying to forecast and allows us to still easily generate fitted nowcasts.
  mean_t = mean(t[1:N]);
  sd_t = sd(t[1:N]);
  
  // z-score
  t_scale = (t - mean_t) / sd_t;
  
  // Number of comparisons to base case (for multinomial)
  ncomp = ncat - 1;
}

parameters {
 
  ///////////////////////////
  // Growth rate (i.e., slope)
  ///////////////////////////
  
  // Means of growth rates
  vector[ncomp] mu_r;  // vector of average relative variant growth rate across countries
  
  // Hierarchical parameters over distribution of relative growth rate means
  real mu_hier_raw;
  real<lower = 0> sigma_hierarchical;

  // Covariance/correlation of growth rates
  // https://mc-stan.org/docs/2_29/stan-users-guide/multivariate-hierarchical-priors.html
  cholesky_factor_corr[ncomp] L_Omega; // Cholesky factor of correlation matrix
  vector<lower=0>[ncomp] tau; // scale of var-covar mat. -- currently het. variance
  matrix[ncomp, K] z_omega; // matrix of iid normal, to non-center the MVN
  
  ///////////////////////////
  // Intercept
  ///////////////////////////
  
  // Unlike slopes, intercept distributions are independent normal 
  // distributions. Each variant has its own mean and standard deviation for 
  // its intercept. This allows variant introduction patterns to be different 
  // for each variant.
  vector[ncomp] mu_b0;                // Mean of the intercept vector
  vector<lower = 0>[ncomp] sigma_b0;  // Sd of the intercept vector
  matrix[ncomp, K] z_b0;              // iid std normal to non-center b0
  
}

transformed parameters {

  ///////////////////////////
  // Parameters
  ///////////////////////////

  real mu_hier_scaled;       // mean relative growth rate over all variants
  matrix[ncomp, K] b0;  // intercepts matrix, one intercept per variant-country
  matrix[ncomp, K] r_scaled;  // relative variant growth rate in country k

  ///////////////////////////
  // Rescale and non-center
  ///////////////////////////

  // Rescale mu_hier_scaled so that sampler works on unit scale
  mu_hier_scaled = mu_hier_raw * 0.1;
  
  // Non-centered parameterization of intercepts, vectorized over countries
  // implies b0[k] ~ N(mu_b0[k], sigma_b0[k]) where k \in 1:K
  // use element-wise matrix multiplication between sigma_b0 and z_b0
  b0 = rep_matrix(mu_b0, K) + rep_matrix(sigma_b0, K) .* z_b0;

  // Multivariate non-centered paramaterization, implies:
  // r_scaled[,k] ~ MVN(mu_r, Sigma)
  // with Cholesky factor of Sigma, L_Sigma. 
  // Therefore, Sigma = L_sigma * L_Sigma^T
  // and L_Sigma = tau * L_Omega, 
  // where L_Omega is the Cholesky factor of the correlation matrix.
  // Useful approach for Neal's funnel in multivariate case.
  // see: https://mc-stan.org/docs/2_29/reference-manual/cholesky-factors-of-correlation-matrices-1.html
  // and https://mc-stan.org/docs/2_18/stan-users-guide/multivariate-hierarchical-priors-section.html
  r_scaled = rep_matrix(mu_r, K) + diag_pre_multiply(tau, L_Omega) * z_omega;

}

model {

  ///////////////////////////
  // Evaluate likelihood
  ///////////////////////////

  // Likelihood, w/ switch to simulate from priors/actually evaluate likelihood
  if (!prior_only) {
      
    // Initialize country-specific vectors from matrices of params
    vector[ncomp] b0_k; // country-specific intercepts (1 for each variant)
    vector[ncomp] r_k; // country-specific growth advantage (1 for each variant)
    // Initialize linear predictor term
    // Fitted relative log odds of observing variant i in country k
    // equivalent to eta = log(growth rate of X) - log(growth rate of reference)
    vector[ncat] eta;
      
    // Iterate over countries
    for(k in 1:K){
      
      // Pull out country-specific params as vectors
      b0_k = to_vector(b0[,k]);
      r_k = to_vector(r_scaled[,k]);
        
      // Within country, iterate over days
      for(i in 1:N){
        
        // Don't evaluate days with no sequences
        if(trials[k, i] > 0){
    
          // Base case for reference variant
          eta[1] = 0;
          
          // Linear predictor, relative to base case eta[1]
          eta[2:ncat] =  b0_k + (r_k * t_scale[i]);
        
          // Softmax maps $eta \in R^{ncat}$ to an ncat-length simplex -- in 
          // other words, the estimated variant prevalences.  
          // https://mc-stan.org/docs/2_26/functions-reference/softmax.html
          Y[k, i] ~ multinomial(softmax(eta));
        
        }
      }
    }
  }
  
  ///////////////////////////
  // Priors
  ///////////////////////////
  
  // Vectorized priors on independent hierarchical intercept means and sds.
  // This approach implies that different variants are introduced at 
  // independent times (no structure in the means) and that the variants 
  // spread differently between countries (no structure in variances).
  mu_b0 ~ student_t(3, -5, 5);
  sigma_b0 ~ normal(2, 1);
  to_vector(z_b0) ~ std_normal();
  
  // Hierarchical distribution over average variant slopes.
  // Shrinks growth rates for variants with only a few global sequences to 
  // global mean relative growth rate to avoid overestimating new variants.
  mu_r ~ normal(mu_hier_scaled, sigma_hierarchical); 
  
  // Hierarchical mean and sd of variant slope means
  // Prior implies that a given variant has ~97% prob. of being less fit than 
  // the current dominant variant (the base case).
  mu_hier_raw ~ normal(0, .5); // scale mu_hier_scaled to be on unit scale
  sigma_hierarchical ~ normal(1, .1); // no scaling done for variance

  // Correlation matrix and scale
  // Prior implies marginal correlations weakly regularized away from extremes
  L_Omega ~ lkj_corr_cholesky(2);
  tau ~ normal(.5, .2);
  // nuisance param, iid normal for ncp of Omega (& therefore MVN var-covar)
  to_vector(z_omega) ~ std_normal();
} 

generated quantities {

  //////////////////////////////
  // Set up generated quantities
  //////////////////////////////

  // Fitted mean probabilities 
  matrix[N+num_days_nowcast, ncat] p_hat[K];
  
  // Posterior predictive distribution, for posterior predictive checking
  int Y_tilde[K, N, ncat];
  
  // Log likelihood, for model diagnostics
  //matrix[K, N] log_lik;
  // Correlation matrix
  matrix[ncomp, ncomp] Omega_hat;
  // Recover mu_r and r on original scale from standardized predictors
  matrix[ncomp, K] r_hat;
  vector[ncomp] mu_hat;
  
  
  ///////////////////////////
  // Generate quantities
  ///////////////////////////

  // Convert Cholesky factor corr to full corr matrix
  Omega_hat = multiply_lower_tri_self_transpose(L_Omega);
  
  // https://stats.stackexchange.com/questions/146754/recovering-original-regression-coefficients-from-standardized
  r_hat = r_scaled / sd_t;
  mu_hat = mu_r / sd_t; 
    
  for(k in 1:K){
    
    // Initialize country-specific vectors from matrices of all params
    vector[ncomp] b0_k;
    vector[ncomp] r_k;
    // Initialize linear predictor term
    vector[ncat] eta;
    
    b0_k = to_vector(b0[,k]);
    r_k = to_vector(r_scaled[,k]);

    for(i in 1:(N+num_days_nowcast)){

      // Base case
      eta[1] = 0;
  
      eta[2:ncat] =  b0_k + (r_k * t_scale[i]);
      
      p_hat[k, i] = to_row_vector(softmax(eta));
      
      // Only generate log lik and predicted count if seqs observed on that day
      if(i <= N){
         
        if(trials[k, i] > 0){ 
            
            Y_tilde[k, i] = multinomial_rng(softmax(eta), trials[k, i]);
        
            //log_lik[k, i] = multinomial_lpmf(Y[k, i] | softmax(eta));
        
        }
      }
    }
  }
}

