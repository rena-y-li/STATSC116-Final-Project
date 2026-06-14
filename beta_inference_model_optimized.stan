// beta_inference_model_optimized.stan

// this is the code used to fit both the 15 and 50 gene models, where the horseshoe is replaced by slightly restricted priors
// a sparse matrix formulation was used to increase computational efficiency

data {
  int<lower=0> N;                       // number of cells (20,000 training)
  int<lower=0> P;                       // number of genes (15 or 50)
  vector<lower=0, upper=1>[N] y;        // pseudotime position
  matrix[N, 4] covariates;              // age, sex, PMI, sequencing depth
  
  // making a sparse matrix for the genes so zero-entries are not computed
  int<lower=0> NNZ;                     // total number of non-zero elements
  vector[NNZ] w;                        // the non-zero values
  array[NNZ] int v;                     // column indices
  array[N + 1] int u;                   // row starting indices
  
  // group sizes
  int<lower=0> N_donor;
  int<lower=0> N_sample;
  int<lower=0> N_region;
  int<lower=0> N_broad;

  // grouping indices
  array[N] int donor_id;
  array[N] int sample_id;
  array[N] int region_id;
  array[N_region] int<lower=1, upper=N_broad> region_broad_map; 
}

transformed data {
  // centered covariates to mitigate divergences and speed up MCMC
  matrix[N, 4] covariates_centered;
  for (j in 1:4) {
    covariates_centered[, j] = covariates[, j] - mean(covariates[, j]);
  }
}

parameters {

  // intercept
  real beta_0;   

  // covariate coefficients
  vector[4] beta_cov;   

  // gene coefficients
  vector[P] beta_genes;           
  
  // region hierarchy
  vector[N_broad] alpha_broad;          // broad region intercepts
  vector[N_region] z_region;            // non-centered subregion
  real<lower=0> sigma_region;           // within-broad-region standard deviation
  
  // donor hierarchy
  vector[N_donor] z_donor;              // non-centered donor intercepts
  real<lower=0> sigma_donor;            // donor standard deviation
  
  // sample hierarchy nested in donor
  vector[N_sample] z_sample;            // non-centered sample intercepts
  real<lower=0> sigma_sample;           // sample standard deviation
  
  // precision parameter for Beta
  real<lower=0> phi; 
}

transformed parameters {
  // vectorized hierarchy calculation for efficiency
  vector[N_region] alpha_region = alpha_broad[region_broad_map] + z_region * sigma_region;

  // donor, sample intercepts
  vector[N_donor] alpha_donor = z_donor * sigma_donor;
  vector[N_sample] alpha_sample = z_sample * sigma_sample;
}

model {

  // gene prior (restrictive normal)
  beta_genes   ~ normal(0, 0.5);

  // intercept and covariate priors (weakly informative defaults)
  beta_0       ~ normal(0, 1);
  beta_cov     ~ normal(0, 1);

  // Beta precision parameter prior
  phi          ~ exponential(1);

  // regional priors
  alpha_broad  ~ normal(0, 1);
  z_region     ~ std_normal();
  sigma_region ~ exponential(1);

  // donor priors
  z_donor      ~ std_normal();
  sigma_donor  ~ exponential(1);

  // sample priors
  z_sample     ~ std_normal();
  sigma_sample ~ exponential(1);
  
  // sparse matrix calculation
  vector[N] gene_effects = csr_matrix_times_vector(N, P, w, v, u, beta_genes);
  vector[N] eta = beta_0 + 
                  covariates_centered[,1] * beta_cov[1] + 
                  covariates_centered[,2] * beta_cov[2] + 
                  covariates_centered[,3] * beta_cov[3] + 
                  covariates_centered[,4] * beta_cov[4] + 
                  gene_effects + 
                  alpha_region[region_id] + 
                  alpha_donor[donor_id] + 
                  alpha_sample[sample_id];
  
  // Likelihood evaluation with vectorized inverse-logit
  y ~ beta_proportion(inv_logit(eta), phi);
}

generated quantities {
  vector[N] log_lik;

  // calculat egene effects
  vector[N] gene_effects_gq = csr_matrix_times_vector(N, P, w, v, u, beta_genes);
  
  // initially compute eta_gq (log-odds)
  vector[N] eta_gq = beta_0 + 
                     covariates_centered[,1] * beta_cov[1] + 
                     covariates_centered[,2] * beta_cov[2] + 
                     covariates_centered[,3] * beta_cov[3] + 
                     covariates_centered[,4] * beta_cov[4] + 
                     gene_effects_gq + 
                     alpha_region[region_id] + 
                     alpha_donor[donor_id] + 
                     alpha_sample[sample_id];
                     
  // compute log likelihood with precomputed eta_gq.
  for (i in 1:N) {
    log_lik[i] = beta_proportion_lpdf(y[i] | inv_logit(eta_gq[i]), phi);
  }
}

