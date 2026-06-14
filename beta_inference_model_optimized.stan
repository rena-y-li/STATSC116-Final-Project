// beta_inference_model_optimized.stan

// this was the code used to fit both the 15 and 50 gene models

data {
  int<lower=0> N;                       // 20,000 cells
  int<lower=0> P;                       // number of genes (15 or 50)
  vector<lower=0, upper=1>[N] y;        // pseudotime outcome
  matrix[N, 4] covariates;              // raw age, sex, pmi, depth
  
  // Sparse Matrix Components for Selected Genes
  int<lower=0> NNZ;                     // total number of non-zero elements
  vector[NNZ] w;                        // the actual non-zero expression values
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
  // Center the covariates exactly once before sampling starts
  matrix[N, 4] covariates_centered;
  for (j in 1:4) {
    covariates_centered[, j] = covariates[, j] - mean(covariates[, j]);
  }
}

parameters {
  real beta_0;                          // global intercept
  vector[4] beta_cov;                   // age, sex, pmi, depth coefficients
  vector[P] beta_genes;                 // coefficients for selected genes
  
  // region hierarchy
  vector[N_broad] alpha_broad;          // broad region intercepts
  vector[N_region] z_region;            // non-centered fine region
  real<lower=0> sigma_region;           // within-broad SD
  
  // donor hierarchy
  vector[N_donor] z_donor;              // non-centered donor intercepts
  real<lower=0> sigma_donor;            // donor SD
  
  // sample hierarchy nested in donor
  vector[N_sample] z_sample;            // non-centered sample intercepts
  real<lower=0> sigma_sample;           // sample SD
  
  // Precision parameter for Beta distribution
  real<lower=0> phi; 
}

transformed parameters {
  // Fully vectorized hierarchical lookups (Instantaneous)
  vector[N_region] alpha_region = alpha_broad[region_broad_map] + z_region * sigma_region;
  vector[N_donor] alpha_donor = z_donor * sigma_donor;
  vector[N_sample] alpha_sample = z_sample * sigma_sample;
}

model {
  // Priors
  beta_genes   ~ normal(0, 0.5);
  beta_0       ~ normal(0, 1);
  beta_cov     ~ normal(0, 1);
  phi          ~ exponential(1);
  
  alpha_broad  ~ normal(0, 1);
  z_region     ~ std_normal();
  sigma_region ~ exponential(1);
  
  z_donor      ~ std_normal();
  sigma_donor  ~ exponential(1);
  
  z_sample     ~ std_normal();
  sigma_sample ~ exponential(1);
  
  // Sparse Matrix Calculation
  vector[N] gene_effects = csr_matrix_times_vector(N, P, w, v, u, beta_genes);
  
  // Fully vectorized and unrolled linear predictor (Bypasses heavy matrix multiplication graph)
  vector[N] eta = beta_0 + 
                  covariates_centered[,1] * beta_cov[1] + 
                  covariates_centered[,2] * beta_cov[2] + 
                  covariates_centered[,3] * beta_cov[3] + 
                  covariates_centered[,4] * beta_cov[4] + 
                  gene_effects + 
                  alpha_region[region_id] + 
                  alpha_donor[donor_id] + 
                  alpha_sample[sample_id];
  
  // Likelihood evaluation using vectorized inverse-logit
  y ~ beta_proportion(inv_logit(eta), phi);
}

generated quantities {
  vector[N] log_lik;

  // 1. Calculate gene effects exactly once for this draw
  vector[N] gene_effects_gq = csr_matrix_times_vector(N, P, w, v, u, beta_genes);
  
  // 2. Reconstruct eta_gq cleanly as a single vectorized operation
  vector[N] eta_gq = beta_0 + 
                     covariates_centered[,1] * beta_cov[1] + 
                     covariates_centered[,2] * beta_cov[2] + 
                     covariates_centered[,3] * beta_cov[3] + 
                     covariates_centered[,4] * beta_cov[4] + 
                     gene_effects_gq + 
                     alpha_region[region_id] + 
                     alpha_donor[donor_id] + 
                     alpha_sample[sample_id];
                     
  // 3. Compute log_lik point-by-point.
  // Because eta_gq is pre-computed, this loop is vastly faster than before.
  for (i in 1:N) {
    log_lik[i] = beta_proportion_lpdf(y[i] | inv_logit(eta_gq[i]), phi);
  }
}

