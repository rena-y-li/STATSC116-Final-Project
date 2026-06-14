// beta_model_sparse.stan (Sparse Matrix Version)

data {
  int<lower=0> N;                      // number of cells
  int<lower=0> P;                      // number of genes
  vector<lower=0, upper=1>[N] y;       // pseudotime outcome
  matrix[N, 4] covariates;             // age, sex, pmi, depth
  
  // Sparse Matrix Components for Genes (Replaces matrix[N, P] X)
  int<lower=0> NNZ;                    // total number of non-zero elements
  vector[NNZ] w;                       // the actual non-zero values
  array[NNZ] int v;                    // column indices
  array[N + 1] int u;                  // row starting indices
  
  
  // group sizes
  int<lower=0> N_donor;
  int<lower=0> N_sample;
  int<lower=0> N_region;
  int<lower=0> N_broad;

  // grouping indices
  array[N] int donor_id;
  array[N] int sample_id;
  array[N] int region_id;
  array[N_region] int<lower=1, upper=N_broad> region_broad_map;      // Number of regions
  
  // horseshoe hyperparameter
  real<lower=0> scale_global;
}

parameters {
  // intercept
  real beta_0;
  
  // covariate coefficients
  vector[4] beta_cov;                  // Hardcoded to 4
  
  // horseshoe for genes
  vector[P] z;                         // non-centered gene coefficients
  vector<lower=0>[P] lambda;           // local scales
  real<lower=1e-5> tau;                // global scale with safety floor
  real<lower=0> caux_global;           // slab global auxiliary
  
  // region hierarchy
  vector[N_broad] alpha_broad;         // broad region intercepts
  vector[N_region] z_region;           // non-centered fine region
  real<lower=0> sigma_region;          // within-broad SD
  
  // donor hierarchy
  vector[N_donor] z_donor;             // non-centered donor intercepts
  real<lower=0> sigma_donor;           // donor SD
  
  // sample hierarchy nested in donor
  vector[N_sample] z_sample;           // non-centered sample intercepts
  real<lower=0> sigma_sample;          // sample SD
  
  // Precision parameter for Beta distribution
  real<lower=0> phi; 
}

transformed parameters {
  // 1. Corrected Regularized Horseshoe (Finnish Horseshoe)
  real slab_scale = 2.0;
  real slab_df    = 4.0;
  real c2 = square(slab_scale) * caux_global; 
  
  vector[P] lambda_tilde;
  vector[P] beta_genes;
  
  // Fully vectorized horseshoe
  vector[P] lam2 = square(lambda);
  lambda_tilde = sqrt( (c2 * lam2) ./ (c2 + square(tau) * lam2) );
  beta_genes = z .* lambda_tilde * tau;
  
  // 2. Fine region intercepts drawn from broad region mean
  vector[N_region] alpha_region;
  for (r in 1:N_region) {
    alpha_region[r] = alpha_broad[region_broad_map[r]] + z_region[r] * sigma_region;
  }
  
  // 3. Donor and Sample intercepts
  vector[N_donor] alpha_donor = z_donor * sigma_donor;
  vector[N_sample] alpha_sample = z_sample * sigma_sample;
}

model {
  // horseshoe priors 
  z           ~ std_normal();
  lambda      ~ student_t(1, 0, 1);
  tau         ~ student_t(1, 0, scale_global);
  caux_global ~ inv_gamma(2.0, 2.0);   
  
  // covariate priors
  beta_0   ~ normal(0, 1);
  beta_cov ~ normal(0, 1);
  
  // precision prior for Beta distribution
  phi ~ exponential(1);
  
  // region hierarchy priors
  alpha_broad  ~ normal(0, 1);
  z_region     ~ std_normal();
  sigma_region ~ exponential(1);
  
  // donor priors
  z_donor     ~ std_normal();
  sigma_donor ~ exponential(1);
  
  // sample priors
  z_sample     ~ std_normal();
  sigma_sample ~ exponential(1);
  
  // likelihood mapping everything to (0,1) via inv_logit
  vector[N] mu;
  
  // SPARSE MATRIX CALCULATION (This replaces X * beta_genes)
  vector[N] gene_effects = csr_matrix_times_vector(N, P, w, v, u, beta_genes);
  
  mu = inv_logit(beta_0 + 
                 covariates * beta_cov + 
                 gene_effects + 
                 alpha_region[region_id] + 
                 alpha_donor[donor_id] + 
                 alpha_sample[sample_id]);
  
  y ~ beta_proportion(mu, phi);
}

