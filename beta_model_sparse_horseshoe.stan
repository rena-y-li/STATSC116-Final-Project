// beta_model_sparse_horseshoe.stan 

// this is the initial hierarchical beta model with the full 358 gene set and horseshoe priors
// a sparse matrix formulation was used to increase computational efficiency

data {
  int<lower=0> N;                      // number of cells (20,000 training)
  int<lower=0> P;                      // number of genes (358)
  vector<lower=0, upper=1>[N] y;       // pseudotime position
  matrix[N, 4] covariates;             // age, sex, PMI, sequencing depth
  
  // making a sparse matrix for the genes so zero-entries are not computed
  int<lower=0> NNZ;                    // total number of non-zero elements
  vector[NNZ] w;                       // the non-zero values
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
  array[N_region] int<lower=1, upper=N_broad> region_broad_map;      // number of regions
  
  // horseshoe hyperparameter
  real<lower=0> scale_global;
}

parameters {
  // intercept
  real beta_0;
  
  // covariate coefficients
  vector[4] beta_cov;                  // hardcoded for 4 covariates
  
  // horseshoe for genes
  vector[P] z;                         // non-centered gene coefficients
  vector<lower=0>[P] lambda;           // local scales
  real<lower=1e-5> tau;                // global scale with lower bound
  real<lower=0> caux_global;           // slab global auxiliary for horseshoe
  
  // region hierarchal structure
  vector[N_broad] alpha_broad;         // broad region intercepts
  vector[N_region] z_region;           // non-centered subregion
  real<lower=0> sigma_region;          // within-broad-region standard deviation
  
  // donor hierarchy
  vector[N_donor] z_donor;             // non-centered donor intercepts
  real<lower=0> sigma_donor;           // donor standard deviation
  
  // sample hierarchy nested in donor
  vector[N_sample] z_sample;           // non-centered sample intercepts
  real<lower=0> sigma_sample;          // sample standard deviation
  
  // precision parameter for Beta
  real<lower=0> phi; 
}

transformed parameters {
  // horseshoe
  real slab_scale = 2.0;
  real slab_df    = 4.0;
  real c2 = square(slab_scale) * caux_global; 
  
  vector[P] lambda_tilde;
  vector[P] beta_genes;
  
  // vectorized horseshoe calculation for efficiency
  vector[P] lam2 = square(lambda);
  lambda_tilde = sqrt( (c2 * lam2) ./ (c2 + square(tau) * lam2) );
  beta_genes = z .* lambda_tilde * tau;
  
  // subregion intercepts, drawn from broad region mean (modeled as random but fixed would have been better)
  vector[N_region] alpha_region;
  for (r in 1:N_region) {
    alpha_region[r] = alpha_broad[region_broad_map[r]] + z_region[r] * sigma_region;
  }
  
  // donor, sample intercepts
  vector[N_donor] alpha_donor = z_donor * sigma_donor;
  vector[N_sample] alpha_sample = z_sample * sigma_sample;
}

model {
  // horseshoe priors 
  z           ~ std_normal();
  lambda      ~ student_t(1, 0, 1);
  tau         ~ student_t(1, 0, scale_global);
  caux_global ~ inv_gamma(2.0, 2.0);   
  
  // intercept and covariate priors (weakly informative defaults)
  beta_0   ~ normal(0, 1);
  beta_cov ~ normal(0, 1);
  
  // Beta precision parameter prior
  phi ~ exponential(1);
  
  // regional priors
  alpha_broad  ~ normal(0, 1);
  z_region     ~ std_normal();
  sigma_region ~ exponential(1);
  
  // donor priors
  z_donor     ~ std_normal();
  sigma_donor ~ exponential(1);
  
  // sample priors
  z_sample     ~ std_normal();
  sigma_sample ~ exponential(1);
  
  // likelihood mapping everything to (0,1) with inv_logit
  vector[N] mu;
  
  // sparse matrix calculation
  vector[N] gene_effects = csr_matrix_times_vector(N, P, w, v, u, beta_genes);
  
  mu = inv_logit(beta_0 + 
                 covariates * beta_cov + 
                 gene_effects + 
                 alpha_region[region_id] + 
                 alpha_donor[donor_id] + 
                 alpha_sample[sample_id]);
  
  y ~ beta_proportion(mu, phi);
}

