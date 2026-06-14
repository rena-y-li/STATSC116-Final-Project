# project_workflow.R. this document provides code used to clean the dataset, as well as run the models and analyze performance

######## DATA CLEANING ########

# load cortex oligodendrocyte lineage pseudotime data
load("cortex_oligos_pseudotime_final.rdata")
library(Seurat)
library(dplyr)
library(rstan)

# constrain pseudotime to between 0 and 1
pdtime_squished <- (cortex_pdtime * (length(cortex_pdtime) - 1) + 0.5) / length(cortex_pdtime)

# find highly variable genes and subset to only those expressed in > 5% of cells
hvgs <- VariableFeatures(sub, assay = "SCT")
counts_matrix <- GetAssayData(sub, assay = "SCT", slot = "counts")
features_5pct <- which((Matrix::rowSums(counts_matrix > 0) / ncol(sub)) > 0.05) |> names()
cortex_genes_final <- intersect(hvgs, features_5pct)
dropped_genes <- setdiff(hvgs, cortex_genes_final)

# known cortex oligodendrocyte marker genes
cortex_markers <- c(
    "ASPA", "BMP4", "CD9", "CD82", "CLDN11", "CNP", "CSPG4", "ENPP4", 
    "GPR17", "MAG", "MAL", "MBP", "MOBP", "MOG", "MYRF", "MYT1", 
    "NCAM1", "NKX2-2", "NKX6-2", "OLIG1", "OLIG2", "OPALIN", 
    "PDGFRA", "PLP1", "PTPRZ1", "SMARCA4", "SOX10", "SOX17", 
    "TF", "ZFP488", "ZNF24", "ZNF536"
)

# add back dropped key marker genes
dropped_cortex_markers <- intersect(cortex_markers, dropped_genes)
cortex_genes_final <- unique(c(cortex_genes_final, dropped_cortex_markers))

deng_sub <- GetAssayData(sub, assay = "SCT", slot = "data")
expr_mat <- t(as.matrix(deng_sub[cortex_genes_final, ]))

# create dataframe
dat <- data.frame(
  pseudotime         = pdtime_squished,
  PatientID          = meta$PatientID,
  sample_id          = meta$sample_id,
  Brain_Region       = meta$Brain_Region,
  Broad_Brain_Region = "Cortex",
  Age_z              = scale(meta$Age),
  Sex                = meta$Sex,
  PMI_z              = scale(meta$PMI),
  nCount_SCT_z       = scale(log(meta$nCount_SCT))
) |> cbind(as.data.frame(expr_mat))

# split data by 80/20 train/test 
set.seed(116)
all_idx <- seq_len(nrow(dat))
train_all <- sample(all_idx, size = floor(0.8 * nrow(dat)))
test_idx <- setdiff(all_idx, train_all)

set.seed(116)
train_idx <- sample(train_all, size = 20000)

# create training and testing datasets
dat_train <- dat[train_idx, ]
dat_test <- dat[test_idx, ]


####### SIMPLE MODELS FOR COMPARISON ########

lm_simple <- lm(pseudotime ~ Brain_Region + Age_z + Sex + PMI_z + nCount_SCT_z, data = dat_train)
summary(lm_simple)

library(lme4)
library(performance)

lmer_simple <- lmer(pseudotime ~ Brain_Region + Age_z + Sex + PMI_z + nCount_SCT_z + (1 | PatientID) + (1 | PatientID:sample_id),
                    data = dat_train, REML = TRUE)
r2(lmer_simple)


######## MODEL 1 - BETA WITH HORSESHOE PRIORS ########

# numerically code sex
dat_train$Sex_int <- ifelse(dat_train$Sex == "M", 1, 0)

# create covariates
region_levels <- levels(factor(dat_train$Brain_Region))
broad_levels <- levels(factor(dat_train$Broad_Brain_Region))

donor_id <- as.integer(factor(dat_train$PatientID))
sample_id <- as.integer(factor(dat_train$sample_id))
region_id <- as.integer(factor(dat_train$Brain_Region, levels = region_levels))
broad_id <- as.integer(factor(dat_train$Broad_Brain_Region, levels = broad_levels))

# map regions to numbers
region_broad_map <- dat_train |>
  distinct(Brain_Region, Broad_Brain_Region) |>
  mutate(
    region_id = as.integer(factor(Brain_Region, levels = region_levels)),
    broad_id  = as.integer(factor(Broad_Brain_Region, levels = broad_levels))
  ) |>
  arrange(region_id) |>
  pull(broad_id)

# sparse matrix for computation
X_dense <- as.matrix(dat_train[, cortex_genes_final])
sparse_parts <- rstan::extract_sparse_parts(X_dense)

# data to use for the .stan model
stan_dat <- list(
  N                = nrow(dat_train),
  P                = length(cortex_genes_final),
  y                = dat_train$pseudotime,
  covariates       = cbind(dat_train$Age_z, dat_train$Sex_int, dat_train$PMI_z, dat_train$nCount_SCT_z),
  donor_id         = donor_id,
  sample_id        = sample_id,
  region_id        = region_id,
  broad_id         = broad_id,
  region_broad_map = region_broad_map,
  N_donor          = length(unique(donor_id)),
  N_sample         = length(unique(sample_id)),
  N_region         = length(unique(region_id)),
  N_broad          = length(unique(broad_id)),
  w                = sparse_parts$w,
  v                = sparse_parts$v,
  u                = sparse_parts$u,
  NNZ              = length(sparse_parts$w)
)
rm(X_dense, sparse_parts); gc()
stan_dat$scale_global <- 0.0015

# run model with 4 chains
options(mc.cores = 4)
rstan_options(auto_write = TRUE)

beta_fit_sparse <- stan(
  file    = "beta_model_sparse_horseshoe.stan",
  data    = stan_dat,
  chains  = 4,
  cores   = 4,
  iter    = 2000,
  warmup  = 1000,
  control = list(adapt_delta = 0.95, max_treedepth = 11),
  pars    = c("beta_0", "beta_cov", "beta_genes", "phi", "alpha_region", "alpha_broad")
)
saveRDS(beta_fit_sparse, file = "cortex_horseshoe_trajectory_beta_fit_sparse.rds")


######## SECONDARY MODEL DATA PREPARATION ########

library(Matrix)

donor_fac <- as.factor(dat_train$PatientID)
sample_fac <- as.factor(dat_train$sample_id)
region_fac <- as.factor(dat_train$Brain_Region)
broad_fac <- as.factor(rep("Cortex", nrow(dat_train)))

covariates_mat_20k <- cbind(
  dat_train$Age_z,
  as.numeric(as.factor(dat_train$Sex)) - 1,
  dat_train$PMI_z,
  dat_train$nCount_SCT_z
)

region_broad_map_20k <- data.frame(
  region = unique(dat_train$Brain_Region),
  broad  = "Cortex"
) |> arrange(region) |>
  mutate(broad_id = as.integer(as.factor(broad))) |>
  pull(broad_id)

make_stan_sparse_data_20k <- function(gene_list) {
  X <- as.matrix(dat_train[, gene_list])
  mat_csr <- as(X, "RsparseMatrix")
  list(
    P   = length(gene_list),
    NNZ = length(mat_csr@x),
    w   = mat_csr@x,
    v   = mat_csr@j + 1,
    u   = mat_csr@p + 1
  )
}

base_stan_data_20k <- list(
  N                = nrow(dat_train),
  y                = dat_train$pseudotime,
  covariates       = covariates_mat_20k,
  N_donor          = length(levels(donor_fac)),
  N_sample         = length(levels(sample_fac)),
  N_region         = length(levels(region_fac)),
  N_broad          = length(levels(broad_fac)),
  donor_id         = as.integer(donor_fac),
  sample_id        = as.integer(sample_fac),
  region_id        = as.integer(region_fac),
  region_broad_map = region_broad_map_20k
)

# creating objects for 15 and 50 gene runs
sparse_strict_15_20k <- make_stan_sparse_data_20k(genes_strict_15)
sparse_moderate_50_20k <- make_stan_sparse_data_20k(genes_moderate_50)

stan_data_strict_20k <- c(base_stan_data_20k, sparse_strict_15_20k)
stan_data_moderate_20k <- c(base_stan_data_20k, sparse_moderate_50_20k)

save(
  stan_data_strict_20k, stan_data_moderate_20k,
  dat, train_all, test_idx,
  genes_strict_15, genes_moderate_50,
  file = "cortex_inference_20k_standata.rdata"
)


######## MODEL 2 - BETA WITH 15 GENES ########

library(rstan)
library(dplyr)

load("cortex_inference_20k_standata.rdata")

options(mc.cores = 4)
rstan_options(auto_write = TRUE)

compiled_model <- stan_model("beta_inference_model_optimized.stan")

cat("Starting 20k strict model (15 genes)...\n")
fit_strict_20k <- sampling(
  compiled_model,
  data    = stan_data_strict_20k,
  chains  = 4,
  cores   = 4,
  iter    = 2000,
  warmup  = 1000,
  seed    = 116,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

saveRDS(fit_strict_20k, file = "cortex_beta_strict_15_20k_fit.rds")

sampler_params <- get_sampler_params(fit_strict_20k, inc_warmup = FALSE)
num_divergences <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
num_treedepths <- sum(sapply(sampler_params, function(x) sum(x[, "treedepth__"] >= 15)))

summary_strict_20k <- summary(fit_strict_20k, probs = c(0.1, 0.5, 0.9))$summary
gene_rows_strict <- grep("^beta_genes\\[", rownames(summary_strict_20k))
summary_genes_strict <- summary_strict_20k[gene_rows_strict, , drop = FALSE]

df_strict_20k <- data.frame(
  Gene      = genes_strict_15,
  Mean      = summary_genes_strict[, "mean"],
  SD        = summary_genes_strict[, "sd"],
  Lower_80  = summary_genes_strict[, "10%"],
  Median_50 = summary_genes_strict[, "50%"],
  Upper_80  = summary_genes_strict[, "90%"],
  Rhat      = summary_genes_strict[, "Rhat"],
  n_eff     = summary_genes_strict[, "n_eff"]
)

write.csv(df_strict_20k, file = "cortex_beta_summary_strict_15_20k_optimized_.99_15.csv", row.names = FALSE)
cat("Done.\n")


######## MODEL 2 - BETA WITH 50 GENES ########

library(rstan)
library(dplyr)

load("cortex_inference_20k_standata.rdata")

options(mc.cores = 4)
rstan_options(auto_write = TRUE)

compiled_model <- stan_model("beta_inference_model_optimized.stan")

cat("Starting 20k moderate model (50 genes)...\n")
fit_moderate_20k <- sampling(
  compiled_model,
  data    = stan_data_moderate_20k,
  chains  = 4,
  cores   = 4,
  iter    = 2000,
  warmup  = 1000,
  seed    = 116,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

saveRDS(fit_moderate_20k, file = "cortex_beta_moderate_50_20k_fit.rds")

sampler_params <- get_sampler_params(fit_moderate_20k, inc_warmup = FALSE)
num_divergences <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
num_treedepths <- sum(sapply(sampler_params, function(x) sum(x[, "treedepth__"] >= 15)))

summary_moderate_20k <- summary(fit_moderate_20k, probs = c(0.1, 0.5, 0.9))$summary
gene_rows_moderate <- grep("^beta_genes\\[", rownames(summary_moderate_20k))
summary_genes_moderate <- summary_moderate_20k[gene_rows_moderate, , drop = FALSE]

df_moderate_20k <- data.frame(
  Gene      = genes_moderate_50,
  Mean      = summary_genes_moderate[, "mean"],
  SD        = summary_genes_moderate[, "sd"],
  Lower_80  = summary_genes_moderate[, "10%"],
  Median_50 = summary_genes_moderate[, "50%"],
  Upper_80  = summary_genes_moderate[, "90%"],
  Rhat      = summary_genes_moderate[, "Rhat"],
  n_eff     = summary_genes_moderate[, "n_eff"]
)

write.csv(df_moderate_20k, file = "cortex_beta_summary_moderate_50_20k_optimized_.99_15.csv", row.names = FALSE)
cat("Done.\n")


######## LOO MODEL COMPARISON ########

# fit objects
fit_strict_20k <- readRDS("cortex_beta_strict_15_20k_fit.rds")
fit_moderate_20k <- readRDS("cortex_beta_moderate_50_20k_fit.rds")

# LOO
library(loo)

log_lik_strict <- extract_log_lik(fit_strict_20k,   parameter_name = "log_lik", merge_chains = FALSE)
log_lik_moderate <- extract_log_lik(fit_moderate_20k, parameter_name = "log_lik", merge_chains = FALSE)

loo_strict <- loo(log_lik_strict,   r_eff = relative_eff(exp(log_lik_strict)))
loo_moderate <- loo(log_lik_moderate, r_eff = relative_eff(exp(log_lik_moderate)))

loo_compare(loo_strict, loo_moderate)


######## 50 GENE PPC ########

library(bayesplot)

y_rep_moderate <- as.matrix(extract(fit_moderate_20k)$y_rep)
y_train <- dat[train_idx, "pseudotime"]

png("ppc_density_moderate_50.png", width = 800, height = 600)
print(ppc_dens_overlay(y_train, y_rep_moderate[1:100, ]))
dev.off()

png("ppc_stat_sd_moderate.png", width = 800, height = 600)
print(ppc_stat(y_train, y_rep_moderate, stat = "sd"))
dev.off()

stat_upper <- function(y) mean(y > 0.9)
stat_lower <- function(y) mean(y < 0.3)

png("ppc_stat_upper_moderate.png", width = 800, height = 600)
print(ppc_stat(y_train, y_rep_moderate, stat = stat_upper))
dev.off()

png("ppc_stat_lower_moderate.png", width = 800, height = 600)
print(ppc_stat(y_train, y_rep_moderate, stat = stat_lower))
dev.off()


######## TEST SET PERFORMANCE and RMSE ########
                              
posterior_strict <- extract(fit_strict_20k)
posterior_moderate <- extract(fit_moderate_20k)
n_draws <- length(posterior_strict$beta_0)

mu_test_strict <- matrix(0, nrow = n_draws, ncol = nrow(dat_test))
for (i in 1:n_draws) {
  eta <- posterior_strict$beta_0[i] +
         covariates_test_centered %*% posterior_strict$beta_cov[i,] +
         X_test_strict %*% posterior_strict$beta_genes[i,] +
         posterior_strict$alpha_region[i, region_id_test] +
         posterior_strict$alpha_donor[i, donor_id_test] +
         posterior_strict$alpha_sample[i, sample_id_test]
  mu_test_strict[i,] <- plogis(eta)
}

mu_test_moderate <- matrix(0, nrow = n_draws, ncol = nrow(dat_test))
for (i in 1:n_draws) {
  eta <- posterior_moderate$beta_0[i] +
         covariates_test_centered %*% posterior_moderate$beta_cov[i,] +
         X_test_moderate %*% posterior_moderate$beta_genes[i,] +
         posterior_moderate$alpha_region[i, region_id_test] +
         posterior_moderate$alpha_donor[i, donor_id_test] +
         posterior_moderate$alpha_sample[i, sample_id_test]
  mu_test_moderate[i,] <- plogis(eta)
}

mu_pred_strict <- colMeans(mu_test_strict)
mu_pred_moderate <- colMeans(mu_test_moderate)
y_test <- dat_test$pseudotime

rmse_strict <- sqrt(mean((mu_pred_strict - y_test)^2))
rmse_moderate <- sqrt(mean((mu_pred_moderate - y_test)^2))

cat("RMSE strict (15 genes):", rmse_strict, "\n")
cat("RMSE moderate (50 genes):", rmse_moderate, "\n")


######## CONVERGENCE DIAGNOSTICS ########
                             
summary(fit_strict_20k)$summary |> as.data.frame() |> 
  summarise(max_rhat = max(Rhat, na.rm=TRUE), min_ess = min(n_eff, na.rm=TRUE))

summary(fit_moderate_20k)$summary |> as.data.frame() |> 
  summarise(max_rhat = max(Rhat, na.rm=TRUE), min_ess = min(n_eff, na.rm=TRUE))

summary(beta_fit_sparse)$summary |> as.data.frame() |> 
  summarise(max_rhat = max(Rhat, na.rm=TRUE), min_ess = min(n_eff, na.rm=TRUE))

