# PURPOSE: Generate the AdaPT CV results for Positional

# WARNING: This code was run on a server

# Access necessary packages
library(tidyverse)
library(adaptMT)
library(parallel)

# ------------------------------------------------------------------------------

# r^2 = .25 ------------------------------------------------------------

# Load the model dataset
model_data <-
  read_csv("data/adapt_model_data/positional/rsquared25_model_data.csv")

# Get the variables:

# GTEx slopes
gtex_eqtl_slope_cols <- stringr::str_subset(
  colnames(model_data), "slope") %>%
  # Only use the subset that are the averages:
  stringr::str_subset("ave")

# BrainVar slopes
brainvar_eqtl_slope_cols <- stringr::str_subset(
  colnames(model_data), "beta") %>%
  # Only use the subset that are the averages:
  stringr::str_subset("ave")

# eqtl indicators:
eqtl_ind_cols <- stringr::str_subset(colnames(model_data),
                                     "_eqtl$")

# biotype indicators:
biotype_cols <- stringr::str_subset(colnames(model_data),
                                    "biotype_")

# LOEUF - use the minimum
loeuf_cols <- c("min_loeuf")

# Size related variables:
size_cols <- c("n_snps")

# Now each of the wgcna_label columns:
wgcna_cols <- stringr::str_subset(colnames(model_data),
                                  "any_gene_")

# Other phenotypes:
asd_pheno_cols <- c("cap_ea_quad_z", "cap_scz_quad_z")
scz_pheno_cols <- c("cap_ea_quad_z", "cap_asd_quad_z")
ea_pheno_cols <- c("cap_scz_quad_z", "cap_asd_quad_z")

# Make a vector with all of these variables:
asd_variables <- c(asd_pheno_cols, gtex_eqtl_slope_cols, brainvar_eqtl_slope_cols,
                   eqtl_ind_cols, size_cols, loeuf_cols, biotype_cols, wgcna_cols)
scz_variables <- c(scz_pheno_cols, gtex_eqtl_slope_cols, brainvar_eqtl_slope_cols,
                   eqtl_ind_cols, size_cols, loeuf_cols, biotype_cols, wgcna_cols)
ea_variables <- c(ea_pheno_cols, gtex_eqtl_slope_cols, brainvar_eqtl_slope_cols,
                  eqtl_ind_cols, size_cols, loeuf_cols, biotype_cols, wgcna_cols)

model_args <- list("nrounds450md1eta06gamma0" = list("nrounds" = 450,
                                                     "max_depth" = 1,
                                                     "eta" = 0.06,
                                                     "gamma" = 0,
                                                     "verbose" = 0,
                                                     "nthread" = 5),
                   "nrounds100md2eta04gamma0" = list("nrounds" = 100,
                                                     "max_depth" = 2,
                                                     "eta" = 0.04,
                                                     "gamma" = 0,
                                                     "verbose" = 0,
                                                     "nthread" = 5))

# ASD:
asd_adapt_results <- adapt_xgboost_cv(as.matrix(
  model_data[,asd_variables]),
  model_data$asd_quad_pval,
  verbose = list(print = FALSE,
                 fit = FALSE,
                 ms = FALSE),
  piargs = model_args,
  muargs = model_args,
  s0 = rep(0.05, nrow(model_data)),
  n_folds = 5,
  niter_ms = 10,
  nms = as.integer(2))
print(length(which(asd_adapt_results$qvals <= 0.05)))

saveRDS(asd_adapt_results,
        "data/adapt_results/positional/asd/rsquared25_xgb_cv.rds")

# SCZ:
scz_adapt_results <- adapt_xgboost_cv(as.matrix(
  model_data[,scz_variables]),
  model_data$scz_quad_pval,
  verbose = list(print = FALSE,
                 fit = FALSE,
                 ms = FALSE),
  piargs = model_args,
  muargs = model_args,
  s0 = rep(0.05, nrow(model_data)),
  n_folds = 5,
  niter_ms = 10,
  nms = as.integer(2))
print(length(which(scz_adapt_results$qvals <= 0.05)))

saveRDS(scz_adapt_results,
        "data/adapt_results/positional/scz/rsquared25_xgb_cv.rds")

# EA:
ea_adapt_results <- adapt_xgboost_cv(as.matrix(
  model_data[,ea_variables]),
  model_data$ea_quad_pval,
  verbose = list(print = FALSE,
                 fit = FALSE,
                 ms = FALSE),
  piargs = model_args,
  muargs = model_args,
  s0 = rep(0.05, nrow(model_data)),
  n_folds = 5,
  niter_ms = 10,
  nms = as.integer(2))
print(length(which(ea_adapt_results$qvals <= 0.05)))

saveRDS(ea_adapt_results,
        "data/adapt_results/positional/ea/rsquared25_xgb_cv.rds")

# r^2 = .50 ------------------------------------------------------------

# Load the model dataset
model_data <-
  read_csv("data/adapt_model_data/positional/rsquared50_model_data.csv")

# Get the variables:

# GTEx slopes
gtex_eqtl_slope_cols <- stringr::str_subset(
  colnames(model_data), "slope") %>%
  # Only use the subset that are the averages:
  stringr::str_subset("ave")

# BrainVar slopes
brainvar_eqtl_slope_cols <- stringr::str_subset(
  colnames(model_data), "beta") %>%
  # Only use the subset that are the averages:
  stringr::str_subset("ave")

# eqtl indicators:
eqtl_ind_cols <- stringr::str_subset(colnames(model_data),
                                     "_eqtl$")

# biotype indicators:
biotype_cols <- stringr::str_subset(colnames(model_data),
                                    "biotype_")

# LOEUF - use the minimum
loeuf_cols <- c("min_loeuf")

# Size related variables:
size_cols <- c("n_snps")

# Now each of the wgcna_label columns:
wgcna_cols <- stringr::str_subset(colnames(model_data),
                                  "any_gene_")

# Other phenotypes:
asd_pheno_cols <- c("cap_ea_quad_z", "cap_scz_quad_z")
scz_pheno_cols <- c("cap_ea_quad_z", "cap_asd_quad_z")
ea_pheno_cols <- c("cap_scz_quad_z", "cap_asd_quad_z")

# Make a vector with all of these variables:
asd_variables <- c(asd_pheno_cols, gtex_eqtl_slope_cols, brainvar_eqtl_slope_cols,
                   eqtl_ind_cols, size_cols, loeuf_cols, biotype_cols, wgcna_cols)
scz_variables <- c(scz_pheno_cols, gtex_eqtl_slope_cols, brainvar_eqtl_slope_cols,
                   eqtl_ind_cols, size_cols, loeuf_cols, biotype_cols, wgcna_cols)
ea_variables <- c(ea_pheno_cols, gtex_eqtl_slope_cols, brainvar_eqtl_slope_cols,
                  eqtl_ind_cols, size_cols, loeuf_cols, biotype_cols, wgcna_cols)

model_args <- list("nrounds150md1eta06gamma0" = list("nrounds" = 150,
                                                     "max_depth" = 1,
                                                     "eta" = 0.06,
                                                     "gamma" = 0,
                                                     "verbose" = 0,
                                                     "nthread" = 5),
                   "nrounds250md2eta06gamma0" = list("nrounds" = 250,
                                                     "max_depth" = 2,
                                                     "eta" = 0.06,
                                                     "gamma" = 0,
                                                     "verbose" = 0,
                                                     "nthread" = 5))

# ASD:
asd_adapt_results <- adapt_xgboost_cv(as.matrix(
  model_data[,asd_variables]),
  model_data$asd_quad_pval,
  verbose = list(print = FALSE,
                 fit = FALSE,
                 ms = FALSE),
  piargs = model_args,
  muargs = model_args,
  s0 = rep(0.05, nrow(model_data)),
  n_folds = 5,
  niter_ms = 10,
  nms = as.integer(2))
print(length(which(asd_adapt_results$qvals <= 0.05)))

saveRDS(asd_adapt_results,
        "data/adapt_results/positional/asd/rsquared50_xgb_cv.rds")

# SCZ:
scz_adapt_results <- adapt_xgboost_cv(as.matrix(
  model_data[,scz_variables]),
  model_data$scz_quad_pval,
  verbose = list(print = FALSE,
                 fit = FALSE,
                 ms = FALSE),
  piargs = model_args,
  muargs = model_args,
  s0 = rep(0.05, nrow(model_data)),
  n_folds = 5,
  niter_ms = 10,
  nms = as.integer(2))
print(length(which(scz_adapt_results$qvals <= 0.05)))

saveRDS(scz_adapt_results,
        "data/adapt_results/positional/scz/rsquared50_xgb_cv.rds")

# EA:
ea_adapt_results <- adapt_xgboost_cv(as.matrix(
  model_data[,ea_variables]),
  model_data$ea_quad_pval,
  verbose = list(print = FALSE,
                 fit = FALSE,
                 ms = FALSE),
  piargs = model_args,
  muargs = model_args,
  s0 = rep(0.05, nrow(model_data)),
  n_folds = 5,
  niter_ms = 10,
  nms = as.integer(2))
print(length(which(ea_adapt_results$qvals <= 0.05)))

saveRDS(ea_adapt_results,
        "data/adapt_results/positional/ea/rsquared50_xgb_cv.rds")

# r^2 = .75 ------------------------------------------------------------

# Load the model dataset
model_data <-
  read_csv("data/adapt_model_data/positional/rsquared75_model_data.csv")

# Get the variables:

# GTEx slopes
gtex_eqtl_slope_cols <- stringr::str_subset(
  colnames(model_data), "slope") %>%
  # Only use the subset that are the averages:
  stringr::str_subset("ave")

# BrainVar slopes
brainvar_eqtl_slope_cols <- stringr::str_subset(
  colnames(model_data), "beta") %>%
  # Only use the subset that are the averages:
  stringr::str_subset("ave")

# eqtl indicators:
eqtl_ind_cols <- stringr::str_subset(colnames(model_data),
                                     "_eqtl$")

# biotype indicators:
biotype_cols <- stringr::str_subset(colnames(model_data),
                                    "biotype_")

# LOEUF - use the minimum
loeuf_cols <- c("min_loeuf")

# Size related variables:
size_cols <- c("n_snps")

# Now each of the wgcna_label columns:
wgcna_cols <- stringr::str_subset(colnames(model_data),
                                  "any_gene_")

# Other phenotypes:
asd_pheno_cols <- c("cap_ea_quad_z", "cap_scz_quad_z")
scz_pheno_cols <- c("cap_ea_quad_z", "cap_asd_quad_z")
ea_pheno_cols <- c("cap_scz_quad_z", "cap_asd_quad_z")

# Make a vector with all of these variables:
asd_variables <- c(asd_pheno_cols, gtex_eqtl_slope_cols, brainvar_eqtl_slope_cols,
                   eqtl_ind_cols, size_cols, loeuf_cols, biotype_cols, wgcna_cols)
scz_variables <- c(scz_pheno_cols, gtex_eqtl_slope_cols, brainvar_eqtl_slope_cols,
                   eqtl_ind_cols, size_cols, loeuf_cols, biotype_cols, wgcna_cols)
ea_variables <- c(ea_pheno_cols, gtex_eqtl_slope_cols, brainvar_eqtl_slope_cols,
                  eqtl_ind_cols, size_cols, loeuf_cols, biotype_cols, wgcna_cols)

model_args <- list("nrounds450md1eta05gamma0" = list("nrounds" = 450,
                                                     "max_depth" = 1,
                                                     "eta" = 0.05,
                                                     "gamma" = 0,
                                                     "verbose" = 0,
                                                     "nthread" = 5),
                   "nrounds200md2eta05gamma0" = list("nrounds" = 200,
                                                     "max_depth" = 2,
                                                     "eta" = 0.05,
                                                     "gamma" = 0,
                                                     "verbose" = 0,
                                                     "nthread" = 5))

# ASD:
asd_adapt_results <- adapt_xgboost_cv(as.matrix(
  model_data[,asd_variables]),
  model_data$asd_quad_pval,
  verbose = list(print = FALSE,
                 fit = FALSE,
                 ms = FALSE),
  piargs = model_args,
  muargs = model_args,
  s0 = rep(0.05, nrow(model_data)),
  n_folds = 5,
  niter_ms = 10,
  nms = as.integer(2))
print(length(which(asd_adapt_results$qvals <= 0.05)))

saveRDS(asd_adapt_results,
        "data/adapt_results/positional/asd/rsquared75_xgb_cv.rds")

# SCZ:
scz_adapt_results <- adapt_xgboost_cv(as.matrix(
  model_data[,scz_variables]),
  model_data$scz_quad_pval,
  verbose = list(print = FALSE,
                 fit = FALSE,
                 ms = FALSE),
  piargs = model_args,
  muargs = model_args,
  s0 = rep(0.05, nrow(model_data)),
  n_folds = 5,
  niter_ms = 10,
  nms = as.integer(2))
print(length(which(scz_adapt_results$qvals <= 0.05)))

saveRDS(scz_adapt_results,
        "data/adapt_results/positional/scz/rsquared75_xgb_cv.rds")

# EA:
ea_adapt_results <- adapt_xgboost_cv(as.matrix(
  model_data[,ea_variables]),
  model_data$ea_quad_pval,
  verbose = list(print = FALSE,
                 fit = FALSE,
                 ms = FALSE),
  piargs = model_args,
  muargs = model_args,
  s0 = rep(0.05, nrow(model_data)),
  n_folds = 5,
  niter_ms = 10,
  nms = as.integer(2))
print(length(which(ea_adapt_results$qvals <= 0.05)))

saveRDS(ea_adapt_results,
        "data/adapt_results/positional/ea/rsquared75_xgb_cv.rds")
