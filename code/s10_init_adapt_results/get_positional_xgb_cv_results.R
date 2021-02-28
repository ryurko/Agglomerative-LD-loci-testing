# PURPOSE: Initialize the AdaPT-CV ASD, SCZ, and EA results for Positional data

library(tidyverse)
library(adaptMT)

# r^2 = 0.25 --------------------------------------------------------------

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
scz_pheno_cols <- c("cap_ea_quad_z", "cap_asd_quad_z")
asd_pheno_cols <- c("cap_ea_quad_z", "cap_scz_quad_z")
ea_pheno_cols <- c("cap_scz_quad_z", "cap_asd_quad_z")

# Make a vector with all of these variables:
asd_variables <- c(asd_pheno_cols, gtex_eqtl_slope_cols, brainvar_eqtl_slope_cols,
                   eqtl_ind_cols, size_cols, loeuf_cols, biotype_cols, wgcna_cols)
scz_variables <- c(scz_pheno_cols, gtex_eqtl_slope_cols, brainvar_eqtl_slope_cols,
                   eqtl_ind_cols, size_cols, loeuf_cols, biotype_cols, wgcna_cols)
ea_variables <- c(ea_pheno_cols, gtex_eqtl_slope_cols, brainvar_eqtl_slope_cols,
                  eqtl_ind_cols, size_cols, loeuf_cols, biotype_cols, wgcna_cols)

# Use the best settings for alpha = 0.01
model_args <- list("nrounds100md1eta02gamma0" = list("nrounds" = 100,
                                                     "max_depth" = 1,
                                                     "eta" = 0.02,
                                                     "gamma" = 0,
                                                     "verbose" = 0,
                                                     "nthread" = 5),
                   "nrounds100md2eta03gamma0" = list("nrounds" = 100,
                                                     "max_depth" = 2,
                                                     "eta" = 0.03,
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
        "data/adapt_results/positional/asd_rsquared25_xgb_cv.rds")

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
        "data/adapt_results/positional/scz_rsquared25_xgb_cv.rds")

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
        "data/adapt_results/positional/ea_rsquared25_xgb_cv.rds")



# Next for r^2 = 0.5 ------------------------------------------------------

model_data <-
  read_csv("data/adapt_model_data/positional/rsquared50_model_data.csv")

# Update the settings
model_args <- list("nrounds100md1eta02gamma0" = list("nrounds" = 100,
                                                     "max_depth" = 1,
                                                     "eta" = 0.02,
                                                     "gamma" = 0,
                                                     "verbose" = 0,
                                                     "nthread" = 5),
                   "nrounds100md2eta02gamma0" = list("nrounds" = 100,
                                                     "max_depth" = 2,
                                                     "eta" = 0.02,
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
        "data/adapt_results/positional/asd_rsquared50_xgb_cv.rds")

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
        "data/adapt_results/positional/scz_rsquared50_xgb_cv.rds")

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
        "data/adapt_results/positional/ea_rsquared50_xgb_cv.rds")

