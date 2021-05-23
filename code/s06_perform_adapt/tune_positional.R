# PURPOSE: Initialize fake SCZ tuning results the Positional data

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

# Other phenotypes:
pheno_cols <- c("cap_ea_quad_z", "cap_asd_quad_z")

# Now each of the wgcna_label columns:
wgcna_cols <- stringr::str_subset(colnames(model_data),
                                  "any_gene_")

# Make a vector with all of these variables:
scz_variables <- c(pheno_cols, gtex_eqtl_slope_cols, brainvar_eqtl_slope_cols,
                   eqtl_ind_cols, size_cols, loeuf_cols, biotype_cols, wgcna_cols)


# Create a table of XGBoost settings
model_settings <- expand.grid("nrounds" = seq(100, 450, by = 50),
                              "max_depth" = c(1, 2),
                              "eta" = seq(0.03, 0.06, by = 0.01),
                              "gamma" = c(0)) %>%
  as_tibble()

# Now create the results in parallel for SCZ
scz_results <- mclapply(1:nrow(model_settings), mc.cores = 16,
                        FUN = function(model_i) {

                          # Initialize the model args:
                          model_args <- list("nrounds" = model_settings$nrounds[model_i],
                                             "max_depth" = model_settings$max_depth[model_i],
                                             "eta" = model_settings$eta[model_i],
                                             "gamma" = model_settings$gamma[model_i],
                                             "verbose" = 0,
                                             "nthread" = 3)

                          adapt_results <- adapt_xgboost(as.matrix(
                            model_data[,scz_variables]),
                            model_data$fake_scz_quad_pval,
                            verbose = list(print = FALSE,
                                           fit = FALSE,
                                           ms = FALSE),
                            piargs = model_args,
                            muargs = model_args,
                            s0 = rep(0.05, nrow(model_data)))

                          # Now return a dataframe of results over
                          # target alpha values:
                          target_alphas <- c(.01, .05, .1, .15, .2)

                          # Generate a data frame of results for the vector of alpha values:
                          do.call(rbind,
                                  lapply(target_alphas,
                                         function(alpha) {
                                           # Access the discoveries for alpha:
                                           adapt_disc <- which(adapt_results$qvals <= alpha)

                                           # Return the fdp and power:
                                           return(data.frame("phenotype" = "fake_scz",
                                                             "ld_loci_type" = "pos_rsquared25",
                                                             "method" = "adapt_xgboost",
                                                             "variables" = "all",
                                                             "alpha" = alpha,
                                                             "n_disc" = length(adapt_disc),
                                                             "nrounds" = model_settings$nrounds[model_i],
                                                             "max_depth" = model_settings$max_depth[model_i],
                                                             "eta" = model_settings$eta[model_i],
                                                             "gamma" = model_settings$gamma[model_i]))
                                         }))

                        }) %>%
  bind_rows()

# Save the dataframe summarizing results:
readr::write_csv(scz_results,
                 "data/adapt_results/positional/fake_scz_tune/rsquared25_tune_results.csv")


# rsquared .50 ------------------------------------------------------------

# Repeat the same steps above for 0.50 threshold dataset

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

# Other phenotypes:
pheno_cols <- c("cap_ea_quad_z", "cap_asd_quad_z")

# Now each of the wgcna_label columns:
wgcna_cols <- stringr::str_subset(colnames(model_data),
                                  "any_gene_")

# Make a vector with all of these variables:
scz_variables <- c(pheno_cols, gtex_eqtl_slope_cols, brainvar_eqtl_slope_cols,
                   eqtl_ind_cols, size_cols, loeuf_cols, biotype_cols, wgcna_cols)


model_settings <- expand.grid("nrounds" = seq(100, 450, by = 50),
                              "max_depth" = c(1, 2),
                              "eta" = seq(0.03, 0.06, by = 0.01),
                              "gamma" = c(0)) %>%
  as_tibble()

# Now create the results in parallel for SCZ
scz_results <- mclapply(1:nrow(model_settings), mc.cores = 16,
                        FUN = function(model_i) {

                          # Initialize the model args:
                          model_args <- list("nrounds" = model_settings$nrounds[model_i],
                                             "max_depth" = model_settings$max_depth[model_i],
                                             "eta" = model_settings$eta[model_i],
                                             "gamma" = model_settings$gamma[model_i],
                                             "verbose" = 0,
                                             "nthread" = 3)

                          adapt_results <- adapt_xgboost(as.matrix(
                            model_data[,scz_variables]),
                            model_data$fake_scz_quad_pval,
                            verbose = list(print = FALSE,
                                           fit = FALSE,
                                           ms = FALSE),
                            piargs = model_args,
                            muargs = model_args,
                            s0 = rep(0.05, nrow(model_data)))

                          # Now return a dataframe of results over
                          # target alpha values:
                          target_alphas <- c(.01, .05, .1, .15, .2)

                          # Generate a data frame of results for the vector of alpha values:
                          do.call(rbind,
                                  lapply(target_alphas,
                                         function(alpha) {
                                           # Access the discoveries for alpha:
                                           adapt_disc <- which(adapt_results$qvals <= alpha)

                                           # Return the fdp and power:
                                           return(data.frame("phenotype" = "fake_scz",
                                                             "ld_loci_type" = "pos_rsquared50",
                                                             "method" = "adapt_xgboost",
                                                             "variables" = "all",
                                                             "alpha" = alpha,
                                                             "n_disc" = length(adapt_disc),
                                                             "nrounds" = model_settings$nrounds[model_i],
                                                             "max_depth" = model_settings$max_depth[model_i],
                                                             "eta" = model_settings$eta[model_i],
                                                             "gamma" = model_settings$gamma[model_i]))
                                         }))

                        }) %>%
  bind_rows()

# Save the dataframe summarizing results:
readr::write_csv(scz_results,
                 "data/adapt_results/positional/fake_scz_tune/rsquared50_tune_results.csv")


# rsquared .75 ------------------------------------------------------------

# Repeat the same steps above for 0.75 threshold dataset

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

# Other phenotypes:
pheno_cols <- c("cap_ea_quad_z", "cap_asd_quad_z")

# Now each of the wgcna_label columns:
wgcna_cols <- stringr::str_subset(colnames(model_data),
                                  "any_gene_")

# Make a vector with all of these variables:
scz_variables <- c(pheno_cols, gtex_eqtl_slope_cols, brainvar_eqtl_slope_cols,
                   eqtl_ind_cols, size_cols, loeuf_cols, biotype_cols, wgcna_cols)


model_settings <- expand.grid("nrounds" = seq(100, 450, by = 50),
                              "max_depth" = c(1, 2),
                              "eta" = seq(0.03, 0.06, by = 0.01),
                              "gamma" = c(0)) %>%
  as_tibble()

# Now create the results in parallel for SCZ
scz_results <- mclapply(1:nrow(model_settings), mc.cores = 16,
                        FUN = function(model_i) {

                          # Initialize the model args:
                          model_args <- list("nrounds" = model_settings$nrounds[model_i],
                                             "max_depth" = model_settings$max_depth[model_i],
                                             "eta" = model_settings$eta[model_i],
                                             "gamma" = model_settings$gamma[model_i],
                                             "verbose" = 0,
                                             "nthread" = 3)

                          adapt_results <- adapt_xgboost(as.matrix(
                            model_data[,scz_variables]),
                            model_data$fake_scz_quad_pval,
                            verbose = list(print = FALSE,
                                           fit = FALSE,
                                           ms = FALSE),
                            piargs = model_args,
                            muargs = model_args,
                            s0 = rep(0.05, nrow(model_data)))

                          # Now return a dataframe of results over
                          # target alpha values:
                          target_alphas <- c(.01, .05, .1, .15, .2)

                          # Generate a data frame of results for the vector of alpha values:
                          do.call(rbind,
                                  lapply(target_alphas,
                                         function(alpha) {
                                           # Access the discoveries for alpha:
                                           adapt_disc <- which(adapt_results$qvals <= alpha)

                                           # Return the fdp and power:
                                           return(data.frame("phenotype" = "fake_scz",
                                                             "ld_loci_type" = "pos_rsquared75",
                                                             "method" = "adapt_xgboost",
                                                             "variables" = "all",
                                                             "alpha" = alpha,
                                                             "n_disc" = length(adapt_disc),
                                                             "nrounds" = model_settings$nrounds[model_i],
                                                             "max_depth" = model_settings$max_depth[model_i],
                                                             "eta" = model_settings$eta[model_i],
                                                             "gamma" = model_settings$gamma[model_i]))
                                         }))

                        }) %>%
  bind_rows()

# Save the dataframe summarizing results:
readr::write_csv(scz_results,
                 "data/adapt_results/positional/fake_scz_tune/rsquared75_tune_results.csv")

