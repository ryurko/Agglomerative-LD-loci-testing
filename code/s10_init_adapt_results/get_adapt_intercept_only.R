# PURPOSE: Initialize the intercept-only results for each dataset and phenotype

library(tidyverse)
library(adaptMT)

# Load the AdaPT datasets -------------------------------------------------

pos_esnps_25_model_data <-
  read_csv("data/adapt_model_data/positional_esnps/rsquared25_model_data.csv")

pos_25_model_data <-
  read_csv("data/adapt_model_data/positional/rsquared25_model_data.csv")

pos_esnps_50_model_data <-
  read_csv("data/adapt_model_data/positional_esnps/rsquared50_model_data.csv")

pos_50_model_data <-
  read_csv("data/adapt_model_data/positional/rsquared50_model_data.csv")

pos_esnps_75_model_data <-
  read_csv("data/adapt_model_data/positional_esnps/rsquared75_model_data.csv")

pos_75_model_data <-
  read_csv("data/adapt_model_data/positional/rsquared75_model_data.csv")

# Initialize ASD results --------------------------------------------------

# r^2 = .25

# Positional
asd_pos_25_int_only <-
  adapt_glm(x = pos_25_model_data,
            pvals = pos_25_model_data$asd_vegas_pval,
            pi_formulas = "1",
            mu_formulas = "1",
            verbose = list(print = FALSE,
                           fit = FALSE,
                           ms = FALSE),
            s0 = rep(0.05, nrow(pos_25_model_data)))
length(which(asd_pos_25_int_only$qvals <= 0.05))
saveRDS(asd_pos_25_int_only,
        "data/adapt_results/positional/asd_rsquared25_int_only.rds")

# Positional + eSNPs
asd_pos_esnps_25_int_only <-
  adapt_glm(x = pos_esnps_25_model_data,
            pvals = pos_esnps_25_model_data$asd_vegas_pval,
            pi_formulas = "1",
            mu_formulas = "1",
            verbose = list(print = FALSE,
                           fit = FALSE,
                           ms = FALSE),
            s0 = rep(0.05, nrow(pos_esnps_25_model_data)))
length(which(asd_pos_esnps_25_int_only$qvals <= 0.05))
saveRDS(asd_pos_esnps_25_int_only,
        "data/adapt_results/positional_esnps/asd_rsquared25_int_only.rds")

# r^2 = .50

# Positional
asd_pos_50_int_only <-
  adapt_glm(x = pos_50_model_data,
            pvals = pos_50_model_data$asd_vegas_pval,
            pi_formulas = "1",
            mu_formulas = "1",
            verbose = list(print = FALSE,
                           fit = FALSE,
                           ms = FALSE),
            s0 = rep(0.05, nrow(pos_50_model_data)))
length(which(asd_pos_50_int_only$qvals <= 0.05))
saveRDS(asd_pos_50_int_only,
        "data/adapt_results/positional/asd_rsquared50_int_only.rds")

# Positional + eSNPs
asd_pos_esnps_50_int_only <-
  adapt_glm(x = pos_esnps_50_model_data,
            pvals = pos_esnps_50_model_data$asd_vegas_pval,
            pi_formulas = "1",
            mu_formulas = "1",
            verbose = list(print = FALSE,
                           fit = FALSE,
                           ms = FALSE),
            s0 = rep(0.05, nrow(pos_esnps_50_model_data)))
length(which(asd_pos_esnps_50_int_only$qvals <= 0.05))
saveRDS(asd_pos_esnps_50_int_only,
        "data/adapt_results/positional_esnps/asd_rsquared50_int_only.rds")

# r^2 = .75

# Positional
asd_pos_75_int_only <-
  adapt_glm(x = pos_75_model_data,
            pvals = pos_75_model_data$asd_vegas_pval,
            pi_formulas = "1",
            mu_formulas = "1",
            verbose = list(print = FALSE,
                           fit = FALSE,
                           ms = FALSE),
            s0 = rep(0.05, nrow(pos_75_model_data)))
length(which(asd_pos_75_int_only$qvals <= 0.05))
saveRDS(asd_pos_75_int_only,
        "data/adapt_results/positional/asd_rsquared75_int_only.rds")

# Positional + eSNPs
asd_pos_esnps_75_int_only <-
  adapt_glm(x = pos_esnps_75_model_data,
            pvals = pos_esnps_75_model_data$asd_vegas_pval,
            pi_formulas = "1",
            mu_formulas = "1",
            verbose = list(print = FALSE,
                           fit = FALSE,
                           ms = FALSE),
            s0 = rep(0.05, nrow(pos_esnps_75_model_data)))
length(which(asd_pos_esnps_75_int_only$qvals <= 0.05))
saveRDS(asd_pos_esnps_75_int_only,
        "data/adapt_results/positional_esnps/asd_rsquared75_int_only.rds")

# Initialize SCZ results --------------------------------------------------

# r^2 = .25

# Positional
scz_pos_25_int_only <-
  adapt_glm(x = pos_25_model_data,
            pvals = pos_25_model_data$scz_vegas_pval,
            pi_formulas = "1",
            mu_formulas = "1",
            verbose = list(print = FALSE,
                           fit = FALSE,
                           ms = FALSE),
            s0 = rep(0.05, nrow(pos_25_model_data)))
length(which(scz_pos_25_int_only$qvals <= 0.05))
saveRDS(scz_pos_25_int_only,
        "data/adapt_results/positional/scz_rsquared25_int_only.rds")

# Positional + eSNPs
scz_pos_esnps_25_int_only <-
  adapt_glm(x = pos_esnps_25_model_data,
            pvals = pos_esnps_25_model_data$scz_vegas_pval,
            pi_formulas = "1",
            mu_formulas = "1",
            verbose = list(print = FALSE,
                           fit = FALSE,
                           ms = FALSE),
            s0 = rep(0.05, nrow(pos_esnps_25_model_data)))
length(which(scz_pos_esnps_25_int_only$qvals <= 0.05))
saveRDS(scz_pos_esnps_25_int_only,
        "data/adapt_results/positional_esnps/scz_rsquared25_int_only.rds")

# r^2 = .50

# Positional
scz_pos_50_int_only <-
  adapt_glm(x = pos_50_model_data,
            pvals = pos_50_model_data$scz_vegas_pval,
            pi_formulas = "1",
            mu_formulas = "1",
            verbose = list(print = FALSE,
                           fit = FALSE,
                           ms = FALSE),
            s0 = rep(0.05, nrow(pos_50_model_data)))
length(which(scz_pos_50_int_only$qvals <= 0.05))
saveRDS(scz_pos_50_int_only,
        "data/adapt_results/positional/scz_rsquared50_int_only.rds")

# Positional + eSNPs
scz_pos_esnps_50_int_only <-
  adapt_glm(x = pos_esnps_50_model_data,
            pvals = pos_esnps_50_model_data$scz_vegas_pval,
            pi_formulas = "1",
            mu_formulas = "1",
            verbose = list(print = FALSE,
                           fit = FALSE,
                           ms = FALSE),
            s0 = rep(0.05, nrow(pos_esnps_50_model_data)))
length(which(scz_pos_esnps_50_int_only$qvals <= 0.05))
saveRDS(scz_pos_esnps_50_int_only,
        "data/adapt_results/positional_esnps/scz_rsquared50_int_only.rds")

# r^2 = .75

# Positional
scz_pos_75_int_only <-
  adapt_glm(x = pos_75_model_data,
            pvals = pos_75_model_data$scz_vegas_pval,
            pi_formulas = "1",
            mu_formulas = "1",
            verbose = list(print = FALSE,
                           fit = FALSE,
                           ms = FALSE),
            s0 = rep(0.05, nrow(pos_75_model_data)))
length(which(scz_pos_75_int_only$qvals <= 0.05))
saveRDS(scz_pos_75_int_only,
        "data/adapt_results/positional/scz_rsquared75_int_only.rds")

# Positional + eSNPs
scz_pos_esnps_75_int_only <-
  adapt_glm(x = pos_esnps_75_model_data,
            pvals = pos_esnps_75_model_data$scz_vegas_pval,
            pi_formulas = "1",
            mu_formulas = "1",
            verbose = list(print = FALSE,
                           fit = FALSE,
                           ms = FALSE),
            s0 = rep(0.05, nrow(pos_esnps_75_model_data)))
length(which(scz_pos_esnps_75_int_only$qvals <= 0.05))
saveRDS(scz_pos_esnps_75_int_only,
        "data/adapt_results/positional_esnps/scz_rsquared75_int_only.rds")

# Initialize EA results ---------------------------------------------------

# r^2 = .25

# Positional
ea_pos_25_int_only <-
  adapt_glm(x = pos_25_model_data,
            pvals = pos_25_model_data$ea_vegas_pval,
            pi_formulas = "1",
            mu_formulas = "1",
            verbose = list(print = FALSE,
                           fit = FALSE,
                           ms = FALSE),
            s0 = rep(0.05, nrow(pos_25_model_data)))
length(which(ea_pos_25_int_only$qvals <= 0.05))
saveRDS(ea_pos_25_int_only,
        "data/adapt_results/positional/ea_rsquared25_int_only.rds")

# Positional + eSNPs
ea_pos_esnps_25_int_only <-
  adapt_glm(x = pos_esnps_25_model_data,
            pvals = pos_esnps_25_model_data$ea_vegas_pval,
            pi_formulas = "1",
            mu_formulas = "1",
            verbose = list(print = FALSE,
                           fit = FALSE,
                           ms = FALSE),
            s0 = rep(0.05, nrow(pos_esnps_25_model_data)))
length(which(ea_pos_esnps_25_int_only$qvals <= 0.05))
saveRDS(ea_pos_esnps_25_int_only,
        "data/adapt_results/positional_esnps/ea_rsquared25_int_only.rds")

# r^2 = .50

# Positional
ea_pos_50_int_only <-
  adapt_glm(x = pos_50_model_data,
            pvals = pos_50_model_data$ea_vegas_pval,
            pi_formulas = "1",
            mu_formulas = "1",
            verbose = list(print = FALSE,
                           fit = FALSE,
                           ms = FALSE),
            s0 = rep(0.05, nrow(pos_50_model_data)))
length(which(ea_pos_50_int_only$qvals <= 0.05))
saveRDS(ea_pos_50_int_only,
        "data/adapt_results/positional/ea_rsquared50_int_only.rds")

# Positional + eSNPs
ea_pos_esnps_50_int_only <-
  adapt_glm(x = pos_esnps_50_model_data,
            pvals = pos_esnps_50_model_data$ea_vegas_pval,
            pi_formulas = "1",
            mu_formulas = "1",
            verbose = list(print = FALSE,
                           fit = FALSE,
                           ms = FALSE),
            s0 = rep(0.05, nrow(pos_esnps_50_model_data)))
length(which(ea_pos_esnps_50_int_only$qvals <= 0.05))
saveRDS(ea_pos_esnps_50_int_only,
        "data/adapt_results/positional_esnps/ea_rsquared50_int_only.rds")

# r^2 = .75

# Positional
ea_pos_75_int_only <-
  adapt_glm(x = pos_75_model_data,
            pvals = pos_75_model_data$ea_vegas_pval,
            pi_formulas = "1",
            mu_formulas = "1",
            verbose = list(print = FALSE,
                           fit = FALSE,
                           ms = FALSE),
            s0 = rep(0.05, nrow(pos_75_model_data)))
length(which(ea_pos_75_int_only$qvals <= 0.05))
saveRDS(ea_pos_75_int_only,
        "data/adapt_results/positional/ea_rsquared75_int_only.rds")

# Positional + eSNPs
ea_pos_esnps_75_int_only <-
  adapt_glm(x = pos_esnps_75_model_data,
            pvals = pos_esnps_75_model_data$ea_vegas_pval,
            pi_formulas = "1",
            mu_formulas = "1",
            verbose = list(print = FALSE,
                           fit = FALSE,
                           ms = FALSE),
            s0 = rep(0.05, nrow(pos_esnps_75_model_data)))
length(which(ea_pos_esnps_75_int_only$qvals <= 0.05))
saveRDS(ea_pos_esnps_75_int_only,
        "data/adapt_results/positional_esnps/ea_rsquared75_int_only.rds")

