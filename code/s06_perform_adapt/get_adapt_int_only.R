# PURPOSE: Generate the AdaPT intercept-only results

library(tidyverse)
library(adaptMT)

# Helper function to paste for results:
get_model_data_file_path <- function(assign_type, ld_threshold) {
  paste0("data/adapt_model_data/", assign_type,
         "/rsquared", ld_threshold * 100, "_model_data.csv")
}

# Load the AdaPT datasets -------------------------------------------------

# r^2 = 0.25
pos_esnps_25_model_data <-
  read_csv(get_model_data_file_path("positional_esnps", 0.25))
pos_25_model_data <-
  read_csv(get_model_data_file_path("positional", 0.25))

# r^2 = 0.50
pos_esnps_50_model_data <-
  read_csv(get_model_data_file_path("positional_esnps", 0.5))
pos_50_model_data <-
  read_csv(get_model_data_file_path("positional", 0.5))

# r^2 = 0.75
pos_esnps_75_model_data <-
  read_csv(get_model_data_file_path("positional_esnps", 0.75))
pos_75_model_data <-
  read_csv(get_model_data_file_path("positional", 0.75))

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
        "data/adapt_results/positional/asd/rsquared25_int_only.rds")

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
        "data/adapt_results/positional_esnps/asd/rsquared25_int_only.rds")

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
        "data/adapt_results/positional/asd/rsquared50_int_only.rds")

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
        "data/adapt_results/positional_esnps/asd/rsquared50_int_only.rds")

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
        "data/adapt_results/positional/asd/rsquared75_int_only.rds")

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
        "data/adapt_results/positional_esnps/asd/rsquared75_int_only.rds")

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
        "data/adapt_results/positional/scz/rsquared25_int_only.rds")

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
        "data/adapt_results/positional_esnps/scz/rsquared25_int_only.rds")

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
        "data/adapt_results/positional/scz/rsquared50_int_only.rds")

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
        "data/adapt_results/positional_esnps/scz/rsquared50_int_only.rds")

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
        "data/adapt_results/positional/scz/rsquared75_int_only.rds")

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
        "data/adapt_results/positional_esnps/scz/rsquared75_int_only.rds")


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
        "data/adapt_results/positional/ea/rsquared25_int_only.rds")

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
        "data/adapt_results/positional_esnps/ea/rsquared25_int_only.rds")

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
        "data/adapt_results/positional/ea/rsquared50_int_only.rds")

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
        "data/adapt_results/positional_esnps/ea/rsquared50_int_only.rds")

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
        "data/adapt_results/positional/ea/rsquared75_int_only.rds")

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
        "data/adapt_results/positional_esnps/ea/rsquared75_int_only.rds")


