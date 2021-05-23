# PURPOSE: Compute the Monte Carlo based p-values using the quadratic test statistics

# Positional + eSNPs version

library(tidyverse)
library(data.table)
library(Rcpp)
library(RcppArmadillo)

# ------------------------------------------------------------------------------

# Source the Rcpp functions:
sourceCpp("code/s05_create_gene_locus_level_data/sim_quad_test_stat_fast.cpp")

# Will lazily do this for each threshold

# R-squared .25 ----------------------------------------------------------------

# Load the tidy SNP gene-set data:
tidy_snp_locus_data <-
  read_csv("data/tidy_snp_locus/positional_esnps/agglom_rsquared25_ld_loci.csv")

# Define the flip it function now:
flip_matrix <- function(x) {
  zero2 <- which(x == 0)
  two0 <- which(x == 2)
  x[zero2] <- 2
  x[two0] <- 0
  return(x)
}

# Will loop through each chromosome
all_chr <- unique(tidy_snp_locus_data$chr)


for (chr_i in all_chr) {

  # Get the genotype data for the chromosome
  numeric_snp_ref_genotypes <- as.matrix(
    data.table::fread(
      paste0("data/ld_reference/chr_level/chr",
             chr_i, "_ref_genotypes.csv")))

  chr_tidy_snp_locus_data <- tidy_snp_locus_data %>%
    dplyr::filter(chr == chr_i) %>%
    dplyr::select(-chr)

  # The gene set IDs
  group_ids <- unique(chr_tidy_snp_locus_data$ld_locus_id)

  # Now generate results in parallel for different gene-sets in the CHR
  tidy_locus_summary <-
    mclapply(group_ids, mc.cores = 20,
             function(group_x) {

               # Filter to the SNPs belonging to group X in the dataset:
               group_x_snp_rows <- chr_tidy_snp_locus_data %>%
                 filter(ld_locus_id == group_x)

               # If only 1 SNP then just make a dataset with
               # that SNP's single p-value and test statistics
               # just for reference:
               if (nrow(group_x_snp_rows) == 1) {

                 group_results <- group_x_snp_rows %>%
                   dplyr::select(ld_locus_id,
                                 scz_z_squared, scz_p,
                                 asd_z_squared, asd_p,
                                 ea_z_squared, ea_p) %>%
                   rename(scz_quad_stat = scz_z_squared,
                          scz_quad_pval = scz_p,
                          asd_quad_stat = asd_z_squared,
                          asd_quad_pval = asd_p,
                          ea_quad_stat = ea_z_squared,
                          ea_quad_pval = ea_p) %>%
                   mutate(n_snps = 1)

               } else {

                 # Pull the SNPs:
                 group_x_snps <- pull(group_x_snp_rows,
                                      hg19_id)

                 # Only look at these SNPs' genotype data - with the genotype
                 # columns in the correct order:
                 group_geno_ref_matrix <- numeric_snp_ref_genotypes[, group_x_snps]
                 # Flip the matrix:
                 group_geno_ref_matrix <- flip_matrix(group_geno_ref_matrix)

                 # Compute the group's correlation matrix:
                 group_cor_matrix <- compute_cor_matrix(group_geno_ref_matrix)

                 # Next compute the eigenvalues:
                 group_ev <-  compute_eig_vals(group_cor_matrix)

                 # Get a truncated version above a threshold (just 1e-15 for now)
                 trunc_group_ev <- group_ev[which(group_ev > 1e-15)]

                 group_results <- group_x_snp_rows %>%
                   group_by(ld_locus_id) %>%
                   summarize(scz_quad_stat = sum(scz_z_squared),
                             asd_quad_stat = sum(asd_z_squared),
                             ea_quad_stat = sum(ea_z_squared),
                             n_snps = n()) %>%
                   ungroup()

                 # Generate the test stats as a vector, using 1 million
                 quad_test_stat_sims <- sim_ev_quad_stats(trunc_group_ev,
                                                            1000000)

                 # Now many sims are greater than the observed for each test_stat:
                 n_sims_greater_1_scz <-
                   sum(as.numeric(quad_test_stat_sims > group_results$scz_quad_stat))

                 n_sims_greater_1_asd <-
                   sum(as.numeric(quad_test_stat_sims > group_results$asd_quad_stat))

                 n_sims_greater_1_ea <-
                   sum(as.numeric(quad_test_stat_sims > group_results$ea_quad_stat))


                 # Do it again:
                 quad_test_stat_sims <- sim_ev_quad_stats(trunc_group_ev,
                                                            1000000)
                 n_sims_greater_2_scz <-
                   sum(as.numeric(quad_test_stat_sims > group_results$scz_quad_stat))

                 n_sims_greater_2_asd <-
                   sum(as.numeric(quad_test_stat_sims > group_results$asd_quad_stat))

                 n_sims_greater_2_ea <-
                   sum(as.numeric(quad_test_stat_sims > group_results$ea_quad_stat))

                 # Get the group results with this dataset:
                 group_results <- group_results %>%
                   mutate(scz_quad_pval = (1 + n_sims_greater_1_scz + n_sims_greater_2_scz) / (2000000 + 1),
                          asd_quad_pval = (1 + n_sims_greater_1_asd + n_sims_greater_2_asd) / (2000000 + 1),
                          ea_quad_pval = (1 + n_sims_greater_1_ea + n_sims_greater_2_ea) / (2000000 + 1))  %>%
                   dplyr::select(ld_locus_id,
                                 scz_quad_stat,
                                 scz_quad_pval,
                                 asd_quad_stat,
                                 asd_quad_pval,
                                 ea_quad_stat,
                                 ea_quad_pval,
                                 n_snps)


               }
               return(group_results)

             }) %>%
    bind_rows() %>%
    mutate(chr = chr_i)
  # Save

  write_csv(tidy_locus_summary,
            paste0("data/locus_level_pvals/positional_esnps/rsquared25/chr",
                   chr_i, "_locus_summary.csv"))

  print(paste0("Finished CHR", chr_i, "!"))


}


# R-squared .50 ----------------------------------------------------------------

# Load the tidy SNP gene-set data:
tidy_snp_locus_data <-
  read_csv("data/tidy_snp_locus/positional_esnps/agglom_rsquared50_ld_loci.csv")

# Will loop through each chromosome
all_chr <- unique(tidy_snp_locus_data$chr)


for (chr_i in all_chr) {

  # Get the genotype data for the chromosome
  numeric_snp_ref_genotypes <- as.matrix(
    data.table::fread(
      paste0("data/ld_reference/chr_level/chr",
             chr_i, "_ref_genotypes.csv")))

  chr_tidy_snp_locus_data <- tidy_snp_locus_data %>%
    dplyr::filter(chr == chr_i) %>%
    dplyr::select(-chr)

  # The gene set IDs
  group_ids <- unique(chr_tidy_snp_locus_data$ld_locus_id)

  # Now generate results in parallel for different gene-sets in the CHR
  tidy_locus_summary <-
    mclapply(group_ids, mc.cores = 20,
             function(group_x) {

               # Filter to the SNPs belonging to group X in the dataset:
               group_x_snp_rows <- chr_tidy_snp_locus_data %>%
                 filter(ld_locus_id == group_x)

               # If only 1 SNP then just make a dataset with
               # that SNP's single p-value and test statistics
               # just for reference:
               if (nrow(group_x_snp_rows) == 1) {

                 group_results <- group_x_snp_rows %>%
                   dplyr::select(ld_locus_id,
                                 scz_z_squared, scz_p,
                                 asd_z_squared, asd_p,
                                 ea_z_squared, ea_p) %>%
                   rename(scz_quad_stat = scz_z_squared,
                          scz_quad_pval = scz_p,
                          asd_quad_stat = asd_z_squared,
                          asd_quad_pval = asd_p,
                          ea_quad_stat = ea_z_squared,
                          ea_quad_pval = ea_p) %>%
                   mutate(n_snps = 1)

               } else {

                 # Pull the SNPs:
                 group_x_snps <- pull(group_x_snp_rows,
                                      hg19_id)

                 # Only look at these SNPs' genotype data - with the genotype
                 # columns in the correct order:
                 group_geno_ref_matrix <- numeric_snp_ref_genotypes[, group_x_snps]
                 # Flip the matrix:
                 group_geno_ref_matrix <- flip_matrix(group_geno_ref_matrix)

                 # Compute the group's correlation matrix:
                 group_cor_matrix <- compute_cor_matrix(group_geno_ref_matrix)

                 # Next compute the eigenvalues:
                 group_ev <-  compute_eig_vals(group_cor_matrix)

                 # Get a truncated version above a threshold (just 1e-15 for now)
                 trunc_group_ev <- group_ev[which(group_ev > 1e-15)]

                 group_results <- group_x_snp_rows %>%
                   group_by(ld_locus_id) %>%
                   summarize(scz_quad_stat = sum(scz_z_squared),
                             asd_quad_stat = sum(asd_z_squared),
                             ea_quad_stat = sum(ea_z_squared),
                             n_snps = n()) %>%
                   ungroup()

                 # Generate the test stats as a vector, using 1 million
                 quad_test_stat_sims <- sim_ev_quad_stats(trunc_group_ev,
                                                            1000000)

                 # Now many sims are greater than the observed for each test_stat:
                 n_sims_greater_1_scz <-
                   sum(as.numeric(quad_test_stat_sims > group_results$scz_quad_stat))

                 n_sims_greater_1_asd <-
                   sum(as.numeric(quad_test_stat_sims > group_results$asd_quad_stat))

                 n_sims_greater_1_ea <-
                   sum(as.numeric(quad_test_stat_sims > group_results$ea_quad_stat))


                 # Do it again:
                 quad_test_stat_sims <- sim_ev_quad_stats(trunc_group_ev,
                                                            1000000)
                 n_sims_greater_2_scz <-
                   sum(as.numeric(quad_test_stat_sims > group_results$scz_quad_stat))

                 n_sims_greater_2_asd <-
                   sum(as.numeric(quad_test_stat_sims > group_results$asd_quad_stat))

                 n_sims_greater_2_ea <-
                   sum(as.numeric(quad_test_stat_sims > group_results$ea_quad_stat))

                 # Get the group results with this dataset:
                 group_results <- group_results %>%
                   mutate(scz_quad_pval = (1 + n_sims_greater_1_scz + n_sims_greater_2_scz) / (2000000 + 1),
                          asd_quad_pval = (1 + n_sims_greater_1_asd + n_sims_greater_2_asd) / (2000000 + 1),
                          ea_quad_pval = (1 + n_sims_greater_1_ea + n_sims_greater_2_ea) / (2000000 + 1))  %>%
                   dplyr::select(ld_locus_id,
                                 scz_quad_stat,
                                 scz_quad_pval,
                                 asd_quad_stat,
                                 asd_quad_pval,
                                 ea_quad_stat,
                                 ea_quad_pval,
                                 n_snps)


               }
               return(group_results)

             }) %>%
    bind_rows() %>%
    mutate(chr = chr_i)
  # Save

  write_csv(tidy_locus_summary,
            paste0("data/locus_level_pvals/positional_esnps/rsquared50/chr",
                   chr_i, "_locus_summary.csv"))

  print(paste0("Finished CHR", chr_i, "!"))


}

# R-squared .75 ----------------------------------------------------------------

# Load the tidy SNP gene-set data:
tidy_snp_locus_data <-
  read_csv("data/tidy_snp_locus/positional_esnps/agglom_rsquared75_ld_loci.csv")

# Define the flip it function now:
flip_matrix <- function(x) {
  zero2 <- which(x == 0)
  two0 <- which(x == 2)
  x[zero2] <- 2
  x[two0] <- 0
  return(x)
}

# Will loop through each chromosome
all_chr <- unique(tidy_snp_locus_data$chr)


for (chr_i in all_chr) {

  # Get the genotype data for the chromosome
  numeric_snp_ref_genotypes <- as.matrix(
    data.table::fread(
      paste0("data/ld_reference/asd/chr_level/chr",
             chr_i, "_ref_genotypes.csv")))

  chr_tidy_snp_locus_data <- tidy_snp_locus_data %>%
    dplyr::filter(chr == chr_i) %>%
    dplyr::select(-chr)

  # The gene set IDs
  group_ids <- unique(chr_tidy_snp_locus_data$ld_locus_id)

  # Now generate results in parallel for different gene-sets in the CHR
  tidy_locus_summary <-
    mclapply(group_ids, mc.cores = 20,
             function(group_x) {

               # Filter to the SNPs belonging to group X in the dataset:
               group_x_snp_rows <- chr_tidy_snp_locus_data %>%
                 filter(ld_locus_id == group_x)

               # If only 1 SNP then just make a dataset with
               # that SNP's single p-value and test statistics
               # just for reference:
               if (nrow(group_x_snp_rows) == 1) {

                 group_results <- group_x_snp_rows %>%
                   dplyr::select(ld_locus_id,
                                 scz_z_squared, scz_p,
                                 asd_z_squared, asd_p,
                                 ea_z_squared, ea_p) %>%
                   rename(scz_quad_stat = scz_z_squared,
                          scz_quad_pval = scz_p,
                          asd_quad_stat = asd_z_squared,
                          asd_quad_pval = asd_p,
                          ea_quad_stat = ea_z_squared,
                          ea_quad_pval = ea_p) %>%
                   mutate(n_snps = 1)

               } else {

                 # Pull the SNPs:
                 group_x_snps <- pull(group_x_snp_rows,
                                      hg19_id)

                 # Only look at these SNPs' genotype data - with the genotype
                 # columns in the correct order:
                 group_geno_ref_matrix <- numeric_snp_ref_genotypes[, group_x_snps]
                 # Flip the matrix:
                 group_geno_ref_matrix <- flip_matrix(group_geno_ref_matrix)

                 # Compute the group's correlation matrix:
                 group_cor_matrix <- compute_cor_matrix(group_geno_ref_matrix)

                 # Next compute the eigenvalues:
                 group_ev <-  compute_eig_vals(group_cor_matrix)

                 # Get a truncated version above a threshold (just 1e-15 for now)
                 trunc_group_ev <- group_ev[which(group_ev > 1e-15)]

                 group_results <- group_x_snp_rows %>%
                   group_by(ld_locus_id) %>%
                   summarize(scz_quad_stat = sum(scz_z_squared),
                             asd_quad_stat = sum(asd_z_squared),
                             ea_quad_stat = sum(ea_z_squared),
                             n_snps = n()) %>%
                   ungroup()

                 # Generate the test stats as a vector, using 1 million
                 quad_test_stat_sims <- sim_ev_quad_stats(trunc_group_ev,
                                                            1000000)

                 # Now many sims are greater than the observed for each test_stat:
                 n_sims_greater_1_scz <-
                   sum(as.numeric(quad_test_stat_sims > group_results$scz_quad_stat))

                 n_sims_greater_1_asd <-
                   sum(as.numeric(quad_test_stat_sims > group_results$asd_quad_stat))

                 n_sims_greater_1_ea <-
                   sum(as.numeric(quad_test_stat_sims > group_results$ea_quad_stat))


                 # Do it again:
                 quad_test_stat_sims <- sim_ev_quad_stats(trunc_group_ev,
                                                            1000000)
                 n_sims_greater_2_scz <-
                   sum(as.numeric(quad_test_stat_sims > group_results$scz_quad_stat))

                 n_sims_greater_2_asd <-
                   sum(as.numeric(quad_test_stat_sims > group_results$asd_quad_stat))

                 n_sims_greater_2_ea <-
                   sum(as.numeric(quad_test_stat_sims > group_results$ea_quad_stat))

                 # Get the group results with this dataset:
                 group_results <- group_results %>%
                   mutate(scz_quad_pval = (1 + n_sims_greater_1_scz + n_sims_greater_2_scz) / (2000000 + 1),
                          asd_quad_pval = (1 + n_sims_greater_1_asd + n_sims_greater_2_asd) / (2000000 + 1),
                          ea_quad_pval = (1 + n_sims_greater_1_ea + n_sims_greater_2_ea) / (2000000 + 1))  %>%
                   dplyr::select(ld_locus_id,
                                 scz_quad_stat,
                                 scz_quad_pval,
                                 asd_quad_stat,
                                 asd_quad_pval,
                                 ea_quad_stat,
                                 ea_quad_pval,
                                 n_snps)


               }
               return(group_results)

             }) %>%
    bind_rows() %>%
    mutate(chr = chr_i)
  # Save

  write_csv(tidy_locus_summary,
            paste0("data/locus_level_pvals/positional_esnps/rsquared75/chr",
                   chr_i, "_locus_summary.csv"))

  print(paste0("Finished CHR", chr_i, "!"))


}

