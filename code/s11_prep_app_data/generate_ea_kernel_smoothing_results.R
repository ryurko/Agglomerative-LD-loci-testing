# PURPOSE: Generate the kernel smoothing results for the positional SNPs at
#          both the gene and LD loci level, with null results generated at
#          the LD loci level. Results are generated across all types of
#          SNP-to-gene assignments

# EA VERSION

library(tidyverse)
library(np)

# Write function for kernel zooming interpolation -------------------------

# Define a function that takes in a dataset of SNP results to generate kernel
# smoothing for, and returns a dataset with the results interpolated along
# a grid of evenly spaced points

#' Perform kernel smoothing interpolation
#'
#' @param bp Vector of base pair positions
#' @param obs_zsquared Observed z-squared stats to smooth over
#' @param h Kernel bandwidth to use in kernel regression
#' @param n_points Number of points to generate for interpolation

get_smooth_interpolation <- function(bp, obs_zsquared, h, n_points) {
  
  # Create grid points to evaluate at
  bp_grid <- seq(min(bp), max(bp), length.out = n_points)
  
  # Round to nearest digit and take unique:
  bp_grid <- unique(round(bp_grid))
  
  # Combine this with the actual data points:
  bp_grid <- sort(unique(c(bp_grid, bp)))
  
  # Return a tibble with the smooth values along this grid using the
  # actual values to train it:
  tibble(smooth_zsquared = npreg(txdat = bp,
                                 tydat = obs_zsquared,
                                 exdat = bp_grid,
                                 bws = h)$mean) %>%
    mutate(bp = bp_grid)
  
}


# Write a function for LD loci level smoothing results -------------------


get_ld_loci_level_smooth_results <- function(snp_ld_loci_data, bw_data, stat_column,
                                             n_ld_loci_points = 1000) {
  
  # Assume the columns are there...
  
  map_dfr(unique(snp_ld_loci_data$ld_loci_id),
          function(ld_loci_i) {
            
            # Get the LD loci data and bandwidth:
            ld_loci_data <- snp_ld_loci_data %>%
              filter(ld_loci_id == ld_loci_i)
            ld_loci_bw <- bw_data %>%
              filter(ld_loci_id == ld_loci_i) %>%
              pull(gcv_bw)
            
            # If there is a bandwidth value then perform the interpolation
            # otherwise just return the observed values (this is so that
            # none of the LD locis are removed from the data for display
            # in the app)
            if (is.na(ld_loci_bw)) {
              smooth_results <- tibble(bp = ld_loci_data$bp,
                                       smooth_zsquared = ld_loci_data[[stat_column]])
            } else {
              smooth_results <- get_smooth_interpolation(ld_loci_data$bp,
                                                         ld_loci_data[[stat_column]],
                                                         ld_loci_bw, n_ld_loci_points)
            }
            
            smooth_results %>%
              mutate(ld_loci_id = ld_loci_i) %>%
              return()
            
          })
  
}



# Write a similar function for gene-level ---------------------------------


get_gene_level_smooth_results <- function(snp_gene_ld_loci_data, bw_data,
                                          stat_column, n_gene_points = 100) {
  
  # Assume the columns are there...
  
  map_dfr(unique(snp_gene_ld_loci_data$ld_loci_id),
          function(ld_loci_i) {
            
            # Get the LD loci data
            ld_loci_data <- snp_gene_ld_loci_data %>%
              filter(ld_loci_id == ld_loci_i)
            
            # Get the LD loci's bandwidth for all:
            ld_loci_bw <- bw_data %>%
              filter(ld_loci_id == ld_loci_i) %>%
              pull(gcv_bw)
            
            # What are the unique genes:
            ld_loci_genes <- ld_loci_data %>%
              pull(ensembl_id) %>%
              unique()
            
            if (is.na(ld_loci_bw)) {
              
              smooth_ld_loci_gene_results <-
                tibble(bp = ld_loci_data$bp, smooth_zsquared = ld_loci_data[[stat_column]],
                       ensembl_id = ld_loci_data$ensembl_id)
              
              
            } else {
              
              # Proceed through each gene:
              smooth_ld_loci_gene_results <-
                map_dfr(ld_loci_genes,
                        function(gene_i) {
                          
                          gene_data <- ld_loci_data %>%
                            filter(ensembl_id == gene_i)
                          
                          # Only generate the results if more than 2 SNPs:
                          if (nrow(gene_data) > 2) {
                            smooth_gene_results <-
                              get_smooth_interpolation(gene_data$bp,
                                                       gene_data[[stat_column]],
                                                       ld_loci_bw, n_gene_points)
                          } else {
                            smooth_gene_results <-
                              tibble(bp = gene_data$bp,
                                     smooth_zsquared = gene_data[[stat_column]])
                          }
                          
                          smooth_gene_results %>%
                            mutate(ensembl_id = gene_i) %>%
                            return()
                          
                        })
              
            }
            
            smooth_ld_loci_gene_results %>%
              mutate(ld_loci_id = ld_loci_i) %>%
              return()
            
          })
}


# Load the input positional SNP to LD loci and bandwidth datasets --------

# Positional + eSNPS
pos_esnps_loc_snp_ld_loci_results <-
  read_csv("data/kernel_smoothing/ea/input/positional_esnps/positional_snp_ld_loci_data.csv")
pos_esnps_loc_snp_gene_results <-
  read_csv("data/kernel_smoothing/ea/input/positional_esnps/positional_snp_gene_ld_loci_data.csv")
pos_esnps_bw_data <-
  read_csv("data/kernel_smoothing/ea/input/positional_esnps/ld_loci_gcv_bw.csv")

# Positional
pos_loc_snp_ld_loci_results <-
  read_csv("data/kernel_smoothing/ea/input/positional/positional_snp_ld_loci_data.csv")
pos_loc_snp_gene_results <-
  read_csv("data/kernel_smoothing/ea/input/positional/positional_snp_gene_ld_loci_data.csv")
pos_bw_data <-
  read_csv("data/kernel_smoothing/ea/input/positional/ld_loci_gcv_bw.csv")


# Generate the LD loci and gene-level results ----------------------------

# Positional + eSNPs
pos_esnps_ld_loci_level_smoothing <-
  get_ld_loci_level_smooth_results(pos_esnps_loc_snp_ld_loci_results, pos_esnps_bw_data,
                                   "ea_z_squared", n_ld_loci_points = 1000)

pos_esnps_gene_level_smoothing <-
  get_gene_level_smooth_results(pos_esnps_loc_snp_gene_results, pos_esnps_bw_data,
                                "ea_z_squared", n_gene_points = 100)
# Join info to gene-level and save:
pos_esnps_gene_level_smoothing <- pos_esnps_gene_level_smoothing %>%
  dplyr::rename(smooth_ea_zsquared = smooth_zsquared) %>%
  dplyr::left_join(pos_esnps_loc_snp_gene_results,
                   by = c("bp", "ld_loci_id", "ensembl_id")) %>%
  mutate(is_fake_bp = as.numeric(is.na(hg38_id))) %>%
  # For ease later on in app, save the snp id with the version's hg coordinates
  dplyr::rename(snp_id = hg38_id)
# And save
write_csv(pos_esnps_gene_level_smoothing,
          "data/kernel_smoothing/ea/output/positional_esnps/positional_gene_level_smoothing.csv")

# Repeat these steps for the remaining assignment types

# Positional
pos_ld_loci_level_smoothing <-
  get_ld_loci_level_smooth_results(pos_loc_snp_ld_loci_results, pos_bw_data,
                                   "ea_z_squared", n_ld_loci_points = 1000)

pos_gene_level_smoothing <-
  get_gene_level_smooth_results(pos_loc_snp_gene_results, pos_bw_data,
                                "ea_z_squared", n_gene_points = 100)
pos_gene_level_smoothing <- pos_gene_level_smoothing %>%
  dplyr::rename(smooth_ea_zsquared = smooth_zsquared) %>%
  dplyr::left_join(pos_loc_snp_gene_results,
                   by = c("bp", "ld_loci_id", "ensembl_id")) %>%
  mutate(is_fake_bp = as.numeric(is.na(hg38_id))) %>%
  dplyr::rename(snp_id = hg38_id)
write_csv(pos_gene_level_smoothing,
          "data/kernel_smoothing/ea/output/positional/positional_gene_level_smoothing.csv")


# Write a function to generate null LD loci level results ----------------

library(parallel)
library(snpcombineR)

# Define the flip it function now:
flip_matrix <- function(x) {
  zero2 <- which(x == 0)
  two0 <- which(x == 2)
  x[zero2] <- 2
  x[two0] <- 0
  return(x)
}

get_ld_loci_level_null_results <- function(snp_ld_loci_data, bw_data,
                                           n_ld_loci_points = 1000, n_null_sims = 1000,
                                           sim_percentiles = c(0.50, 0.75, 0.95),
                                           n_cores = 5) {
  
  snp_ld_loci_data <- snp_ld_loci_data %>%
    separate(ld_loci_id, c("chr", "set_id"), sep = "_", remove = FALSE) %>%
    mutate(chr = as.numeric(str_remove(chr, "chr"))) %>%
    dplyr::select(-set_id)
  
  # Generate results over the unique CHR:
  lapply(unique(snp_ld_loci_data$chr),
         function(chr_i) {
           chr_ld_loci_data <- snp_ld_loci_data %>%
             filter(chr == chr_i)
           
           # Load the chromosome genotype matrix with only these SNPs:
           chr_ref_genotypes <- as.matrix(
             data.table::fread(
               paste0("data/ld_reference/chr_level/chr",
                      chr_i, "_ref_genotypes.csv")))[, unique(chr_ld_loci_data$snp_id)]
           
           # Now proceed to generate the results in parallel for each LD loci:
           ld_loci_ids <- unique(chr_ld_loci_data$ld_loci_id)
           
           chr_ld_loci_summary <-
             mclapply(ld_loci_ids, mc.cores = n_cores,
                      function(ld_loci_i) {
                        
                        # Get the LD loci data:
                        ld_loci_data <- chr_ld_loci_data %>%
                          filter(ld_loci_id == ld_loci_i)
                        
                        # Get the LD loci's bandwidth:
                        ld_loci_bw <- bw_data %>%
                          filter(ld_loci_id == ld_loci_i) %>%
                          pull(gcv_bw)
                        
                        # If there is a bandwidth, proceed to generate its smooth results:
                        if (!is.na(ld_loci_bw)) {
                          
                          # Next generate null data
                          ld_loci_null_sims <- t(
                            mvtnorm::rmvnorm(n_null_sims,
                                             sigma =
                                               corpcor::make.positive.definite(
                                                 compute_cor_matrix(
                                                   flip_matrix(
                                                     chr_ref_genotypes[, ld_loci_data$snp_id]
                                                   ))), method = "chol"))
                          #colnames(ld_loci_null_sims) <- paste0("sim", c(1:ncol(ld_loci_null_sims)))
                          
                          # Generate the smooth results for each simulation
                          ld_loci_null_results <-
                            map_dfr(1:ncol(ld_loci_null_sims),
                                    function(sim_i) {
                                      get_smooth_interpolation(
                                        ld_loci_data$bp,
                                        # Make sure to square the z-stats!
                                        (ld_loci_null_sims[, sim_i])^2,
                                        ld_loci_bw, n_ld_loci_points) %>%
                                        mutate(sim_index = sim_i)
                                    }) %>%
                            group_by(bp) %>%
                            summarise(x = list(enframe(quantile(smooth_zsquared,
                                                                probs = sim_percentiles),
                                                       "percentile", "null_zsquared"))) %>%
                            unnest(x) %>%
                            pivot_wider(names_from = percentile,
                                        values_from = null_zsquared,
                                        names_prefix = "null_") %>%
                            janitor::clean_names() %>%
                            mutate(ld_loci_id = ld_loci_i)
                          
                        } else {
                          # Otherwise do not set up positional null results
                          ld_loci_null_results <- NA
                        }
                        
                        return(list("null_results" = ld_loci_null_results))
                        
                      })
           return(chr_ld_loci_summary)
           
           
         }) %>%
    # Then finally return the non-missinld_loci
    map_dfr(function(chr_null_list) {
      result_list <-
        map(1:length(chr_null_list),
            function(ld_loci_i) {
              chr_null_list[[ld_loci_i]]$null_results
            })
      bind_rows(result_list[!is.na(result_list)])
    })
}


# Generate the null results for each  -------------------------------------

pos_esnps_ld_loci_null_results <-
  get_ld_loci_level_null_results(pos_esnps_loc_snp_ld_loci_results, pos_esnps_bw_data,
                                 n_ld_loci_points = 1000, n_null_sims = 1000,
                                 sim_percentiles = c(0.50, 0.75, 0.95),
                                 n_cores = 5)
pos_ld_loci_null_results <-
  get_ld_loci_level_null_results(pos_loc_snp_ld_loci_results, pos_bw_data,
                                 n_ld_loci_points = 1000, n_null_sims = 1000,
                                 sim_percentiles = c(0.50, 0.75, 0.95),
                                 n_cores = 5)

# Join to the LD loci results and save -----------------------------------

# Join the LD loci level null results and LD loci level info 

# Start with Positional + eSNPs
pos_esnps_ld_loci_level_smoothing <- pos_esnps_ld_loci_level_smoothing %>%
  dplyr::rename(smooth_ea_zsquared = smooth_zsquared) %>%
  # Join the null results
  dplyr::left_join(pos_esnps_ld_loci_null_results,
                   by = c("bp", "ld_loci_id")) %>%
  # Join the remaining LD loci info
  dplyr::left_join(pos_esnps_loc_snp_ld_loci_results,
                   by = c("bp", "ld_loci_id")) %>%
  dplyr::rename(snp_id = hg38_id) %>%
  mutate(is_fake_bp = as.numeric(is.na(snp_id)))
write_csv(pos_esnps_ld_loci_level_smoothing,
          "data/kernel_smoothing/ea/output/positional_esnps/positional_ld_loci_level_smoothing.csv")


# Repeat for Positional
pos_ld_loci_level_smoothing <- pos_ld_loci_level_smoothing %>%
  dplyr::rename(smooth_ea_zsquared = smooth_zsquared) %>%
  dplyr::left_join(pos_ld_loci_null_results,
                   by = c("bp", "ld_loci_id")) %>%
  dplyr::left_join(pos_loc_snp_ld_loci_results,
                   by = c("bp", "ld_loci_id")) %>%
  dplyr::rename(snp_id = hg38_id) %>%
  mutate(is_fake_bp = as.numeric(is.na(snp_id)))
write_csv(pos_ld_loci_level_smoothing,
          "data/kernel_smoothing/ea/output/positional/positional_ld_loci_level_smoothing.csv")
