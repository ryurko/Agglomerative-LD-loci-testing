# PURPOSE: Initialize scripts for generating kernel smoothing results

library(np)
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

# Now make a new helper to smooth within the locus clusters
get_locus_cluster_level_smooth_results <-
  function(snp_ld_loci_data, bw_data, stat_column, n_cluster_points = 1000) {

    # Assume the columns are there...

    map_dfr(unique(snp_ld_loci_data$ld_locus_id),
            function(ld_loci_i) {

              # Get the LD loci data:
              ld_loci_data <- snp_ld_loci_data %>%
                filter(ld_locus_id == ld_loci_i)
              ld_loci_bw <- bw_data %>%
                filter(ld_locus_id == ld_loci_i) %>%
                pull(gcv_bw)

              # Now iterate through each cluster:
              map_dfr(unique(ld_loci_data$intra_locus_cluster),
                      function(cluster_i) {

                        # Get the cluster data and bandwidth
                        cluster_data <- ld_loci_data %>%
                          filter(intra_locus_cluster == cluster_i)

                        # If there is a bandwidth value then perform the interpolation
                        # otherwise just return the observed values (this is so that
                        # none of the LD locis are removed from the data for display
                        # in the app)
                        if (is.na(ld_loci_bw) | nrow(cluster_data) < 3) {
                          # Will make these appear as peaks - by adding fake SNPs
                          # to both ends with values of 0
                          smooth_results <- tibble(bp = cluster_data$bp,
                                                   smooth_zsquared = cluster_data[[stat_column]]) %>%
                            bind_rows(tibble(
                              bp = c(min(cluster_data$bp) - 1, max(cluster_data$bp) + 1),
                              smooth_zsquared = rep(0, 2)
                            ))
                        } else {
                          smooth_results <-
                            get_smooth_interpolation(cluster_data$bp,
                                                     cluster_data[[stat_column]],
                                                     ld_loci_bw, n_cluster_points)
                        }
                        smooth_results  %>%
                          mutate(intra_locus_cluster = cluster_i) %>%
                          return()
                      }) %>%
                mutate(ld_locus_id = ld_loci_i)
            })
  }


# Write a similar function for gene-level ---------------------------------


get_gene_level_smooth_results <- function(snp_gene_ld_loci_data, bw_data,
                                          stat_column, n_gene_points = 100) {

  # Assume the columns are there...

  map_dfr(unique(snp_gene_ld_loci_data$ld_locus_id),
          function(ld_loci_i) {

            # Get the LD loci data
            ld_loci_data <- snp_gene_ld_loci_data %>%
              filter(ld_locus_id == ld_loci_i)

            # Get the LD loci's bandwidth for all:
            ld_loci_bw <- bw_data %>%
              filter(ld_locus_id == ld_loci_i) %>%
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
                                     smooth_zsquared = gene_data[[stat_column]]) %>%
                              bind_rows(tibble(
                                bp = c(min(gene_data$bp) - 1, max(gene_data$bp) + 1),
                                smooth_zsquared = rep(0, 2)
                              ))
                          }

                          smooth_gene_results %>%
                            mutate(ensembl_id = gene_i) %>%
                            return()

                        })

            }

            smooth_ld_loci_gene_results %>%
              mutate(ld_locus_id = ld_loci_i) %>%
              return()

          })
}


# Cluster SNPs ------------------------------------------------------------


#' Cluster SNPs within a locus
#'
#' @param bp Vector of base pair positions
#' @param min_locus_size Minimum locus size required to cluster
#' @param cutoff_prop Proportion of locus size (based on start and end) to control
#' the resulting intra locus clusters
#' @return Vector of cluster assignments using single linkage

get_intra_loci_groups <- function(bp, min_locus_size, cutoff_prop = 0.05) {

  locus_size <- max(bp) - min(bp)

  # If this exceeds the min_locus_size then use the clustering

  if (locus_size >= min_locus_size) {

    # Perform hclust
    intra_locus_hclust <- hclust(dist(bp), method = "single")
    # Get the label with the cutoff_prop:
    results <- cutree(intra_locus_hclust,
                      h = ceiling(locus_size * cutoff_prop))

  } else {
    results <- rep(1, length(bp))
  }

  return(results)

}



# Background  -------------------------------------------------------------

# First write a helper function to generate background values for SCZ and EA
get_locus_cluster_level_scz_ea_background <-
  function(snp_ld_loci_data, bw_data, n_cluster_points = 1000) {

    # First compute the versions of the background SCZ and EA so that they are on
    # the same scale as ASD:

    snp_ld_loci_data <- snp_ld_loci_data %>%
      # First compute the LD loci level test-stats only using the positional SNPs
      group_by(ld_locus_id) %>%
      mutate(asd_test_stat = sum(asd_z_squared),
             scz_test_stat = sum(scz_z_squared),
             ea_test_stat = sum(ea_z_squared)) %>%
      ungroup() %>%
      # Next compute the adjustment factors for scaling the phenotypes
      mutate(scz_adj_factor = asd_test_stat / scz_test_stat,
             ea_adj_factor = asd_test_stat / ea_test_stat,
             # Finally compute the background values
             adj_scz_z_squared = scz_z_squared * scz_adj_factor,
             adj_ea_z_squared = ea_z_squared * ea_adj_factor)

    # Assume the columns are there...
    map_dfr(unique(snp_ld_loci_data$ld_locus_id),
            function(ld_loci_i) {

              # Get the LD loci data:
              ld_loci_data <- snp_ld_loci_data %>%
                filter(ld_locus_id == ld_loci_i)
              ld_loci_bw <- bw_data %>%
                filter(ld_locus_id == ld_loci_i) %>%
                pull(gcv_bw)

              # Now iterate through each cluster:
              map_dfr(unique(ld_loci_data$intra_locus_cluster),
                      function(cluster_i) {

                        # Get the cluster data and bandwidth
                        cluster_data <- ld_loci_data %>%
                          filter(intra_locus_cluster == cluster_i)

                        # If there is a bandwidth value then perform the interpolation
                        # otherwise just return the observed values (this is so that
                        # none of the LD locis are removed from the data for display
                        # in the app)
                        if (is.na(ld_loci_bw) | nrow(cluster_data) < 3) {
                          # Will make these appear as peaks - by adding fake SNPs
                          # to both ends with values of 0
                          smooth_results <- tibble(bp = cluster_data$bp,
                                                   smooth_scz_background = cluster_data[["adj_scz_z_squared"]],
                                                   smooth_ea_background = cluster_data[["adj_ea_z_squared"]]) %>%
                            bind_rows(tibble(
                              bp = c(min(cluster_data$bp) - 1, max(cluster_data$bp) + 1),
                              smooth_scz_background = rep(0, 2),
                              smooth_ea_background = rep(0, 2)
                            ))
                        } else {

                          scz_smooth_results <-
                            get_smooth_interpolation(cluster_data$bp,
                                                     cluster_data[["adj_scz_z_squared"]],
                                                     ld_loci_bw, n_cluster_points) %>%
                            dplyr::rename(smooth_scz_background = smooth_zsquared)
                          ea_smooth_results <-
                            get_smooth_interpolation(cluster_data$bp,
                                                     cluster_data[["adj_ea_z_squared"]],
                                                     ld_loci_bw, n_cluster_points) %>%
                            dplyr::rename(smooth_ea_background = smooth_zsquared)

                          smooth_results <- scz_smooth_results %>%
                            dplyr::left_join(ea_smooth_results, by = "bp")

                        }
                        smooth_results  %>%
                          mutate(intra_locus_cluster = cluster_i) %>%
                          return()
                      }) %>%
                mutate(ld_locus_id = ld_loci_i)
            })
  }


# Write a function to generate null LD loci level results ----------------

get_ld_loci_level_null_results <- function(snp_ld_loci_data, bw_data,
                                           n_cluster_points = 1000, n_null_sims = 1000,
                                           sim_percentiles = c(0.50, 0.75, 0.95),
                                           n_cores = 5) {


  # Generate results over the unique CHR:
  lapply(unique(snp_ld_loci_data$chr),
         function(chr_i) {
           chr_ld_loci_data <- snp_ld_loci_data %>%
             filter(chr == chr_i)

           # Load the chromosome genotype matrix with only these SNPs:
           chr_ref_genotypes <- as.matrix(
             data.table::fread(
               paste0("clean_data/ld_reference/asd/chr_level/chr",
                      chr_i, "_ref_genotypes.csv")))[, unique(chr_ld_loci_data$hg19_id)]

           # Now proceed to generate the results in parallel for each LD loci:
           ld_locus_ids <- unique(chr_ld_loci_data$ld_locus_id)

           chr_ld_loci_summary <-
             mclapply(ld_locus_ids, mc.cores = n_cores,
                      function(ld_loci_i) {
                        #browser()

                        # Get the LD loci data:
                        ld_loci_data <- chr_ld_loci_data %>%
                          filter(ld_locus_id == ld_loci_i)

                        # Get the LD loci's bandwidth:
                        ld_loci_bw <- bw_data %>%
                          filter(ld_locus_id == ld_loci_i) %>%
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
                                                     chr_ref_genotypes[, ld_loci_data$hg19_id]
                                                   ))), method = "chol"))

                          # Generate the smooth results for each sim & cluster:
                          ld_loci_null_results <-
                            map_dfr(unique(ld_loci_data$intra_locus_cluster),
                                    function(cluster_i) {

                                      # Get the cluster data and bandwidth
                                      cluster_data <- ld_loci_data %>%
                                        filter(intra_locus_cluster == cluster_i)

                                      # Get the cluster SNPs' indices:
                                      snp_i <- which(ld_loci_data$hg19_id %in%
                                                       cluster_data$hg19_id)

                                      map_dfr(1:ncol(ld_loci_null_sims),
                                              function(sim_i) {

                                                if (nrow(cluster_data) < 3) {

                                                  tibble(smooth_zsquared =
                                                           (ld_loci_null_sims[snp_i, sim_i])^2,
                                                         bp = cluster_data$bp) %>%
                                                    mutate(sim_index = sim_i)


                                                } else {

                                                  get_smooth_interpolation(
                                                    cluster_data$bp,
                                                    # Make sure to square the z-stats!
                                                    (ld_loci_null_sims[snp_i, sim_i])^2,
                                                    ld_loci_bw, n_cluster_points) %>%
                                                    mutate(sim_index = sim_i)


                                                }

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
                                        mutate(intra_locus_cluster = cluster_i)
                                    }) %>%
                            mutate(ld_locus_id = ld_loci_i)

                        } else {
                          # Otherwise do not set up positional null results
                          ld_loci_null_results <- NA
                        }

                        return(list("null_results" = ld_loci_null_results))

                      })
           return(chr_ld_loci_summary)


         }) %>%
    # Then finally return the non-missing ld_loci
    map_dfr(function(chr_null_list) {
      result_list <-
        map(1:length(chr_null_list),
            function(ld_loci_i) {
              chr_null_list[[ld_loci_i]]$null_results
            })
      bind_rows(result_list[!is.na(result_list)])
    })
}

