# PURPOSE: Initialize helper functions for generating the kernel smoothing
#          results for the upload shiny app data

#library(np, lib.loc = "/home/ryurko/Rpackages")
library(np)

# Required input for the app will be two tables:

# 1) Gene table containing at minimum the following columns:
#     - gene_id: e.g., Ensembl ID
#     - start
#     - end
#     - chr
#     - loci_id: an identifier for the gene's loci

# 2) SNP-gene table containing at minimum the following columns:
#     - snp_id: e.g., chr:bp-bp - some unique SNP identifier)
#     - gene_id (the gene the SNP is assigned to)
#     - bp: the SNP's genomic position
#     - snp_signal: SNP-level signal column to provide smoothing over
#                   (e.g., squared z-stats or -log10(p-val))

# And finally there will be input for the number of points in the interpolation

# Helper function for tuning bandwidth with GCV ---------------------------

#' Get GCV bandwidth for each loci
#'
#' @param snp_loci_results SNP-loci pair dataset that is assumed to have
#' `loci_id`, `snp_signal`, and `bp` columns
#' @param min_n_snps Minimum number of SNPs for smoothing (default is 3)
#' @return Dataset with each loci's GCV bandwidth

get_gcv_bw <- function(snp_loci_results, min_n_snps = 3) {

  # Assume input dataset has the loci_id, snp_signal, and bp columns:
  map_dfr(unique(snp_loci_results$loci_id),
          function(loci_i) {

            loci_data <- snp_loci_results %>%
              filter(loci_id == loci_i)

            # If the gene-set has at least the min number of SNPs then proceed
            # with GCV to pick bandwidth:
            if (nrow(loci_data) >= min_n_snps) {
              gs_bw0 <- npregbw(snp_signal ~ bp, data = loci_data)$bw

            } else {
              # Otherwise its NA:
              gs_bw0 <- NA
            }

            # Return a tibble with these bandwidths:
            tibble(loci_id = loci_i, gcv_bw = gs_bw0)

          })

}


# Helper functions for performing kernel smoothing w/ interpolation ------

# First the general function that is called for performing the interpolation

#' Perform kernel smoothing interpolation
#'
#' @param snp_data Dataset of SNPs to a particular group for generating the
#' kernel smoothing results for (assumed to have `bp` and `snp_signal`)
#' @param h Kernel bandwidth to use in kernel regression
#' @param n_points Number of points to generate for interpolation
#' @return Dataset with the interpolated grid of points with the smooth signal

get_smooth_interpolation <- function(snp_data, h, n_points) {

  # Create grid points to evaluate at
  bp_grid <- seq(min(snp_data$bp), max(snp_data$bp), length.out = n_points)

  # Round to nearest digit and take unique:
  bp_grid <- unique(round(bp_grid))

  # Combine this with the actual data points:
  bp_grid <- sort(unique(c(bp_grid, snp_data$bp)))

  # Return a tibble with the smooth values along this grid using the
  # actual values to train it:
  tibble(bp = bp_grid,
         smooth_signal = npreg(txdat = snp_data$bp,
                                 tydat = snp_data$snp_signal,
                                 exdat = bp_grid,
                                 bws = h)$mean) %>%
    # Add an indicator denoting if it is a fake BP:
    mutate(is_fake_bp = 1 - as.numeric(bp %in% snp_data$bp))

}


# Next define the function to generate the LD loci level results

#' Generate loci level kernel smoothing results
#'
#' @param snp_loci_data Dataset of SNPs assigned to loci that is
#' assumed to have `loci_id`, `snp_signal`, and `bp` columns
#' @param loci_gcv_data Dataset of loci GCV bandwidths
#' @param n_points Number of points to generate for interpolation (default is 1000)
#' @return Dataset with the interpolated, smooth signal for every loci

get_loci_level_smoothing <- function(snp_loci_data, loci_gcv_data,
                                       n_points = 1000) {

  snp_loci_data %>%
    # Apply to each loci id:
    pull(loci_id) %>%
    unique() %>%
    map_dfr(function(loci_i) {

      # Get the loci's data
      loci_data <- snp_loci_data %>%
        filter(loci_id == loci_i)

      # Get the loci's bandwidth:
      loci_bw <- loci_gcv_data %>%
        filter(loci_id == loci_i) %>%
        pull(gcv_bw)

      if (!is.na(loci_bw)) {
        smooth_results <- get_smooth_interpolation(loci_data, loci_bw,
                                                   n_points)
      } else {
        smooth_results <- tibble(bp = NA, smooth_signal = NA, is_fake_bp = NA)
      }

      smooth_results %>%
        mutate(loci_id = loci_i) %>%
        return()
    })
}


# Next define the gene-level smoothing equivalent which differs in that it
# also requires a minimum number of SNPs for performing the smoothing

#' Generate gene-level kernel smoothing results
#'
#' @param snp_gene_loci_data Dataset of SNPs assigned to genes that is
#' assumed to have `gene_id`, `loci_id`, `snp_signal`, and `bp` columns
#' @param loci_gcv_data Dataset of loci GCV bandwidths
#' @param n_points Number of points to generate for interpolation (default is 100)
#' @param min_n_snps Minimum number of SNPs for smoothing (default is 3)
#' @return Dataset with the interpolated, smooth signal for every gene

get_gene_level_smoothing <- function(snp_gene_loci_data, loci_gcv_data,
                                       n_points = 100, min_n_snps = 3) {

  # Apply to each loci id:
  snp_gene_loci_data %>%
    pull(loci_id) %>%
    unique() %>%
    map_dfr(function(loci_i) {

      # Get the loci's data:
      loci_data <- snp_gene_loci_data %>%
        filter(loci_id == loci_i)

      # Get the loci's bandwidth:
      loci_bw <- loci_gcv_data %>%
        filter(loci_id == loci_i) %>%
        pull(gcv_bw)

      # What are the loci's genes:
      loci_genes <- loci_data %>%
        pull(gene_id) %>%
        unique()

      # If the loci was sufficiently large enough to generate the smoothing
      # for then proceed at the gene-level:
      if (!is.na(loci_bw)) {

        # Proceed through each gene:
        smooth_gene_results <-
          map_dfr(loci_genes,
                  function(gene_i) {
                    gene_data <- loci_data %>%
                      filter(gene_id == gene_i)

                    # Only generate the results if it has the min number SNPs
                    if (nrow(gene_data) >= min_n_snps) {
                      smooth_gene_results <-
                        get_smooth_interpolation(gene_data, loci_bw, n_points)
                    } else {
                      smooth_gene_results <-
                        tibble(bp = NA, smooth_signal = NA, is_fake_bp = NA)
                    }

                    smooth_gene_results %>%
                      mutate(gene_id = gene_i) %>%
                      return()
                  })

      } else {
        smooth_gene_results <- tibble(bp = NA, smooth_signal = NA,
                                      is_fake_bp = NA, gene_id = NA)
      }

      smooth_gene_results %>%
        mutate(loci_id = loci_i) %>%
        return()

    }) %>%
    return()

}

