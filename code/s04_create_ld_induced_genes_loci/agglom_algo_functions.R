# PURPOSE: Initialize the functions for performing the agglomerative algorithm


# Define top function for computing correlation between genes -------------

#' Compute correlation between genes
#'
#' Walk through the provided table of genes corresponding to a single chr and
#' compute the LD induced correlation between each gene given a window size for
#' genes to consider.
#'
#' @param snp_gene_data A tibble where each row is a SNP-gene pair assumed to have
#' information regarding the SNP and gene ids along with the gene start/end positions.
#' (NOTE: will assume certain column names for now)
#' @param snp_genotype_matrix A matrix of the SNP genotype data to use for
#' computing the correlation matrix (NOTE: assume column names match SNP ids)
#' @param window_size Window size (in Mb) to use for padding to determine the
#' which genes to compute correlations between.
#'
#' @return A table of pairwise correlation values between genes

compute_inter_gene_cor <- function(snp_gene_data, snp_genotype_matrix,
                                   window_size = 6) {

  # ASSUME snp_gene_data CONTAINS THE FOLLOWING COLUMNS:
  # - snp_id
  # - gene_id
  # - start
  # - end

  # Create a version of the data dropping SNP ids and sort genes by end, creating
  # an index column to use for the gene order:
  sorted_gene_data <- snp_gene_data %>%
    dplyr::select(gene_id, start, end) %>%
    dplyr::distinct() %>%
    dplyr::arrange(end) %>%
    # Create an index column for the order by end position:
    dplyr::mutate(gene_order = 1:dplyr::n()) %>%
    # Create a window column adding the padding to the gene's end position:
    # (note: window_size is assumed to be an integer in Mb while start / end are in bp)
    dplyr::mutate(end_window = end + (window_size * 1e6))

  # Initialize a table of values to construct:
  gene_cor_table <- tibble()

  # Initialize empty vector of window genes:
  window_genes <- c()

  # Now proceed through each of the genes in order:
  for (g in 1:nrow(sorted_gene_data)) {

    # Whats the current gene:
    gene_g <- sorted_gene_data$gene_id[g]

    # Make old copy of window genes:
    old_window_genes <- window_genes

    # Update the window genes based on the current gene:
    gene_order_g <- sorted_gene_data$gene_order[g]
    window_end_g <- sorted_gene_data$end_window[g]
    window_genes <- sorted_gene_data %>%
      # Conditions for the window - must end after gene g
      dplyr::filter(gene_order >= gene_order_g,
                    # But fall within the window (either its start or end)
                    (start <= window_end_g) | (end <= window_end_g)) %>%
      pull(gene_id)

    # What are the new genes (excluding gene g as well)
    new_genes <- dplyr::setdiff(window_genes,
                                unique(c(gene_g, old_window_genes)))

    # Next if there are any new genes then compute the pairwise correlations for
    # the old window genes with the new genes, or if the old window was empty then
    # just compute the correlations between all
    if (length(new_genes) > 0) {

      # Compute the updated SNP-level correlation matrix:
      window_snp_cor_matrix <-
        compute_snp_cor_matrix(return_gene_set_snps(window_genes, snp_gene_data),
                               snp_genotype_matrix)

      if (length(old_window_genes) == 0) {


        # Create a table with pairwise correlations:
        window_gene_cor_table <-
          map_dfr(1:(length(window_genes) - 1),
                  function(gene_i) {
                    map_dfr((gene_i + 1):length(window_genes),
                            function(gene_j) {
                              tibble(gene_1 = window_genes[gene_i],
                                     gene_2 = window_genes[gene_j],
                                     induced_r2 =
                                       compute_induced_cor(
                                         return_gene_set_snps(window_genes[gene_i],
                                                              snp_gene_data),
                                         return_gene_set_snps(window_genes[gene_j],
                                                              snp_gene_data),
                                         window_snp_cor_matrix)^2)
                            })
                  })

      } else {

        # Otherwise just compute the correlations between the new genes and the old
        # window genes:
        # Create a table with pairwise correlations:
        window_gene_cor_table <-
          map_dfr(old_window_genes,
                  function(old_gene) {
                    map_dfr(new_genes,
                            function(new_gene) {
                              tibble(gene_1 = old_gene,
                                     gene_2 = new_gene,
                                     induced_r2 =
                                       compute_induced_cor(
                                         return_gene_set_snps(old_gene, snp_gene_data),
                                         return_gene_set_snps(new_gene, snp_gene_data),
                                         window_snp_cor_matrix)^2)
                            })
                  })

      }

      # Update the gene_cor table so far:
      gene_cor_table <- gene_cor_table %>%
        bind_rows(window_gene_cor_table)

    }

    # Remove the gene from the window:
    window_genes <- dplyr::setdiff(window_genes, gene_g)

  }

  # Return the gene_cor_table
  return(gene_cor_table)

}



# Define helper functions -------------------------------------------------

#' Return SNPs assigned to vector of genes
#'
#' @param genes Vector of genes to returns SNPs for
#' @param snp_gene_data A tibble where each row is a SNP-gene pair assumed to have
#' information regarding the SNP and gene ids along with the chromosome and
#' gene start/end positions. (NOTE: will assume certain column names for now)
#' @return Vector of SNPs assigned to genes.

return_gene_set_snps <- function(genes, snp_gene_data) {
  snp_gene_data %>%
    dplyr::filter(gene_id %in% genes) %>%
    dplyr::pull(snp_id) %>%
    unique() %>%
    return
}

#' Compute SNP correlation matrix with row and col names
#'
#' @param snps Vector of SNP ids to use for the correlation matrix
#' @param genotype_matrix Matrix of genotype data to use for computing the
#' correlation matrix.
#' @return Correlation matrix for the provided SNPs with row and col names.

compute_snp_cor_matrix <- function(snps, genotype_matrix) {

  if (length(snps) > 1) {
    snp_cor_matrix <-
      snpcombineR::compute_cor_matrix(genotype_matrix[, snps])
  } else {
    snp_cor_matrix <- matrix(1)
  }


  # Add the row and column names:
  rownames(snp_cor_matrix) <- snps
  colnames(snp_cor_matrix) <- snps

  return(snp_cor_matrix)

}

#' Compute induced correlation between two sets of SNPs
#'
#' @param g1_snps Vector of SNPs for gene set 1
#' @param g2_snps Vector of SNPs for gene set 2
#' @param snp_cor_matrix SNP correlation matrix assumed to have the SNP names
#' @return Value of the induced correlation between the two sets of SNPs

compute_induced_cor <- function(g1_snps, g2_snps, snp_cor_matrix) {
  compute_induced_cov(g1_snps, g2_snps, snp_cor_matrix) /
    (sqrt(compute_gene_var(g1_snps, snp_cor_matrix)) *
       sqrt(compute_gene_var(g2_snps, snp_cor_matrix)))
}

#' Compute induced covariance between two sets of SNPs
#' @param g1_snps Vector of SNPs for gene set 1
#' @param g2_snps Vector of SNPs for gene set 2
#' @param snp_cor_matrix SNP correlation matrix assumed to have the SNP names
#' @return Value of the induced covariance between the two sets of SNPs

compute_induced_cov <- function(g1_snps, g2_snps, snp_cor_matrix) {
  2 * sum(c(snp_cor_matrix[g1_snps, g2_snps])^2)
}

#' Compute gene variance based on vector of SNPs
#' @param gene_snps Vector of SNPs for gene set
#' @param snp_cor_matrix SNP correlation matrix assumed to have the SNP names
#' @return Value of the gene variance

compute_gene_var <- function(gene_snps, snp_cor_matrix) {

  var_sum <- sum(
    as.vector(
      snp_cor_matrix[gene_snps, gene_snps][
        upper.tri(snp_cor_matrix[gene_snps, gene_snps])])^2
  )

  return(2 * length(gene_snps) + 4 * var_sum)
}


# Define function for agglomerative merging of genes ----------------------

#' Agglomerative algorithm for merging genes based on LD-induced correlation
#'
#' Given a table of pairwise correlations for genes and corresponding genotype
#' data, repeatedly merge genes into loci based until remaining pairs are all
#' below desired threshold(s).
#'
#' @param snp_gene_data A tibble where each row is a SNP-gene pair assumed to have
#' information regarding the SNP and gene ids along with the gene start/end positions.
#' (NOTE: will assume certain column names for now)
#' @param snp_genotype_matrix A matrix of the SNP genotype data to use for
#' computing the correlation matrix (NOTE: assume column names match SNP ids)
#' @param gene_cor_table Table of pairwise correlations between genes
#' @param ld_thresholds Vector of r^2 values to use as thresholds for returning
#' lists of merged genes in loci.
#' @param window_size Window size (in Mb) to use for padding to determine the
#' which genes to compute correlations between.
#' @return List of merged genes into loci for each value in ld_thresholds

merge_correlated_genes <- function(snp_gene_data, snp_genotype_matrix,
                                   gene_cor_table, ld_thresholds = c(0.25, 0.5, 0.75),
                                   window_size = 6) {

  # Make a new version of the gene_cor_table but in sorted order and renamed
  # columns (assumed to have gene_1, gene_2, and induced_r2 as columns)
  sorted_set_cor_table <- gene_cor_table %>%
    arrange(desc(induced_r2)) %>%
    # Rename the columns:
    rename(set_1 = gene_1, set_2 = gene_2)

  # Create a table of just the gene location information for updating with
  # merged loci information, given the snp_gene_data
  loci_location_table <- snp_gene_data %>%
    dplyr::select(gene_id, start, end) %>%
    dplyr::distinct() %>%
    dplyr::arrange(end) %>%
    dplyr::rename(locus_id = gene_id) %>%
    mutate(start_window = start - (window_size * 1e6),
           end_window = end + (window_size * 1e6))

  # Initialize the list of results corresponding to the provided thresholds
  ld_thresholds <- sort(ld_thresholds, decreasing = TRUE) # Sort in descending order
  result_list <- vector(mode = "list", length = length(ld_thresholds))
  names(result_list) <- paste0("LD threshold = ", ld_thresholds)

  # Now initialize the list to use for the merging:
  ld_locus_list <- vector(mode = "list", length = nrow(loci_location_table))
  names(ld_locus_list) <- paste0("ld_locus_", 1:length(ld_locus_list))

  # Keep track of the index for the current threshold:
  ld_threshold_i <- 1
  # Make a boolean indicator to control the merging:
  keep_merging <- TRUE

  # Start the merging procedure, storing results for different thresholds
  # along the way:
  while (keep_merging) {

    # Get the top row:
    top_cor_pair <- sorted_set_cor_table[1,]

    # Now check to see if this pair is below the smallest threshold, if it is
    # then add to the results list and stop merging:
    if (top_cor_pair$induced_r2 < min(ld_thresholds)) {

      # Assuming the ld_threshold_i has been updated to its appropriate index
      # along the search:
      result_list[[ld_threshold_i]] <- purrr::compact(ld_locus_list)
      keep_merging <- FALSE

      # Otherwise check for the current threshold which may not be the minimum,
      # store its results if it was met but then proceed to merge the pair
      # and updating the correlation table for everything within the newly
      # merged locus' window

    } else {

      # First check if the current threshold was met
      if (top_cor_pair$induced_r2 < ld_thresholds[ld_threshold_i]) {

        # Add the current list to its results:
        result_list[[ld_threshold_i]] <- purrr::compact(ld_locus_list)
        # Update the LD threshold index:
        ld_threshold_i <- ld_threshold_i + 1

      }

      # Now continue to merge and update
      set_pair <- c(top_cor_pair$set_1, top_cor_pair$set_2)
      # Can remove the rows in the correlation table belonging to these two
      # loci since they will need to be updated:
      sorted_set_cor_table <- sorted_set_cor_table %>%
        dplyr::filter(!(set_1 %in% set_pair) & !(set_2 %in% set_pair))

      # Get the new start and end window for joining whatever these loci are
      # together, then remove their rows from the loci_location_table to update
      # later with the locus_id decided:
      new_locus_window <- loci_location_table %>%
        dplyr::filter(locus_id %in% set_pair) %>%
        dplyr::summarize(start = min(start),
                         end = max(end),
                         start_window = min(start_window),
                         end_window = max(end_window))
      loci_location_table <- loci_location_table %>%
        dplyr::filter(!(locus_id %in% set_pair))

      # There are 3 possible cases now:
      # 1) Need to merge two genes into a new locus (which happens first)
      # 2) Need to merge one gene with one locus
      # 3) Need to merge two loci
      pair_locus_i <- which(stringr::str_detect(set_pair, "ld_locus_"))

      # Will go through these cases before updating the cor table

      # Case 1 - the base case, merging two genes into a new locus:
      if (length(pair_locus_i) == 0) {

        # Just use these two genes as the new locus genes:
        new_locus_genes <- set_pair

        # Find the index of the next empty set:
        new_locus_i <- which(sapply(ld_locus_list, length) == 0)[1]
        # This inherits the name from the list to use

        # Assign these genes to this set:
        ld_locus_list[[new_locus_i]] <- new_locus_genes

      } else if (length(pair_locus_i) == 1) {

        # Will just add the single gene to the already existing locus:
        new_locus_i <- which(names(ld_locus_list) == set_pair[pair_locus_i])
        # Get the new locus genes:
        new_locus_genes <- c(set_pair[-pair_locus_i],
                             ld_locus_list[[new_locus_i]])

        # Update this locus' genes:
        ld_locus_list[[new_locus_i]] <- new_locus_genes

      } else {

        # Both are loci to be merged, meaning one loci must be emptied. For
        # ease will just use whichever one is first in the loci list as the
        # one to retain. Start by getting the loci indices:
        list_pair_i <- which(names(ld_locus_list) %in% set_pair)
        # Use the first one as the merge locus:
        new_locus_i <- list_pair_i[1]
        old_locus_i <- list_pair_i[2]

        # Get the new locus genes:
        new_locus_genes <- c(ld_locus_list[[new_locus_i]],
                             ld_locus_list[[old_locus_i]])

        # Update this locus' genes:
        ld_locus_list[[new_locus_i]] <- new_locus_genes
        # Empty the old locus:
        ld_locus_list[[old_locus_i]] <- NULL


      }

      # Now that the list is updated, need to update the correlation table for
      # the new locus and the genes/loci within its window

      # Will separate the correlation computation between start and end to keep
      # it smaller in terms of the number of SNPs
      # Now find the candidate loci within the start window (and middle)
      start_window_loci <- loci_location_table %>%
        dplyr::filter(end >= new_locus_window$start_window,
                      end <= new_locus_window$end) %>%
        dplyr::pull(locus_id)

      # Only perform this operation if the window is not empty:
      if (length(start_window_loci) > 0) {

        # Get the genes for these:
        start_window_genes <- lapply(start_window_loci,
                                     function(x) {
                                       if (stringr::str_detect(x, "ld_locus")) {
                                         ld_locus_list[x]
                                       } else {
                                         x
                                       }
                                     })

        # Compute the correlations for these separately, starting with the start
        # window loci candidates:
        new_locus_window_cor_matrix <-
          compute_snp_cor_matrix(return_gene_set_snps(c(unlist(ld_locus_list[[new_locus_i]]),
                                                        unlist(start_window_genes)),
                                                      snp_gene_data),
                                 snp_genotype_matrix)

        # Create a table with pairwise correlations:
        start_window_cor_table <-
          purrr::map_dfr(1:length(start_window_genes),
                         function(start_locus_i) {

                           tibble(set_1 = names(ld_locus_list)[new_locus_i],
                                  set_2 = start_window_loci[start_locus_i],
                                  induced_r2 =
                                    compute_induced_cor(
                                      return_gene_set_snps(unlist(ld_locus_list[[new_locus_i]]),
                                                           snp_gene_data),
                                      return_gene_set_snps(unlist(start_window_genes[[start_locus_i]]),
                                                           snp_gene_data),
                                      new_locus_window_cor_matrix)^2)

                         })

        # Now add these to the sorted table and re-sort by induced_r2:
        sorted_set_cor_table <- sorted_set_cor_table %>%
          dplyr::bind_rows(start_window_cor_table) %>%
          dplyr::arrange(desc(induced_r2))

      }


      # For the end window:
      end_window_loci <- loci_location_table %>%
        dplyr::filter(start >= new_locus_window$start,
                      start <= new_locus_window$end_window,
                      end > new_locus_window$end) %>%
        dplyr::pull(locus_id)

      if (length(end_window_loci) > 0) {

        # Get the genes for these:
        end_window_genes <- lapply(end_window_loci,
                                   function(x) {
                                     if (stringr::str_detect(x, "ld_locus")) {
                                       ld_locus_list[x]
                                     } else {
                                       x
                                     }
                                   })

        # Repeat for the end window:
        new_locus_window_cor_matrix <-
          compute_snp_cor_matrix(return_gene_set_snps(c(unlist(ld_locus_list[[new_locus_i]]),
                                                        unlist(end_window_genes)),
                                                      snp_gene_data),
                                 snp_genotype_matrix)

        # Create a table with pairwise correlations:
        end_window_cor_table <-
          purrr::map_dfr(1:length(end_window_genes),
                         function(end_locus_i) {

                           tibble(set_1 = names(ld_locus_list)[new_locus_i],
                                  set_2 = end_window_loci[end_locus_i],
                                  induced_r2 =
                                    compute_induced_cor(
                                      return_gene_set_snps(unlist(ld_locus_list[[new_locus_i]]),
                                                           snp_gene_data),
                                      return_gene_set_snps(unlist(end_window_genes[[end_locus_i]]),
                                                           snp_gene_data),
                                      new_locus_window_cor_matrix)^2)

                         })

        # Now add these to the sorted table and re-sort by induced_r2:
        sorted_set_cor_table <- sorted_set_cor_table %>%
          dplyr::bind_rows(end_window_cor_table) %>%
          dplyr::arrange(desc(induced_r2))


      }

      # Add the new locus into the loci location table
      new_locus_window$locus_id <- names(ld_locus_list)[new_locus_i]
      loci_location_table <- loci_location_table %>%
        bind_rows(new_locus_window)

    }

  }

  # Return the final list of results:
  return(result_list)
  # return(list("result_list" = result_list,
  #             "sorted_cor_table" = sorted_set_cor_table,
  #             "loci_location_table" = loci_location_table))

}


