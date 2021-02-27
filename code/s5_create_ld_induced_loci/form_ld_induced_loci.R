# PURPOSE: Initialize R functions for performing agglomerative clustering of
#          genes based on induced correlations using quadratic test statistic


# Top level function for creating LD induced loci -------------------------

#' Perform agglomerative construction of LD-induced loci
#'
#' Return a list of LD loci constructed based on agglomerative clustering using
#' their induced correlation for the Gaussian quadratic test statistic, given the
#' SNP to gene assignments, genotype data for calculating SNP correlations, a
#' window size for the genes, and correlation threshold.
#'
#' @param snp_gene_data A tibble where each row is a SNP-gene pair assumed to have
#' information regarding the SNP and gene ids along with the chromosome and
#' gene start/end positions. (NOTE: will assume certain column names for now)
#' @param snp_genotype_matrix A matrix of the SNP genotype data to use for
#' computing the correlation matrix (NOTE: assume column names match SNP ids)
#' @param window_size Window size (in Mb) to use for padding to determine the
#' which genes to compute the correlation matrix with.
#' @param rho_threshold Correlation threshold for deciding to merge genes into
#' same LD-based loci.
#' @return A list of LD induced loci following the agglomerative clustering procedure

form_ld_induced_loci <- function(snp_gene_data, snp_genotype_matrix,
                                 window_size = 6, rho_threshold = 0.25) {
  
  # ASSUME snp_gene_data CONTAINS THE FOLLOWING COLUMNS:
  # - snp_id
  # - gene_id
  # (all genes are in same chr for simplicity)
  # - start
  # - end
  # ASSUME THE COLNAMES OF snp_genotype_matrix ARE THE SNP IDS IN snp_gene_data
  
  # Create a version of the data dropping SNP ids and sort genes by end, creating
  # an index column to use for the gene order:
  sorted_gene_data <- snp_gene_data %>%
    dplyr::select(gene_id, start, end) %>%
    dplyr::distinct() %>%
    dplyr::arrange(end) %>%
    dplyr::mutate(gene_order = 1:dplyr::n()) %>%
    # Create a window column adding the padding to the gene's end position:
    # (note: window_size is assumed to be an integer in Mb while start / end are
    # in bp)
    dplyr::mutate(end_window = end + (window_size * 1e6))
  
  # Extract sorted genes:
  sorted_genes <- dplyr::pull(sorted_gene_data, gene_id)
  
  # Initialize the empty list of gene-sets, vector of genes to merge, and
  # empty vector of window genes:
  genes_to_merge <- sorted_genes
  ld_loci_list <- vector(mode = "list", length = length(sorted_genes))
  window_genes <- c()
  window_snp_cor_matrix <- matrix(nrow = 0, ncol = 0)
  
  # Create a correlation merge log:
  cor_merge_log <- vector(mode = "list", length = length(sorted_genes))
  
  # Now proceed through each of the genes in order:
  for (g in 1:length(sorted_genes)) {
    
    # Whats the current gene:
    gene_g <- sorted_genes[g]
    
    # Make old copy of window genes:
    old_window_genes <- window_genes
    
    # Update the window genes based on the current gene:
    gene_order_g <- sorted_gene_data$gene_order[g]
    window_end_g <- sorted_gene_data$end_window[g]
    window_genes <- sorted_gene_data %>%
      dplyr::filter(gene_order >= gene_order_g,
                    (start <= window_end_g) | (end <= window_end_g)) %>%
      pull(gene_id)
    
    # What are the new genes (excluding gene g as well)
    new_genes <- dplyr::setdiff(window_genes,
                                unique(c(old_window_genes, gene_g)))
    
    # Next pick the gene-set to update along with the candidate genes to consider
    # Check the LD loci membership:
    gene_g_loci_ind <- sapply(ld_loci_list, function(x) gene_g %in% x)
    
    if (any(gene_g_loci_ind)) {
      # Which set has gene g?
      ld_loci_i <- which(gene_g_loci_ind)
      # Grab the genes from that set:
      ld_loci_g <- ld_loci_list[[ld_loci_i]]
      
      # Only use the new genes for the candidates since if it was assigned to
      # a previous set then that set already seen the window genes:
      candidate_genes <- new_genes
      
    } else {
      # Otherwise only use the current gene:
      ld_loci_g <- gene_g
      # Determine the gene-set index:
      ld_loci_i <- ifelse(g == 1, 1,
                           # Otherwise use the first empty gene-set:
                           which(sapply(ld_loci_list, length) == 0)[1])
      # Then everything in the window (excluding gene g) is a candidate:
      candidate_genes <- dplyr::setdiff(window_genes, gene_g)
    }
    
    
    # If there are no candidate_genes genes then remove the current gene
    # from the window and advance to the next gene:
    if (length(candidate_genes) == 0) {
      window_genes <- dplyr::setdiff(window_genes, gene_g)
      next
    }
    
    # Grab the new SNPs
    new_snps <- return_ld_loci_snps(c(ld_loci_g, candidate_genes),
                                     snp_gene_data)
    
    if (any(!(new_snps %in%
              colnames(window_snp_cor_matrix)))) {

      # Next compute the correlation matrix for these SNPs (using a helper
      # function that returns the correlation matrix with names):
      window_snp_cor_matrix <-
        compute_snp_cor_matrix(new_snps,
                               snp_genotype_matrix)
      
    }
    
    # Now proceed to loop through the new genes only merging with the current
    # LD loci if the induced correlation exceeds the given threshold:
    for (new_gene_g in candidate_genes) {
      
      # Make an indicator to handle merging gene-sets:
      is_merge_ld_loci <- FALSE
      
      # First check to see if a candidate gene is already in a set:
      new_gene_g_loci_ind <- sapply(ld_loci_list, function(x) new_gene_g %in% x)
      
      if (any(new_gene_g_loci_ind)) {
        
        is_merge_ld_loci <- TRUE
        
        # Which set has new gene g?
        new_ld_loci_i <- which(new_gene_g_loci_ind)
        
        # If this is the same set then move to the next candidate since it's
        # already joined:
        if (new_ld_loci_i == ld_loci_i) {
          next
        }
        
        # Grab the genes from that set:
        new_ld_loci_g <- ld_loci_list[[new_ld_loci_i]]
        new_total_snps <- return_ld_loci_snps(c(ld_loci_g,
                                                 new_ld_loci_g,
                                                 candidate_genes),
                                               snp_gene_data)
        
        if (any(!(new_total_snps %in%
                  colnames(window_snp_cor_matrix)))) {
          
          # Next compute the correlation matrix for these SNPs (using a helper
          # function that returns the correlation matrix with names):
          window_snp_cor_matrix <-
            compute_snp_cor_matrix(new_total_snps,
                                   snp_genotype_matrix)
          
        }
        
      } else {
        # Otherwise only use the current gene:
        new_ld_loci_g <- new_gene_g
        
      }

      # Compute induced correlation between current LD loci and the new gene g:
      candidate_gene_cor <-
        compute_induced_cor(return_ld_loci_snps(ld_loci_g, snp_gene_data),
                            return_ld_loci_snps(new_ld_loci_g, snp_gene_data),
                            window_snp_cor_matrix)
      
      cor_merge_log[[g]] <- c(cor_merge_log[[g]], candidate_gene_cor)
      
      if (candidate_gene_cor >= rho_threshold) {

        ld_loci_g <- dplyr::union(ld_loci_g, new_ld_loci_g)
        # If merging sets then make the old one empty:
        if (is_merge_ld_loci) {
          ld_loci_list[[new_ld_loci_i]] <- NULL
        }
      }
      
    }
    
    # Update the LD loci list:
    ld_loci_list[[ld_loci_i]] <-
      dplyr::union(ld_loci_list[[ld_loci_i]], ld_loci_g)
    
    # Remove the gene from the window:
    window_genes <- dplyr::setdiff(window_genes, gene_g)
    
  }
  
  # Remove empty gene-sets in list:
  ld_loci_list <- purrr::compact(ld_loci_list)
  cor_merge_log <- purrr::compact(cor_merge_log)
  
  # Now determine which genes were not assigned to a gene-set:
  genes_to_merge <- dplyr::setdiff(genes_to_merge,
                                   unlist(ld_loci_list))
  
  # If there are any remaining need to each individual gene to their own set:
  if (length(genes_to_merge) > 0) {
    ld_loci_list <- c(ld_loci_list, as.list(genes_to_merge))
  }
  
  # Return the LD loci list:
  return(list("ld_loci_list" = ld_loci_list,
              "ld_loci_cor_merge_log" = cor_merge_log))
  
}



# Initialize the lower level helper functions -----------------------------

#' Return SNPs assigned to vector of genes
#'
#' @param genes Vector of genes to returns SNPs for
#' @param snp_gene_data A tibble where each row is a SNP-gene pair assumed to have
#' information regarding the SNP and gene ids along with the chromosome and
#' gene start/end positions. (NOTE: will assume certain column names for now)
#' @return Vector of SNPs assigned to genes.

return_ld_loci_snps <- function(genes, snp_gene_data) {
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
#' @param g1_snps Vector of SNPs for LD loci 1
#' @param g2_snps Vector of SNPs for LD loci 2
#' @param snp_cor_matrix SNP correlation matrix assumed to have the SNP names
#' @return Value of the induced correlation between the two sets of SNPs

compute_induced_cor <- function(g1_snps, g2_snps, snp_cor_matrix) {
  compute_induced_cov(g1_snps, g2_snps, snp_cor_matrix) /
    (sqrt(compute_gene_var(g1_snps, snp_cor_matrix)) *
       sqrt(compute_gene_var(g2_snps, snp_cor_matrix)))
}

#' Compute induced covariance between two sets of SNPs
#' @param g1_snps Vector of SNPs for LD loci 1
#' @param g2_snps Vector of SNPs for LD loci 2
#' @param snp_cor_matrix SNP correlation matrix assumed to have the SNP names
#' @return Value of the induced covariance between the two sets of SNPs

compute_induced_cov <- function(g1_snps, g2_snps, snp_cor_matrix) {
  sum(2 * c(snp_cor_matrix[g1_snps, g2_snps])^2)
}

#' Compute gene variance based on vector of SNPs
#' @param gene_snps Vector of SNPs for LD loci 1
#' @param snp_cor_matrix SNP correlation matrix assumed to have the SNP names
#' @return Value of the gene variance

compute_gene_var <- function(gene_snps, snp_cor_matrix) {
  
  var_sum <- sum(
    2 * (as.vector(
      snp_cor_matrix[gene_snps, gene_snps][
        upper.tri(snp_cor_matrix[gene_snps, gene_snps])])^2)
  )
  
  return(2 * length(gene_snps) + 2 * var_sum)
}






