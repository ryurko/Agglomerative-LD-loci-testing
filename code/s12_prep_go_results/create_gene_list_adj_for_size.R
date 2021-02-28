# PURPOSE: Compare the gene size distributions between selection and null results.
#          Then initialize the matching sets for FUMA input to adjust for gene size.

library(tidyverse)

# Load the tidy SNP gene and GENCODE tables -------------------------------

# First the tidy SNP-to-gene tables which will be used for computing gene-size
# in terms of the number of assigned SNPs instead of just length
tidy_snp_gene_data <-
  read_csv("data/tidy_snp_gene/positional_esnps.csv")

# Next the gene to LD loci data, to exclude any genes in selected LD loci:
tidy_gene_ld_loci_data <-
  read_csv("data/tidy_gene_ld_loci/positional_esnps/agglom_rsquared25_ld_loci.csv")

# Load the GENCODE data
gencode_data <- read_csv("data/gencode/gencode_v21_table.csv")

# Load the results for the selected LD loci -----------------------------

# ASD
asd_gene_ld_loci_results <-
  read_csv("data/adapt_results/positional_esnps/snp_gene_tables/asd_rsquared25_gene_ld_loci_table.csv")

# SCZ
scz_gene_ld_loci_results <-
  read_csv("data/adapt_results/positional_esnps/snp_gene_tables/scz_rsquared25_gene_ld_loci_table.csv")


# Load signal genes -------------------------------------------------------

asd_signals_genes <-
  read_lines("data/adapt_results/gene_lists/signals/positional_esnps/asd_rsquared25_disc_genes.txt")
scz_signals_genes <-
  read_lines("data/adapt_results/gene_lists/signals/positional_esnps/scz_rsquared25_disc_genes.txt")


# Construct tables of sizes by length and SNPs  ---------------------------

# First count the number of SNPs assigned to each gene:
gene_info_data <- tidy_snp_gene_data %>%
  group_by(ensembl_id) %>%
  summarize(n_snps = n()) %>%
  ungroup() %>%
  # join the LD loci id
  dplyr::left_join(tidy_gene_ld_loci_data, by = "ensembl_id") %>%
  # join the gencode info
  dplyr::left_join(gencode_data, by = "ensembl_id") %>%
  # Compute size based on end - start
  mutate(gene_size = end - start,
         # Make an indicator denoting the LD loci selection
         is_asd_ld_loci = ifelse(ld_loci_id %in% asd_gene_ld_loci_results$ld_loci_id,
                            "Yes", "No"),
         is_scz_ld_loci = ifelse(ld_loci_id %in% scz_gene_ld_loci_results$ld_loci_id,
                            "Yes", "No"),
         # Indicators for signal genes
         is_asd_sig = ifelse(ensembl_id %in% asd_signals_genes, "Yes", "No"),
         is_scz_sig = ifelse(ensembl_id %in% scz_signals_genes, "Yes", "No"))

# Convert Ensembl to Entrez ids -------------------------------------------

unique_entrez_ids <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                           keys = gene_info_data$ensembl_id,
                                           column = "ENTREZID",
                                           keytype = "ENSEMBL",
                                           # Only use ids with a unique identifier
                                           multiVals = "asNA")
gene_info_data$entrez_id <- unique_entrez_ids



# Load GENCODE data -------------------------------------------------------

gencode_data <-
  read_csv("data/gencode/gencode_v21_table.csv")

# How many are missing that are also protein coding?
gencode_data %>%
  filter(!(ensembl_id %in% gene_info_data$ensembl_id)) %>%
  pull(gene_biotype) %>%
  table()
# antisense        lnc_rna          other protein_coding
#       605           1194          13117            788

# How many of these missing do not have unique Entrez gene IDs:
gencode_data$entrez_id <-
  AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                        keys = gencode_data$ensembl_id,
                        column = "ENTREZID",
                        keytype = "ENSEMBL",
                        # Only use ids with a unique identifier
                        multiVals = "asNA")
gencode_data %>%
  filter(!(ensembl_id %in% gene_info_data$ensembl_id),
         !is.na(entrez_id))
# 2554 have unique IDs according to this.... but there are probably more

# Perform initial KS tests ------------------------------------------------

# First for ASD:
asd_ks_test_result <-
  ks.test(pull(filter(gene_info_data, is_asd_ld_loci == "No"), n_snps),
          pull(filter(gene_info_data, is_asd_ld_loci == "Yes"), n_snps),
          alternative = "greater")
# For SCZ
scz_ks_test_result <-
  ks.test(pull(filter(gene_info_data, is_scz_ld_loci == "No"), n_snps),
          pull(filter(gene_info_data, is_scz_ld_loci == "Yes"), n_snps),
          alternative = "greater")

# Next with only the signal genes:
asd_signal_ks_test_result <-
  ks.test(pull(filter(gene_info_data, is_asd_ld_loci == "No"), n_snps),
          pull(filter(gene_info_data, is_asd_sig == "Yes"), n_snps),
          alternative = "greater")
# D^+ = 0.16477, p-value = 5.537e-12
# alternative hypothesis: the CDF of x lies above that of y

# For SCZ
scz_signal_ks_test_result <-
  ks.test(pull(filter(gene_info_data, is_scz_ld_loci == "No"), n_snps),
          pull(filter(gene_info_data, is_scz_sig == "Yes"), n_snps),
          alternative = "greater")
# D^+ = 0.15574, p-value < 2.2e-16


# Create matching nulls using genes with Entrez IDs -----------------------

library(optmatch)
find_matching_size_null_genes <-
  function(signal_gene_ids, signal_gene_sizes, null_gene_ids,
           null_gene_sizes, n_matches = 10) {
    
    # Make a version of the data with the each possible pair of selected with
    # every null, computing the difference in their size, first as a long
    # dataset to then convert to wide as a distance matrix:
    pairwise_gene_dist_data <-
      map_dfr(1:length(signal_gene_ids),
              function(signal_gene_i) {
                
                # Init a table of the null genes and their sizes:
                tibble(null_gene_id = null_gene_ids,
                       null_gene_size = null_gene_sizes) %>%
                  # Add the signal gene info:
                  mutate(signal_gene_id = signal_gene_ids[signal_gene_i],
                         signal_gene_size = signal_gene_sizes[signal_gene_i],
                         abs_dist_val = abs(null_gene_size - signal_gene_size)) %>%
                  # Drop the gene size columns:
                  dplyr::select(-null_gene_size, -signal_gene_size)
              }) %>%
      # Next make this a wide dataset to use as the distance matrix
      pivot_wider(names_from = null_gene_id,
                  values_from = abs_dist_val)
    
    # Make the matrix version:
    signal_null_gene_dist_matrix <-
      as.matrix(dplyr::select(pairwise_gene_dist_data, -signal_gene_id))
    # Now just add the row names:
    rownames(signal_null_gene_dist_matrix) <- pairwise_gene_dist_data$signal_gene_id

    # perform matches using the closest null genes in size
    setMaxProblemSize(size = Inf)
    null_size_matches <- pairmatch(signal_null_gene_dist_matrix,
                                   controls = n_matches)
    # Remove the unmatched genes
    null_size_matches <- null_size_matches[!is.na(null_size_matches)]
    # Now just return the names of the matched genes by excluding those in the
    # signal gene vector:
    setdiff(names(null_size_matches),
            pairwise_gene_dist_data$signal_gene_id)
    
  }

# Initialize the datasets of the non-implicated and signal genes for ASD
asd_signal_gene_info <- gene_info_data %>%
  filter(is_asd_sig == "Yes")
asd_null_entrez_gene_info <- gene_info_data %>%
  filter(is_asd_ld_loci == "No", !is.na(entrez_id))
# Repeat for SCZ:
scz_signal_gene_info <- gene_info_data %>%
  filter(is_scz_sig == "Yes")
scz_null_entrez_gene_info <- gene_info_data %>%
  filter(is_scz_ld_loci == "No", !is.na(entrez_id))

# For ASD - use 20 matching nulls, which have Entrez IDs
asd_matching_null_genes  <-
  find_matching_size_null_genes(asd_signal_gene_info$ensembl_id,
                                asd_signal_gene_info$n_snps,
                                asd_null_entrez_gene_info$ensembl_id,
                                asd_null_entrez_gene_info$n_snps,
                                n_matches = 20)

# For SCZ use 2 matching nulls, which have Entrez IDs
scz_matching_null_genes  <-
  find_matching_size_null_genes(scz_signal_gene_info$ensembl_id,
                                scz_signal_gene_info$n_snps,
                                scz_null_entrez_gene_info$ensembl_id,
                                scz_null_entrez_gene_info$n_snps,
                                n_matches = 2)


# Save the matching genes with the signals for FUMA -----------------------

# First make a helper function
save_fuma_gene_list <- function(gene_ids, assign_type, list_type, gene_type,
                                phenotype, is_null = FALSE, is_matched = TRUE,
                                n_match = NULL,
                                null_type = "length") {
  # Next save the text file of just the genes to pass as input:
  gene_file <-
    file(
      paste0("data/adapt_results/fuma/input/",
             assign_type, "/", list_type, "/",
             gene_type, "_", phenotype, "_",
             ifelse(is_null,
                    paste0(ifelse(is_matched,
                                  paste0(n_match, "matched_"), ""),
                           paste0("null_", null_type)), "interesting"),
             "_entrez_genes.txt")
    )
  gene_ids %>%
    paste(collapse = "\n") %>%
    writeLines(gene_file)
  close(gene_file)
  
  return(print("Finished saving gene list"))
  
}


# For ASD:
save_fuma_gene_list(c(asd_signal_gene_info$ensembl_id, asd_matching_null_genes),
                    "positional_esnps", "signals", "all", "asd", is_null = TRUE,
                    is_matched = TRUE, null_type = "n_snps", n_match = 20)
# For SCZ:
save_fuma_gene_list(c(scz_signal_gene_info$ensembl_id, scz_matching_null_genes),
                    "positional_esnps", "signals", "all", "scz", is_null = TRUE,
                    is_matched = TRUE, null_type = "n_snps", n_match = 2)

# Save the null lists without matching, just using all nulls
save_fuma_gene_list(c(asd_signal_gene_info$ensembl_id,
                      asd_null_entrez_gene_info$ensembl_id),
                    "positional_esnps", "signals", "all", "asd", is_null = TRUE,
                    is_matched = FALSE, null_type = "n_snps")
# For SCZ:
save_fuma_gene_list(c(scz_signal_gene_info$ensembl_id,
                      scz_null_entrez_gene_info$ensembl_id),
                    "positional_esnps", "signals", "all", "scz", is_null = TRUE,
                    is_matched = FALSE, null_type = "n_snps")
