# PURPOSE: Compare the gene size distributions between selection and null results.
#          Then initialize the matching sets for FUMA input to adjust for gene size.

library(tidyverse)

# Load the tidy SNP gene and GENCODE tables -------------------------------

# First the tidy SNP-to-gene tables which will be used for computing gene-size
# in terms of the number of assigned SNPs instead of just length
tidy_snp_gene_data <-
  read_csv("data/tidy_snp_gene/positional_esnps.csv")

# Next the gene to LD loci data, to exclude any genes in selected LD loci:
tidy_gene_locus_data <-
  read_csv("data/tidy_gene_locus/positional_esnps/agglom_rsquared25_ld_loci.csv")

# Load the results for the selected LD loci -----------------------------

# ASD
asd_gene_locus_results <-
  read_csv("data/adapt_results/positional_esnps/asd/snp_gene_tables/rsquared25_gene_locus_table.csv")

# SCZ
scz_gene_locus_results <-
  read_csv("data/adapt_results/positional_esnps/scz/snp_gene_tables/rsquared25_gene_locus_table.csv")

# EA
ea_gene_locus_results <-
  read_csv("data/adapt_results/positional_esnps/ea/snp_gene_tables/rsquared25_gene_locus_table.csv")

# Load signal genes -------------------------------------------------------

asd_signals_genes <-
  read_lines("data/adapt_results/positional_esnps/gene_lists/signals/asd_rsquared25_disc_genes.txt")
scz_signals_genes <-
  read_lines("data/adapt_results/positional_esnps/gene_lists/signals/scz_rsquared25_disc_genes.txt")
ea_signals_genes <-
  read_lines("data/adapt_results/positional_esnps/gene_lists/signals/ea_rsquared25_disc_genes.txt")

# Load GENCODE data -------------------------------------------------------

gencode_data <- read_csv("data/gencode/gencode_v21_table.csv")


# Construct tables of sizes by length and SNPs  ---------------------------

# First count the number of SNPs assigned to each gene:
gene_info_data <- tidy_snp_gene_data %>%
  group_by(ensembl_id) %>%
  summarize(n_snps = n()) %>%
  ungroup() %>%
  # join the LD loci id
  dplyr::left_join(tidy_gene_locus_data, by = "ensembl_id") %>%
  # join the gencode info
  dplyr::left_join(gencode_data, by = "ensembl_id") %>%
  # Compute size based on end - start
  mutate(gene_size = end - start,
         # Make an indicator denoting the LD loci selection
         is_asd_locus = ifelse(ld_locus_id %in% asd_gene_locus_results$ld_locus_id,
                                 "Yes", "No"),
         is_scz_locus = ifelse(ld_locus_id %in% scz_gene_locus_results$ld_locus_id,
                                 "Yes", "No"),
         is_ea_locus = ifelse(ld_locus_id %in% ea_gene_locus_results$ld_locus_id,
                                "Yes", "No"),
         # Indicators for signal genes
         is_asd_sig = ifelse(ensembl_id %in% asd_signals_genes, "Yes", "No"),
         is_scz_sig = ifelse(ensembl_id %in% scz_signals_genes, "Yes", "No"),
         is_ea_sig = ifelse(ensembl_id %in% ea_signals_genes, "Yes", "No"))

# Convert Ensembl to Entrez ids -------------------------------------------

unique_entrez_ids <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                           keys = gene_info_data$ensembl_id,
                                           column = "ENTREZID",
                                           keytype = "ENSEMBL",
                                           # Only use ids with a unique identifier
                                           multiVals = "asNA")
gene_info_data$entrez_id <- unique_entrez_ids

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

# Perform initial KS tests ------------------------------------------------

# First for ASD:
asd_ks_test_result <-
  ks.test(pull(filter(gene_info_data, is_asd_locus == "No"), n_snps),
          pull(filter(gene_info_data, is_asd_locus == "Yes"), n_snps),
          alternative = "greater")
# For SCZ
scz_ks_test_result <-
  ks.test(pull(filter(gene_info_data, is_scz_locus == "No"), n_snps),
          pull(filter(gene_info_data, is_scz_locus == "Yes"), n_snps),
          alternative = "greater")

# And EA
ea_ks_test_result <-
  ks.test(pull(filter(gene_info_data, is_ea_locus == "No"), n_snps),
          pull(filter(gene_info_data, is_ea_locus == "Yes"), n_snps),
          alternative = "greater")


# Next with only the signal genes:
asd_signal_ks_test_result <-
  ks.test(pull(filter(gene_info_data, is_asd_locus == "No"), n_snps),
          pull(filter(gene_info_data, is_asd_sig == "Yes"), n_snps),
          alternative = "greater")
# D^+ = 0.14304, p-value = 4.842e-10
# alternative hypothesis: the CDF of x lies above that of y

# For SCZ
scz_signal_ks_test_result <-
  ks.test(pull(filter(gene_info_data, is_scz_locus == "No"), n_snps),
          pull(filter(gene_info_data, is_scz_sig == "Yes"), n_snps),
          alternative = "greater")
# D^+ = 0.17718, p-value < 2.2e-16

# For EA
ea_signal_ks_test_result <-
  ks.test(pull(filter(gene_info_data, is_ea_locus == "No"), n_snps),
          pull(filter(gene_info_data, is_ea_sig == "Yes"), n_snps),
          alternative = "greater")
# D^+ = 0.17655, p-value < 2.2e-16

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
  filter(is_asd_locus == "No", !is.na(entrez_id))
# Repeat for SCZ:
scz_signal_gene_info <- gene_info_data %>%
  filter(is_scz_sig == "Yes")
scz_null_entrez_gene_info <- gene_info_data %>%
  filter(is_scz_locus == "No", !is.na(entrez_id))
# And for SCZ:
ea_signal_gene_info <- gene_info_data %>%
  filter(is_ea_sig == "Yes")
ea_null_entrez_gene_info <- gene_info_data %>%
  filter(is_ea_locus == "No", !is.na(entrez_id))


# Save the matching genes with the signals for FUMA -----------------------

# First make a helper function
save_fuma_gene_list <- function(gene_ids, assign_type, phenotype,
                                is_null = FALSE, is_matched = TRUE,
                                n_match = NULL, null_type = "length") {
  # Next save the text file of just the genes to pass as input:
  gene_file <-
    file(
      paste0("data/fuma/input/",
             assign_type, "/",
             phenotype, "_",
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


# ASD matching ------------------------------------------------------------


asd_matching_null_genes_n_snps  <-
  find_matching_size_null_genes(asd_signal_gene_info$ensembl_id,
                                asd_signal_gene_info$n_snps,
                                asd_null_entrez_gene_info$ensembl_id,
                                asd_null_entrez_gene_info$n_snps,
                                n_matches = 20)
save_fuma_gene_list(c(asd_signal_gene_info$ensembl_id,
                      asd_matching_null_genes_n_snps),
                    "positional_esnps", "asd", is_null = TRUE,
                    is_matched = TRUE, null_type = "n_snps", n_match = 20)


# SCZ matching ------------------------------------------------------------


# For SCZ use 2 matching nulls, which have Entrez IDs
scz_matching_null_genes_n_snps  <-
  find_matching_size_null_genes(scz_signal_gene_info$ensembl_id,
                                scz_signal_gene_info$n_snps,
                                scz_null_entrez_gene_info$ensembl_id,
                                scz_null_entrez_gene_info$n_snps,
                                n_matches = 2)
save_fuma_gene_list(c(scz_signal_gene_info$ensembl_id, scz_matching_null_genes_n_snps),
                    "positional_esnps", "scz", is_null = TRUE,
                    is_matched = TRUE, null_type = "n_snps", n_match = 2)


# EA matching -------------------------------------------------------------


# For EA use 1 matching null, which have Entrez IDs
ea_matching_null_genes_n_snps  <-
  find_matching_size_null_genes(ea_signal_gene_info$ensembl_id,
                                ea_signal_gene_info$n_snps,
                                ea_null_entrez_gene_info$ensembl_id,
                                ea_null_entrez_gene_info$n_snps,
                                n_matches = 1)
save_fuma_gene_list(c(ea_signal_gene_info$ensembl_id, ea_matching_null_genes_n_snps),
                    "positional_esnps", "ea", is_null = TRUE,
                    is_matched = TRUE, null_type = "n_snps", n_match = 1)
