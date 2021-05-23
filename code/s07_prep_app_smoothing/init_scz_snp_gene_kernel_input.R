# PURPOSE: Initialize input datasets for the kernel smoothing - separating the
#          SNPs that are positionally assigned to genes versus the ones that
#          are assigned based on functional status.

# SCZ VERSION

library(tidyverse)

# Load the LD loci tables ------------------------------------------------

# Starting with r^2 = 0.25 for scz using Positional + eSNPs:
pos_esnps_gene_ld_loci_table <-
  read_csv("data/adapt_results/positional_esnps/scz/snp_gene_tables/rsquared25_gene_locus_table.csv")

# For Positional:
pos_gene_ld_loci_table <-
  read_csv("data/adapt_results/positional/scz/snp_gene_tables/rsquared25_gene_locus_table.csv")


# Load the GENCODE table --------------------------------------------------

gene_type_table <-
  read_csv("data/gencode/gencode_v21_table.csv")

# Load the tidy SNP tables ------------------------------------------------

# Load the SNP to gene with p-values:
tidy_pos_snp_to_gene_table <-
  read_csv("data/tidy_snp_gene/positional.csv")

tidy_esnps_to_gene_table <-
  read_csv("data/tidy_snp_gene/esnps.csv")



# Create positional versus function SNP tables ----------------------------

# For both sets of results, need to create separate SNP tables for both to
# LD loci and then gene assignments, along with separating the SNPs that
# are positionally assigned to genes versus ones that are only functional

# Initialize the pipeline that will be used regardless of the type of result:
add_snp_data_cols <- . %>%
  # Add the squared z stats columns:
  mutate(asd_z_squared = (log(asd_or) / asd_se)^2,
         scz_z_squared = (log(scz_or) / scz_se)^2,
         ea_z_squared = (ea_beta / ea_se)^2) %>%
  separate(hg38_id, into = c("chr38", "bp"), sep = ":", remove = FALSE) %>%
  separate(bp, into = c("bp", "bp2"), sep = "-", remove = TRUE) %>%
  dplyr::select(-chr38, -bp2) %>%
  mutate(bp = as.numeric(bp))

# Pipeline to join gene info:
join_gene_info <- . %>%
  dplyr::left_join(dplyr::select(gene_type_table, ensembl_id, gene_name,
                                 gene_chr, start, end, strand),
                   by = "ensembl_id")


# Create and save Positional results --------------------------------------

# Start with the easy case of the positional results which only include the
# SNPs are assigned to genes based on the positional location
pos_snp_gene_results <- tidy_pos_snp_to_gene_table %>%
  filter(ensembl_id %in% unique(pos_gene_ld_loci_table$ensembl_id)) %>%
  add_snp_data_cols() %>%
  join_gene_info() %>%
  # And join the LD loci assignments:
  left_join(dplyr::select(pos_gene_ld_loci_table, -chr), by = "ensembl_id")

# Save this table with gene info:
write_csv(pos_snp_gene_results,
          "data/kernel_smoothing/scz/input/positional/positional_snp_gene_ld_loci_data.csv")

# Now make a version without the gene assignments:
pos_snp_ld_loci_results <- pos_snp_gene_results %>%
  dplyr::select(-c(ensembl_id, gene_name, gene_chr, start, end, strand)) %>%
  distinct()
# And save:
write_csv(pos_snp_ld_loci_results,
          "data/kernel_smoothing/scz/input/positional/positional_snp_ld_loci_data.csv")

# Next save the positional gene info:
write_csv(filter(gene_type_table,
                 ensembl_id %in% unique(pos_gene_ld_loci_table$ensembl_id)) %>%
            dplyr::left_join(dplyr::select(pos_gene_ld_loci_table, -chr),
                             by = "ensembl_id"),
          "data/kernel_smoothing/scz/input/positional/gene_info.csv")

# Create and save Positional + eSNPs results ------------------------------

# First save the gene info:
write_csv(filter(gene_type_table,
                 ensembl_id %in% unique(pos_esnps_gene_ld_loci_table$ensembl_id)) %>%
            dplyr::left_join(dplyr::select(pos_esnps_gene_ld_loci_table, -chr),
                             by = "ensembl_id"),
          "data/kernel_smoothing/scz/input/positional_esnps/gene_info.csv")


# Next first make a table of the positional only SNPs for the Positional + eSNPs
# results, to then stack with the eSNPs to find which ones are in both
# (will use loc versus fun for location versus functional abbreviation)
pos_esnps_loc_snp_gene_results <- tidy_pos_snp_to_gene_table %>%
  filter(ensembl_id %in% unique(pos_esnps_gene_ld_loci_table$ensembl_id)) %>%
  mutate(is_positional = 1, is_esnp = 0)

# Repeat for the functional eSNPs:
pos_esnps_fun_snp_gene_results <- tidy_esnps_to_gene_table %>%
  filter(ensembl_id %in% unique(pos_esnps_gene_ld_loci_table$ensembl_id)) %>%
  mutate(is_positional = 0, is_esnp = 1)

# Stack to get the assignment types for the SNP gene pairs:
pos_esnps_pairs <- bind_rows(pos_esnps_loc_snp_gene_results,
                             pos_esnps_fun_snp_gene_results) %>%
  dplyr::select(ensembl_id, hg38_id, is_positional, is_esnp) %>%
  unite(gene_snp_id, c("ensembl_id", "hg38_id")) %>%
  group_by(gene_snp_id) %>%
  summarize(is_positional = max(is_positional),
            is_esnp = max(is_esnp)) %>%
  ungroup()

# Now check that all of the positional SNP gene pairs are included already in
# the pos_esnps_loc results:
loc_pairs <- pos_esnps_loc_snp_gene_results %>%
  unite(gene_snp_id, c("ensembl_id", "hg38_id")) %>%
  pull(gene_snp_id) %>%
  unique()
all(loc_pairs %in% pull(filter(pos_esnps_pairs, is_positional == 1),
                        gene_snp_id))
# [1] TRUE

# Good can proceed to save the tables without needing to do any procedure
# other than dropping the is_ columns

# Create and save the positional SNPs with gene info table:
pos_esnps_loc_snp_gene_results <- pos_esnps_loc_snp_gene_results %>%
  dplyr::select(-is_positional, -is_esnp) %>%
  add_snp_data_cols() %>%
  join_gene_info() %>%
  # And join the LD loci assignments:
  left_join(dplyr::select(pos_esnps_gene_ld_loci_table, -chr),
            by = "ensembl_id")

# Save this table with gene info:
write_csv(pos_esnps_loc_snp_gene_results,
          "data/kernel_smoothing/scz/input/positional_esnps/positional_snp_gene_ld_loci_data.csv")

# Now make a version without the gene assignments:
pos_esnps_loc_snp_ld_loci_results <- pos_esnps_loc_snp_gene_results %>%
  dplyr::select(-c(ensembl_id, gene_name, gene_chr, start, end, strand)) %>%
  distinct()
# And save:
write_csv(pos_esnps_loc_snp_ld_loci_results,
          "data/kernel_smoothing/scz/input/positional_esnps/positional_snp_ld_loci_data.csv")

# Finally save the eSNPs only data with the associated gene assignments:
pos_esnps_fun_snp_gene_results <- pos_esnps_fun_snp_gene_results %>%
  dplyr::select(-is_positional, -is_esnp) %>%
  # Join the update SNP type status:
  unite(gene_snp_id, c("ensembl_id", "hg38_id"), remove = FALSE) %>%
  dplyr::left_join(pos_esnps_pairs, by = "gene_snp_id") %>%
  dplyr::select(-gene_snp_id) %>%
  add_snp_data_cols() %>%
  join_gene_info() %>%
  # And join the LD loci assignments:
  left_join(dplyr::select(pos_esnps_gene_ld_loci_table, -chr),
            by = "ensembl_id")
write_csv(pos_esnps_fun_snp_gene_results,
          "data/kernel_smoothing/scz/input/positional_esnps/functional_snp_gene_ld_loci_data.csv")
