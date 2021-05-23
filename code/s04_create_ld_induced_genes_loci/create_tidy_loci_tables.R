# PURPOSE: Create the tidy gene to locus, and SNP to gene/locus tables following
#          the merging algorithm

library(tidyverse)

# First load the tidy SNP to gene tables with p-values --------------------

tidy_positional_esnp_gene_data <-
  read_csv("data/tidy_snp_gene/positional_esnps.csv")

tidy_positional_gene_data <-
  read_csv("data/tidy_snp_gene/positional.csv")


# Write helper functions for initializing loci ----------------------------

# The results are saved as lists of the merged genes for the different possible
# r^2 thresholds, first define a function that just takes in a chromosome results
# to make a single table with each threshold split up

convert_loci_list_to_table <- function(merge_results) {

  map_dfr(1:length(merge_results),
          function(threshold_i) {
            threshold_results <- merge_results[[threshold_i]]
            # Update the names:
            names(threshold_results) <- paste0("ld_locus_", 1:length(threshold_results))

            map_dfr(1:length(threshold_results),
                    function(locus_i) {
                      tibble(ensembl_id = threshold_results[[locus_i]],
                             locus_id = names(threshold_results)[locus_i])
                    }) %>%
              mutate(ld_threshold = names(merge_results)[threshold_i])
        }) %>%
    mutate(locus_type = "locus")
}

# Next a function to then add in the remaining genes that were not merged to
# be in their own separate loci, with new ids for such genes - will just assume
# that the input corresponds to results for a single chr and LD threshold.
# Will then just iterate over chr and thresholds separately

create_gene_locus_table <- function(loci_table, gene_table) {

  # Assume gene table has a column of gene ids, find which ones are not in any
  # of the merged loci:
  non_loci_genes <- gene_table %>%
    filter(!(ensembl_id %in% loci_table$ensembl_id)) %>%
    # Add a column with their ld_locus ids based on the existing number
    mutate(locus_id = paste0("ld_locus_",
                             1:n() + length(unique(loci_table$locus_id)))) %>%
    mutate(locus_type = "gene")

  # Stack the two together and return:
  loci_table %>% bind_rows(non_loci_genes) %>%
    return

}


# Create separate tables for each threshold -------------------------------

# Iterate through each chromosome and then LD threshold to create a list of
# gene to locus tables to save. Start with Positional results:

positional_loci_table <-
  list.files("data/merged_ld_loci/positional", full.names = TRUE) %>%
  map_dfr(function(chr_file) {

    # Get the chr number
    chr_i <- chr_file %>%
      stringr::str_remove_all("data/merged_ld_loci/positional/chr|_results\\.rds") %>%
      as.numeric()

    # Convert the chr file to a table and add the chr index
    read_rds(chr_file) %>%
      convert_loci_list_to_table() %>%
      mutate(chr = chr_i)
  })

# Now create a list of 3 tables corresponding to the genes to loci for each
# threshold, including genes that were not merged

positional_genes_loci_tables <-
  map(unique(positional_loci_table$ld_threshold),
      function(ld_threshold_level) {

        # First create the gene table to use:
        tidy_gene_table <- tidy_positional_gene_data %>%
          dplyr::select(ensembl_id, chr) %>%
          distinct()

        # And now the loci table:
        ld_loci_table <- positional_loci_table %>%
          filter(ld_threshold == ld_threshold_level) %>%
          # Drop this column:
          dplyr::select(-ld_threshold)

        # Now create the gene to locus for each chr:
        map_dfr(unique(ld_loci_table$chr),
                function(chr_i) {
                  create_gene_locus_table(dplyr::filter(ld_loci_table,
                                                        chr == chr_i),
                                          dplyr::filter(tidy_gene_table,
                                                        chr == chr_i))
                }) %>%
          mutate(ld_locus_id = paste0("chr", chr, "_",
                                      str_remove(locus_id, "ld_locus_"))) %>%
          dplyr::select(-locus_id)

      })
names(positional_genes_loci_tables) <- unique(positional_loci_table$ld_threshold)

# Repeat for Positional + eSNPs
positional_esnps_loci_table <-
  list.files("data/merged_ld_loci/positional_esnps", full.names = TRUE) %>%
  map_dfr(function(chr_file) {

    # Get the chr number
    chr_i <- chr_file %>%
      stringr::str_remove_all("data/merged_ld_loci/positional_esnps/chr|_results\\.rds") %>%
      as.numeric()

    # Convert the chr file to a table and add the chr index
    read_rds(chr_file) %>%
      convert_loci_list_to_table() %>%
      mutate(chr = chr_i)
  })

positional_esnps_genes_loci_tables <-
  map(unique(positional_esnps_loci_table$ld_threshold),
      function(ld_threshold_level) {

        # First create the gene table to use:
        tidy_gene_table <- tidy_positional_esnp_gene_data %>%
          dplyr::select(ensembl_id, chr) %>%
          distinct()

        # And now the loci table:
        ld_loci_table <- positional_esnps_loci_table %>%
          filter(ld_threshold == ld_threshold_level) %>%
          # Drop this column:
          dplyr::select(-ld_threshold)

        # Now create the gene to locus for each chr:
        map_dfr(unique(ld_loci_table$chr),
                function(chr_i) {
                  create_gene_locus_table(dplyr::filter(ld_loci_table,
                                                        chr == chr_i),
                                          dplyr::filter(tidy_gene_table,
                                                        chr == chr_i))
                }) %>%
          mutate(ld_locus_id = paste0("chr", chr, "_",
                                      str_remove(locus_id, "ld_locus_"))) %>%
          dplyr::select(-locus_id)

      })
names(positional_esnps_genes_loci_tables) <- unique(positional_esnps_loci_table$ld_threshold)


# Save each table of results separately -----------------------------------

rsquared_thresholds <- c(0.75, 0.50, 0.25)
walk(1:length(rsquared_thresholds),
     function(i) {
       write_csv(positional_genes_loci_tables[[i]],
                 paste0("data/tidy_gene_locus/positional/agglom_rsquared",
                        rsquared_thresholds[i] * 100, "_ld_loci.csv"))
       write_csv(positional_esnps_genes_loci_tables[[i]],
                 paste0("data/tidy_gene_locus/positional_esnps/agglom_rsquared",
                        rsquared_thresholds[i] * 100, "_ld_loci.csv"))
     })


# Create the tidy SNP to locus version ------------------------------------

# Helper function to join gene-sets to SNP-gene data:
create_tidy_snp_locus_data <- function(tidy_snp_gene_data, tidy_gene_locus_data) {
  tidy_snp_gene_data %>%
    dplyr::select(-hg38_id, -hg38_chr, -chr) %>%
    dplyr::left_join(tidy_gene_locus_data, by = "ensembl_id") %>%
    dplyr::select(-ensembl_id) %>%
    # Take the distinct to remove the duplicate SNP-gene rows:
    distinct() %>%
    # Add the squared z stats columns:
    mutate(scz_z_squared = (log(scz_or) / scz_se)^2,
           asd_z_squared = (log(asd_or) / asd_se)^2,
           ea_z_squared = (ea_beta / ea_se)^2)
}

# Apply this to threshold and save:
walk(1:length(rsquared_thresholds),
     function(i) {
       write_csv(create_tidy_snp_locus_data(tidy_positional_gene_data,
                                            positional_genes_loci_tables[[i]]),
                 paste0("data/tidy_snp_locus/positional/agglom_rsquared",
                        rsquared_thresholds[i] * 100, "_ld_loci.csv"))
       write_csv(create_tidy_snp_locus_data(tidy_positional_esnp_gene_data,
                                            positional_esnps_genes_loci_tables[[i]]),
                 paste0("data/tidy_snp_locus/positional_esnps/agglom_rsquared",
                        rsquared_thresholds[i] * 100, "_ld_loci.csv"))
     })

