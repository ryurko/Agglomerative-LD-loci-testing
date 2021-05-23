# PURPOSE: Create tables of results with the SNPs and genes in the selected
#          LD loci across phenotypes and SNP-to-gene assignment

library(tidyverse)

# Write helper functions for loading data ---------------------------------

# First the high-level function for saving the table of genes and SNPs from
# the AdaPT selected LD loci:
save_adapt_snps_gene_tables <-
  function(assign_type, phenotype, rsquared, target_alpha = 0.05) {
    # First load the AdaPT model dataset:
    adapt_model_data <-
      read_csv(paste0("data/adapt_model_data/",
                      assign_type, "/rsquared", rsquared,
                      "_model_data.csv"))

    # Next the AdaPT results:
    adapt_results <-
      readRDS(
        paste0("data/adapt_results/", assign_type,
               "/", phenotype ,"/rsquared", rsquared, "_xgb_cv.rds")
      )

    # Then the gene to LD loci table:
    gene_ld_loci_table <-
      read_csv(
        paste0("data/tidy_gene_locus/",
               assign_type, "/agglom_rsquared", rsquared, "_ld_loci.csv")
      )

    # And the SNP to gene table:
    snp_ld_loci_table <-
      read_csv(
        paste0("data/tidy_snp_gene/", assign_type, ".csv")
      )

    # Get the LD loci for the AdaPT results at the target alpha:
    adapt_ld_loci <- adapt_model_data[which(adapt_results$qvals <= target_alpha),] %>%
      dplyr::select(ld_locus_id) %>%
      pull()

    # Return the two tables of results
    adapt_gene_ld_loci_table <- gene_ld_loci_table %>%
      filter(ld_locus_id %in% adapt_ld_loci)

    adapt_snp_gene_table <- snp_ld_loci_table %>%
      filter(gene_id %in% adapt_gene_ld_loci_table$ensembl_id) %>%
      dplyr::rename(ensembl_id = gene_id)

    # Save both:
    write_csv(adapt_gene_ld_loci_table,
              paste0("data/adapt_results/",
                     assign_type, "/", phenotype, "/snp_gene_tables/",
                     "rsquared", rsquared, "_gene_locus_table.csv"))
    write_csv(adapt_snp_gene_table,
              paste0("data/adapt_results/",
                     assign_type, "/", phenotype, "/snp_gene_tables/",
                     "rsquared", rsquared, "_snp_gene_table.csv"))

    return(paste0("Finished saving results for ",
                  assign_type, " assignment using r^2 = 0.", rsquared))
  }



# Save the results for each combination -----------------------------------

# Positional + eSNPs
save_adapt_snps_gene_tables("positional_esnps", "asd", 25)
save_adapt_snps_gene_tables("positional_esnps", "scz", 25)
save_adapt_snps_gene_tables("positional_esnps", "ea", 25)

save_adapt_snps_gene_tables("positional_esnps", "asd", 50)
save_adapt_snps_gene_tables("positional_esnps", "scz", 50)
save_adapt_snps_gene_tables("positional_esnps", "ea", 50)

save_adapt_snps_gene_tables("positional_esnps", "asd", 75)
save_adapt_snps_gene_tables("positional_esnps", "scz", 75)
save_adapt_snps_gene_tables("positional_esnps", "ea", 75)

# Positional
save_adapt_snps_gene_tables("positional", "asd", 25)
save_adapt_snps_gene_tables("positional", "scz", 25)
save_adapt_snps_gene_tables("positional", "ea", 25)

save_adapt_snps_gene_tables("positional", "asd", 50)
save_adapt_snps_gene_tables("positional", "scz", 50)
save_adapt_snps_gene_tables("positional", "ea", 50)

save_adapt_snps_gene_tables("positional", "asd", 75)
save_adapt_snps_gene_tables("positional", "scz", 75)
save_adapt_snps_gene_tables("positional", "ea", 75)
