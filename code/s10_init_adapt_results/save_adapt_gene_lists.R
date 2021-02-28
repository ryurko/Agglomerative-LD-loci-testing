# PURPOSE: Save the list of genes corresponding to the AdaPT results for the
#          various types of SNP-to-gene assignments. Do this for each phenotype,
#          using all genes in the associated LD loci

library(tidyverse)

# Write a function for saving text file of genes --------------------------

# First the high-level function for saving the table of genes and SNPs from
# the AdaPT selected LD loci:
save_adapt_gene_lists <-
  function(assign_type, phenotype, rsquared) {
    
    # First load the table of genes:
    adapt_gene_ld_loci_table <-
      read_csv(paste0("data/adapt_results/",
                      assign_type, "/snp_gene_tables/",
                      phenotype, "_rsquared", rsquared, "_gene_ld_loci_table.csv"))
    # Next save the text file of just the genes to pass as input:
    adapt_gene_file <-
      file(
        paste0("data/adapt_results/gene_lists/all/",
               assign_type, "/", phenotype, "_rsquared", rsquared,
               "_disc_genes.txt")
      )
    adapt_gene_ld_loci_table %>%
      pull(ensembl_id) %>%
      paste(collapse = "\n") %>%
      writeLines(adapt_gene_file)
    close(adapt_gene_file)
    
    return(paste0("Finished saving gene list for ",
                  assign_type, " assignment using r^2 = 0.", rsquared))
  }


# Save the results for each combination -----------------------------------

# Positional + eSNPs
save_adapt_gene_lists("positional_esnps", "asd", 25)
save_adapt_gene_lists("positional_esnps", "asd", 50)
save_adapt_gene_lists("positional_esnps", "scz", 25)
save_adapt_gene_lists("positional_esnps", "scz", 50)
save_adapt_gene_lists("positional_esnps", "ea", 25)
save_adapt_gene_lists("positional_esnps", "ea", 50)

# Positional
save_adapt_gene_lists("positional", "asd", 25)
save_adapt_gene_lists("positional", "asd", 50)
save_adapt_gene_lists("positional", "scz", 25)
save_adapt_gene_lists("positional", "scz", 50)
save_adapt_gene_lists("positional", "ea", 25)
save_adapt_gene_lists("positional", "ea", 50)

