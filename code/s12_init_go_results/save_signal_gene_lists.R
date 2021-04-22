# PURPOSE: Save the list of signal genes identified within each of the LD loci

library(tidyverse)
library(data.table)

# Write a function for finding the signal genes -----------------------------

save_signal_gene_lists <-
  function(assign_type, phenotype, rsquared, include_functional = TRUE) {
    
    
    # First load the smooth LD loci results ----
    ld_loci_smooth_results <-
      read_csv(paste0("data/kernel_smoothing/", phenotype, "/output/", 
                      assign_type, "/positional_ld_loci_level_smoothing.csv"))
    
    # Find the signals above the null cutoff:
    if (phenotype == "asd") {
      
      ld_loci_signals_data <- ld_loci_smooth_results %>%
        # Just above the null percentile
        filter(smooth_asd_zsquared >= null_95_percent)
      
    } else if (phenotype == "scz") {
      
      ld_loci_signals_data <- ld_loci_smooth_results %>%
        filter(smooth_scz_zsquared >= null_95_percent)
      
    } else {
      
      ld_loci_signals_data <- ld_loci_smooth_results %>%
        filter(smooth_ea_zsquared >= null_95_percent)
      
    }
    ld_loci_signals_data <- ld_loci_signals_data %>%
      # Add the ld_loci_chr column:
      separate(ld_loci_id, c("ld_loci_chr", "set_id"), sep = "_", remove = FALSE) %>%
      mutate(ld_loci_chr = as.numeric(str_remove(ld_loci_chr, "chr"))) %>%
      dplyr::select(-set_id)
    
    # Load the genes to annotate signals to ----
    gene_info_data <-
      read_csv(paste0("data/kernel_smoothing/", phenotype, "/input/",
                      assign_type, "/gene_info.csv"))
    ld_loci_gene_info_data <- gene_info_data %>%
      # Only interested in the signals LD loci
      filter(ld_loci_id %in% unique(ld_loci_signals_data$ld_loci_id))
    
    # Assign the signals to genes ----
    
    # Initialize the SNP dt object for using foverlaps
    snp_signals_dt <- data.table(
      Chr = ld_loci_signals_data$ld_loci_chr,
      start = ld_loci_signals_data$bp,
      end = ld_loci_signals_data$bp,
      id = ld_loci_signals_data$snp_id,
      key = c("Chr", "start", "end")
    )
    
    # Set-up the gene tables using the gencode data:
    genes_dt <- data.table(
      Chr = ld_loci_gene_info_data$gene_chr,
      start = ld_loci_gene_info_data$start,
      end = ld_loci_gene_info_data$end,
      id = ld_loci_gene_info_data$ensembl_id,
      key = c("Chr", "start", "end")
    )
    
    # Make the overlaps table
    signal_gene_assign <- foverlaps(snp_signals_dt, genes_dt)
    
    # Convert to the tidy SNP-gene dataset:
    tidy_signal_gene_data <- signal_gene_assign %>%
      # Remove SNPs unable to map to genes:
      .[!is.na(id),] %>%
      as_tibble() %>%
      dplyr::select(id, i.id, Chr) %>%
      dplyr::rename(ensembl_id = id,
                    peak_id = i.id,
                    chr = Chr)
    
    # What are the signal genes:
    signal_genes <- tidy_signal_gene_data %>% pull(ensembl_id) %>% unique()
    
    # Additionally grab any LD loci that were not included in the signals
    non_signal_genes <- gene_info_data %>%
      filter(ld_loci_id %in%
               setdiff(gene_info_data$ld_loci_id,
                       ld_loci_signals_data$ld_loci_id)) %>%
      pull(ensembl_id)
    
    # Get genes with functional SNPs with noticeable marginal effects:
    if (assign_type != "positional" & include_functional) {
      functional_gene_data <-
        read_csv(paste0("data/kernel_smoothing/", phenotype, "/input/",
                        assign_type, "/functional_snp_gene_ld_loci_data.csv"))
      
      if (phenotype == "asd") {
        
        functional_gene_data <- functional_gene_data %>%
          filter(asd_z_squared >= qchisq(0.95, 1))
        
      } else if (phenotype == "scz") {
        
        functional_gene_data <- functional_gene_data %>%
          filter(scz_z_squared >= qchisq(0.95, 1))
        
      } else {
        
        functional_gene_data <- functional_gene_data %>%
          filter(ea_z_squared >= qchisq(0.95, 1))
        
      }
      functional_genes <- functional_gene_data %>%
        # Just for intergenic SNPs
        filter(is_positional == 0) %>%
        pull(ensembl_id) %>%
        unique()
      
    } else {
      functional_genes <- NULL
    }
    
    # Next save the text file of just the genes to pass as input:
    signal_gene_file <-
      file(
        paste0("data/adapt_results/gene_lists/signals/",
               assign_type, "/", phenotype, "_rsquared", rsquared,
               ifelse(!include_functional, "_pos_only_", "_"),
               "disc_genes.txt")
      )
    unique(c(signal_genes, non_signal_genes, functional_genes)) %>%
      paste(collapse = "\n") %>%
      writeLines(signal_gene_file)
    close(signal_gene_file)
    
    return(paste0("Finished saving signal gene list for ",
                  assign_type, " assignment using r^2 = 0.", rsquared))
  }


# Save the results for each combination -----------------------------------

# Positional + eSNPs
# ASD
save_signal_gene_lists("positional_esnps", "asd", 25, include_functional = TRUE)
save_signal_gene_lists("positional_esnps", "asd", 25, include_functional = FALSE)
# SCZ
save_signal_gene_lists("positional_esnps", "scz", 25, include_functional = TRUE)
save_signal_gene_lists("positional_esnps", "scz", 25, include_functional = FALSE)
# EA
save_signal_gene_lists("positional_esnps", "ea", 25, include_functional = TRUE)
save_signal_gene_lists("positional_esnps", "ea", 25, include_functional = FALSE)

# Positional
save_signal_gene_lists("positional", "asd", 25, include_functional = FALSE)
save_signal_gene_lists("positional", "scz", 25, include_functional = FALSE)
save_signal_gene_lists("positional", "ea", 25, include_functional = FALSE)

