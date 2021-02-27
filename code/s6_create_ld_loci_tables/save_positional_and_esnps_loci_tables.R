# PURPOSE: Initialize the tidy SNP to LD loci datasets following the LD-based
#          LD loci agglomerative algorithm

library(tidyverse)

# Write helper functions for initializing data ----------------------------

# Function to make the tidy LD loci data:
create_tidy_ld_loci_table <- function(ld_loci_results) {
  map_dfr(1:length(ld_loci_results),
          function (chr_i) {
            map_dfr(1:length(ld_loci_results[[chr_i]]$ld_loci_list),
                    function (set_i) {
                      tibble(ensembl_id =
                               ld_loci_results[[chr_i]]$ld_loci_list[[set_i]],
                             ld_loci_id = paste0(names(ld_loci_results)[chr_i], "_", set_i))
                    })
          })
}

# Function to create list of results for given threshold and type:
create_rsquared_ld_loci_data <- function(rsquared, loci_type) {
  results <-
    map(list.files(
      paste0("data/ld_induced_loci/rsquared",
             rsquared * 100, "/", loci_type),
      full.names = TRUE),
      readRDS)
  names(results) <-
    str_remove(list.files(
      paste0("data/ld_induced_loci/rsquared",
             rsquared * 100, "/", loci_type)),
      "_results\\.rds")
  results %>%
    create_tidy_ld_loci_table() %>%
    return
}



# Initialize datasets -----------------------------------------------------

# Add .5 after finished on server
rsquared_vals <- c(0.25, 0.5)

# Create list for Positional + eSNPs:
pos_esnps_ld_loci_list <-
  map(rsquared_vals,
      function(rsquared_cut) create_rsquared_ld_loci_data(rsquared_cut, "positional_esnps"))
names(pos_esnps_ld_loci_list) <- paste0("rsquared = ", rsquared_vals)

# Positional:
pos_ld_loci_list <-
  map(rsquared_vals,
      function(rsquared_cut) create_rsquared_ld_loci_data(rsquared_cut, "positional"))
names(pos_ld_loci_list) <- paste0("rsquared = ", rsquared_vals)

# Save these gene to LD loci datasets:
walk(1:length(pos_esnps_ld_loci_list),
     function(loci_i) {
       write_csv(pos_esnps_ld_loci_list[[loci_i]],
                 paste0(
                   "data/tidy_gene_ld_loci/positional_esnps/agglom_rsquared",
                   rsquared_vals[loci_i] * 100, "_ld_loci.csv"
                 ))
     })

walk(1:length(pos_ld_loci_list),
     function(loci_i) {
       write_csv(pos_ld_loci_list[[loci_i]],
                 paste0(
                   "data/tidy_gene_ld_loci/positional/agglom_rsquared",
                   rsquared_vals[loci_i] * 100, "_ld_loci.csv"
                 ))
     })


# Create and save tidy SNP to LD loci data -------------------------------

# First load the tidy SNP-gene datasets
tidy_positional_esnp_gene_data <-
  read_csv("data/tidy_snp_gene/positional_esnps.csv")
tidy_positional_gene_data <-
  read_csv("data/tidy_snp_gene/positional.csv")

# Helper function to join LD loci to SNP-gene data:
create_tidy_snp_ld_loci_data <- function(tidy_snp_gene_data, tidy_gene_ld_loci_data) {
  tidy_snp_gene_data %>%
    dplyr::left_join(tidy_gene_ld_loci_data, by = "ensembl_id") %>%
    separate(ld_loci_id, into = c("gene_chr", "set_id"), sep = "_", remove = FALSE) %>%
    mutate(gene_chr = as.numeric(str_remove(gene_chr, "chr"))) %>%
    dplyr::select(-set_id, -ensembl_id) %>%
    # Take the distinct to remove the duplicate SNP-gene rows:
    distinct() %>%
    # Add the squared z stats columns:
    mutate(scz_z_squared = (log(scz_or) / scz_se)^2,
           asd_z_squared = (log(asd_or) / asd_se)^2,
           ea_z_squared = (ea_beta / ea_se)^2)
}

# Create list for all_combo:
pos_esnps_tidy_snp_ld_loci_list <-
  map(1:length(pos_esnps_ld_loci_list),
      function(gs_list_i) create_tidy_snp_ld_loci_data(tidy_positional_esnp_gene_data,
                                                  pos_esnps_ld_loci_list[[gs_list_i]]))
names(pos_esnps_tidy_snp_ld_loci_list) <- names(pos_esnps_ld_loci_list)

# Positional:
pos_tidy_snp_ld_loci_list <-
  map(1:length(pos_ld_loci_list),
      function(gs_list_i) create_tidy_snp_ld_loci_data(tidy_positional_gene_data,
                                                  pos_ld_loci_list[[gs_list_i]]))
names(pos_tidy_snp_ld_loci_list) <- names(pos_ld_loci_list)


# Save each dataset separately --------------------------------------------

walk(1:length(pos_esnps_tidy_snp_ld_loci_list),
     function(loci_i) {
       write_csv(pos_esnps_tidy_snp_ld_loci_list[[loci_i]],
                 paste0(
                   "data/tidy_snp_ld_loci/rsquared",
                   rsquared_vals[loci_i] * 100, "/positional_esnps.csv"
                 ))
     })

walk(1:length(pos_tidy_snp_ld_loci_list),
     function(loci_i) {
       write_csv(pos_tidy_snp_ld_loci_list[[loci_i]],
                 paste0(
                   "data/tidy_snp_ld_loci/rsquared",
                   rsquared_vals[loci_i] * 100, "/positional.csv"
                 ))
     })




