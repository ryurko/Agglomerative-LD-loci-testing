# PURPOSE: Prep GO input for REVIGO

library(tidyverse)
library(msigdbr)

# Write helper functions to generate REVIGO input -------------------------

msigdbr_data <- msigdbr(category = "C5") %>%
  dplyr::select(gs_subcat, gs_exact_source, gs_name) %>%
  distinct()

load_fuma_go_results <- function(list_type, assign_type, phenotype) {

  read_tsv(paste0("data/fuma/output/",
                  assign_type, "/", phenotype, "_n_snps_matched/GS.txt"))

}

# Write a helper function to get the GO ID and save the tsv file for GO input
save_revigo_input <- function(assign_type, phenotype,
                              go_type, go_data = msigdbr_data) {
  #browser()
  revigo_input <- load_fuma_go_results(list_type, assign_type, phenotype) %>%
    filter(Category == paste0("GO_", go_type)) %>%
    dplyr::left_join(go_data, by = c("GeneSet" = "gs_name")) %>%
    filter(!is.na(gs_exact_source))

  if (nrow(revigo_input) == 0) {
    return(print("NO GO!"))
  } else {
    revigo_input %>%
      dplyr::select(gs_exact_source, p) %>%
      write_tsv(paste0("data/revigo/input/",
                       assign_type, "/", phenotype, "/",
                       go_type, "_input.tsv"), col_names = FALSE)
    return(print("Saved GO results"))
  }
}



# Create REVIGO input -----------------------------------------------------



# SCZ with 2 matched nulls
save_revigo_input("positional_esnps", "scz",  "bp")
save_revigo_input("positional_esnps", "scz",  "cc")
save_revigo_input("positional_esnps", "scz",  "mf")

# EA with 1 matched null
save_revigo_input("positional_esnps", "ea",  "bp")
save_revigo_input("positional_esnps", "ea",  "cc")
save_revigo_input("positional_esnps", "ea",  "mf")




