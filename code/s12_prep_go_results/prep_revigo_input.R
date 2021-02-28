# PURPOSE: Prep GO input for REVIGO

library(tidyverse)
library(msigdbr)

# Write helper functions to generate REVIGO input -------------------------

msigdbr_data <- msigdbr(category = "C5") %>%
  dplyr::select(gs_subcat, gs_exact_source, gs_name) %>%
  distinct()

load_fuma_go_results <- function(list_type, assign_type, phenotype,
                                 gene_type = "all",
                                 is_matched = FALSE, n_matches = NULL) {
  
  read_tsv(paste0("data/adapt_results/fuma/output/",
                  list_type, "/", assign_type, "/", phenotype,
                  ifelse(gene_type == "all", "_",
                         paste0("_", gene_type, "_")),
                  "rsquared25",
                  ifelse(is_matched, paste0("_", n_matches, "match"),
                         ""),
                  "/GS.txt"))
  
}

# Write a helper function to get the GO ID and save the tsv file for GO input
save_revigo_input <- function(list_type, assign_type, phenotype, gene_type = "all",
                              go_type, is_matched = FALSE, n_matches = NULL,
                              go_data = msigdbr_data) {
  #browser()
  revigo_input <- load_fuma_go_results(list_type, assign_type, phenotype, gene_type,
                                       is_matched, n_matches) %>%
    filter(Category == paste0("GO_", go_type)) %>%
    dplyr::left_join(go_data, by = c("GeneSet" = "gs_name")) %>%
    filter(!is.na(gs_exact_source))
  
  if (nrow(revigo_input) == 0) {
    return(print("NO GO!"))
  } else {
    revigo_input %>%
      dplyr::select(gs_exact_source, p) %>%
      write_tsv(paste0("data/adapt_results/revigo_input/",
                       assign_type, "/", phenotype, "_", list_type,
                       ifelse(gene_type == "all", "_",
                              paste0("_", gene_type, "_")),
                       ifelse(is_matched, paste0(n_matches, "_match_"),
                              ""),
                       go_type, "_input.tsv"), col_names = FALSE)
    return(print("Saved GO results"))
  }
}



# Create REVIGO input -----------------------------------------------------

# ASD BP
save_revigo_input("signals", "positional_esnps", "asd", gene_type = "all", "bp")
# ASD CC
save_revigo_input("signals", "positional_esnps", "asd", gene_type = "all", "cc")

# SCZ BP
save_revigo_input("signals", "positional_esnps", "scz", gene_type = "all", "bp")
# SCZ CC
save_revigo_input("signals", "positional_esnps", "scz", gene_type = "all", "cc")
# SCZ MF
save_revigo_input("signals", "positional_esnps", "scz", gene_type = "all",
                  go_type = "mf")


# SCZ BP and CC with 2 matched nulls
save_revigo_input("signals", "positional_esnps", "scz", gene_type = "all", "bp",
                  is_matched = TRUE, n_matches = 2)
save_revigo_input("signals", "positional_esnps", "scz", gene_type = "all", "cc",
                  is_matched = TRUE, n_matches = 2)
save_revigo_input("signals", "positional_esnps", "scz", gene_type = "all", "mf",
                  is_matched = TRUE, n_matches = 2)




