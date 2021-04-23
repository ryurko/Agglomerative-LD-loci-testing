# PURPOSE: Initialize the ASD AdaPT LD loci level model datasets

library(tidyverse)
library(data.table)

# Load the LD loci p-value results ---------------------------------------

rsquared_thresholds <- c(.25, .50, .75)

# Positional + eSNPs
pos_esnps_ld_loci_pval_list <-
  map(rsquared_thresholds,
      function(rsquared_cut) {
        map_dfr(list.files(
          paste0(
            "data/ld_loci_level_pvals/rsquared", rsquared_cut * 100,
            "/positional_esnps"),
          full.names = TRUE),
          function(file_name) read_csv(file_name))
      })

# Positional
pos_ld_loci_pval_list <-
  map(rsquared_thresholds,
      function(rsquared_cut) {
        map_dfr(list.files(
          paste0(
            "data/ld_loci_level_pvals/rsquared", rsquared_cut * 100,
            "/positional"),
          full.names = TRUE),
          function(file_name) read_csv(file_name))
      })

# Load the gene to LD loci tables ----------------------------------------

# Positional + eSNPs
pos_esnp_gene_ld_loci_list <-
  map(rsquared_thresholds,
      function(rsquared_cut) {
        read_csv(
          paste0(
            "data/tidy_gene_ld_loci/positional_esnps/agglom_rsquared",
            rsquared_cut * 100, "_ld_loci.csv")
        )})

# Positional
pos_gene_ld_loci_list <-
  map(rsquared_thresholds,
      function(rsquared_cut) {
        read_csv(
          paste0(
            "data/tidy_gene_ld_loci/positional/agglom_rsquared",
            rsquared_cut * 100, "_ld_loci.csv")
        )})



# Load gene info ----------------------------------------------------------

gene_type_table <-
  read_csv("data/gencode/gencode_v21_table.csv")

# Load BrainVar eQTL data -------------------------------------------------

# First the SNP-gene pair level data:
brainvar_eqtl_data <-
  read_csv("data/eqtls/brainvar/esnp_eqtl_results.csv")

# since this data is in hg38 - load the hg19_id to join for later:
tidy_esnp_loc_table <- read_csv("data/tidy_snp_gene/esnps.csv")
brainvar_eqtl_data <- brainvar_eqtl_data %>%
  dplyr::left_join(dplyr::select(tidy_esnp_loc_table, hg38_id, hg19_id) %>%
                     distinct(),
                   by = c("esnp_hg38_id" = "hg38_id"))

# Next the WGCNA data:
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(sheet_i) readxl::read_excel(filename,
                                                           sheet = sheet_i))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  return(x)
}
brainvar_wgcna_sheets <- read_excel_allsheets("data/eqtls/brainvar/wgcna.xlsx")


# Create a pipeline for summarizing gene set BrainVar covariates at the LD loci
# level which is joined to this data based on the gene membership:
create_brainvar_eqtl_cov <- . %>%
  group_by(ld_loci_id) %>%
  summarize(n_brainvar_genes = length(unique(ensembl_id)),
            n_brainvar_esnps = length(unique(esnp_hg38_id)),
            ave_abs_pre_beta = mean(abs(Beta_Pre), na.rm = TRUE),
            max_abs_pre_beta = max(abs(Beta_Pre), na.rm = TRUE),
            var_abs_pre_beta = var(abs(Beta_Pre), na.rm = TRUE),
            ave_abs_post_beta = mean(abs(Beta_Post), na.rm = TRUE),
            max_abs_post_beta = max(abs(Beta_Post), na.rm = TRUE),
            var_abs_post_beta = var(abs(Beta_Post), na.rm = TRUE),
            ave_abs_comp_beta = mean(abs(Beta_Comp), na.rm = TRUE),
            max_abs_comp_beta = max(abs(Beta_Comp), na.rm = TRUE),
            var_abs_comp_beta = var(abs(Beta_Comp), na.rm = TRUE)) %>%
  ungroup() %>%
  # Replace missinld_loci for the variances with 0 (these are due to just 1 SNP):
  mutate(var_abs_pre_beta = ifelse(is.na(var_abs_pre_beta),
                                   0, var_abs_pre_beta),
         var_abs_post_beta = ifelse(is.na(var_abs_post_beta),
                                    0, var_abs_post_beta),
         var_abs_comp_beta = ifelse(is.na(var_abs_comp_beta),
                                    0, var_abs_comp_beta),
         brainvar_esnps_to_eqtl_genes = n_brainvar_esnps / n_brainvar_genes)


# Create pipeline to join WGCNA indicators (allowing for missing in case the
# gene set does not include any with sufficient expression coverage)
create_wgcna_indicators <- . %>%
  dplyr::select(ld_loci_id, ensembl_id) %>%
  dplyr::left_join(dplyr::select(brainvar_wgcna_sheets$`Module Genes`,
                                 gene_id, Module),
                   by = c("ensembl_id" = "gene_id")) %>%
  mutate(is_member = 1) %>%
  pivot_wider(names_from = Module, values_from = is_member,
              values_fill = list(is_member = 0),
              names_prefix = "brainvar_any_gene_") %>%
  dplyr::select(-ensembl_id) %>%
  group_by(ld_loci_id) %>%
  summarise_all(max) %>%
  ungroup()



# Load GTEx eQTL data -----------------------------------------------------

# First load the separate files for the two set of eQTLs:
gtex_fcba9_eqtl_data <-
  read_csv("data/eqtls/gtex/fcba9_esnp_eqtl_results.csv")
gtex_accba24_eqtl_data <-
  read_csv("data/eqtls/gtex/accba24_esnp_eqtl_results.csv")

# Next, need to initialize the dataset of values for the eQTLs in one tissue -
# that were not eQTLs in the other. The purpose of this is to be analogous
# to the BrainVar setting where there are eQTLs from one point of the study
# but no in another stage.

# First make new identifiers that are variant-gene IDs based on the style
# used in the GTEx data - with the GTEx alleles to simplify the matching:
gtex_fcba9_eqtl_data <- gtex_fcba9_eqtl_data %>%
  mutate(variant_id = paste(chr, bp, a1, a2, build, sep = "_"),
         variant_gene_id = paste0(variant_id, "_", gene_id))
gtex_accba24_eqtl_data <- gtex_accba24_eqtl_data %>%
  mutate(variant_id = paste(chr, bp, a1, a2, build, sep = "_"),
         variant_gene_id = paste0(variant_id, "_", gene_id))

# First load the Frontal Cortex BA9 dataset with the IDs that are just in the
# accba24 eQTL data - using the variant_gene_id so that this captures alleles:
gtex_fcba9_all_pairs_data <-
  fread("gunzip -c data/eqtls/gtex/Brain_Frontal_Cortex_BA9.allpairs.txt.gz") %>%
  .[, variant_gene_id := paste0(variant_id, "_", gene_id)] %>%
  .[variant_gene_id %in% gtex_accba24_eqtl_data$variant_gene_id,]

# Next for the reverse:
gtex_accba24_all_pairs_data <-
  fread("gunzip -c data/eqtls/gtex/Brain_Anterior_cingulate_cortex_BA24.allpairs.txt.gz") %>%
  .[, variant_gene_id := paste0(variant_id, "_", gene_id)] %>%
  .[variant_gene_id %in% gtex_fcba9_eqtl_data$variant_gene_id,]

# Next, make simplified datasets for both of these, just retaining the various
# identifiers for joining, and the SNP-gene pair slope for summarizing:
gtex_fcba9_eqtl_summary <- gtex_fcba9_eqtl_data %>%
  dplyr::select(variant_gene_id, ensembl_id, hg19_id, slope) %>%
  dplyr::rename(fcba9_slope = slope) %>%
  dplyr::left_join({
    gtex_accba24_all_pairs_data %>%
      as_tibble() %>%
      dplyr::select(variant_gene_id, slope) %>%
      dplyr::rename(accba24_slope = slope)
  }, by = "variant_gene_id") %>%
  # Replace missinld_loci with 0 since that represents no expression:
  mutate(accba24_slope = ifelse(is.na(accba24_slope), 0, accba24_slope))

# Next for accba24:
gtex_accba24_eqtl_summary <- gtex_accba24_eqtl_data %>%
  dplyr::select(variant_gene_id, ensembl_id, hg19_id, slope) %>%
  dplyr::rename(accba24_slope = slope) %>%
  dplyr::left_join({
    gtex_fcba9_all_pairs_data %>%
      as_tibble() %>%
      dplyr::select(variant_gene_id, slope) %>%
      dplyr::rename(fcba9_slope = slope)
  }, by = "variant_gene_id") %>%
  # Replace missinld_loci with 0 since that represents no expression:
  mutate(fcba9_slope = ifelse(is.na(fcba9_slope), 0, fcba9_slope))

# Now make a table denoting the type of eQTL for each SNP-gene pair:
gtex_eqtl_variant_types <-
  bind_rows(mutate(gtex_fcba9_eqtl_summary,
                   is_fcba9_eqtl = 1, is_accba24_eqtl = 0),
            mutate(gtex_accba24_eqtl_summary,
                   is_accba24_eqtl = 1, is_fcba9_eqtl = 0)) %>%
  group_by(variant_gene_id) %>%
  summarize(is_fcba9_eqtl = max(is_fcba9_eqtl),
            is_accba24_eqtl = max(is_accba24_eqtl))

# Finally stack the GTEx eQTL pairs together and join this eQTL type information:
gtex_eqtl_summary <- gtex_fcba9_eqtl_summary %>%
  bind_rows(gtex_accba24_eqtl_summary) %>%
  distinct() %>%
  dplyr::left_join(gtex_eqtl_variant_types, by = "variant_gene_id")

# Create a pipeline for summarizing gene set GTEx covariates:
create_gtex_eqtl_cov <- . %>%
  group_by(ld_loci_id) %>%
  summarize(n_gtex_genes = length(unique(ensembl_id)),
            n_gtex_esnps = length(unique(hg19_id)),
            is_fcba9_eqtl = max(is_fcba9_eqtl, na.rm = TRUE),
            is_accba24_eqtl = max(is_accba24_eqtl, na.rm = TRUE),
            ave_abs_fcba9_slope = mean(abs(fcba9_slope), na.rm = TRUE),
            max_abs_fcba9_slope = max(abs(fcba9_slope), na.rm = TRUE),
            var_abs_fcba9_slope = var(abs(fcba9_slope), na.rm = TRUE),
            ave_abs_accba24_slope = mean(abs(accba24_slope), na.rm = TRUE),
            max_abs_accba24_slope = max(abs(accba24_slope), na.rm = TRUE),
            var_abs_accba24_slope = var(abs(accba24_slope), na.rm = TRUE)) %>%
  mutate_all(function(x) ifelse(is.infinite(x), NA, x)) %>%
  ungroup() %>%
  # Replace missinld_loci for the variances with 0 (these are due to just 1 SNP):
  mutate(var_abs_fcba9_slope = ifelse(is.na(var_abs_fcba9_slope),
                                      0, var_abs_fcba9_slope),
         var_abs_accba24_slope = ifelse(is.na(var_abs_accba24_slope),
                                        0, var_abs_accba24_slope),
         gtex_esnps_to_eqtl_genes = n_gtex_esnps / n_gtex_genes)


# Next join the GTEx WGCNA indicators:
gtex_cortical_wgcna_data <-
  read_csv("data/eqtls/gtex/wgcna/output/cortical_wgcna_results.csv")

create_gtex_wgcna_indicators <- . %>%
  dplyr::select(ld_loci_id, ensembl_id) %>%
  dplyr::left_join(gtex_cortical_wgcna_data, by = "ensembl_id") %>%
  mutate(is_member = 1) %>%
  pivot_wider(names_from = wgcna_label, values_from = is_member,
              values_fill = list(is_member = 0),
              names_prefix = "gtex_any_gene_") %>%
  dplyr::select(-ensembl_id) %>%
  group_by(ld_loci_id) %>%
  summarise_all(max) %>%
  ungroup()


# Load the gnomad plof data -----------------------------------------------

gnomad_plof_data <-
  fread("gunzip -c data/gnomad/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz") %>%
  as_tibble()

# Pipeline to make loeuf summary:
make_loeuf_summary <- . %>%
  dplyr::select(ld_loci_id, ensembl_id) %>%
  left_join(dplyr::select(gnomad_plof_data, gene_id, oe_lof_upper),
            by = c("ensembl_id" = "gene_id")) %>%
  group_by(ld_loci_id) %>%
  summarize(max_loeuf = max(oe_lof_upper, na.rm = TRUE),
            min_loeuf = min(oe_lof_upper, na.rm = TRUE),
            mean_loeuf = mean(oe_lof_upper, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(max_loeuf = ifelse(is.infinite(max_loeuf),
                            NA, max_loeuf),
         min_loeuf = ifelse(is.infinite(min_loeuf),
                            NA, min_loeuf))

# Additional helper functions for pipeline --------------------------------

# Also replace the missinld_loci for eQTL slopes with 0
replace_missing_vals <- function(x) {
  ifelse(is.na(x), 0, x)
}

# Pipeline for addressing missing and denoting the eQTL type:
update_eqtl_vars <- . %>%
  # First make the BrainVar eQTL indicator:
  mutate(any_brainvar_eqtl = as.numeric(!is.na(ave_abs_pre_beta))) %>%
  # Fill in the missing eQTL slope variables with 0:
  mutate_at(vars(ends_with("beta")), replace_missing_vals) %>%
  mutate_at(vars(contains("fcba9")), replace_missing_vals) %>%
  mutate_at(vars(contains("accba24")), replace_missing_vals) %>%
  mutate(any_gtex_eqtl = as.numeric(pmax(is_fcba9_eqtl, is_accba24_eqtl) > 0))

# Add z-stat variables
add_z_stat_vars <- . %>%
  mutate(scz_quad_z = qnorm(1 - scz_quad_pval),
         asd_quad_z = qnorm(1 - asd_quad_pval),
         ea_quad_z = qnorm(1 - ea_quad_pval),
         cap_scz_quad_pval = pmax(1e-15, pmin(1 - 1e-15, scz_quad_pval)),
         cap_asd_quad_pval = pmax(1e-15, pmin(1 - 1e-15, asd_quad_pval)),
         cap_ea_quad_pval = pmax(1e-15, pmin(1 - 1e-15, ea_quad_pval)),
         cap_scz_quad_z = qnorm(1 - cap_scz_quad_pval),
         cap_asd_quad_z = qnorm(1 - cap_asd_quad_pval),
         cap_ea_quad_z = qnorm(1 - cap_ea_quad_pval))

# Add fake SCZ pvalues
make_fake_scz_pval <- . %>%
  dplyr::select(asd_quad_pval) %>%
  arrange(asd_quad_pval) %>%
  mutate(pval_rank = 1:n()) %>%
  dplyr::rename(fake_scz_quad_pval = asd_quad_pval)


# Function for creating model dataset -------------------------------------

create_adapt_model_data <- function(ld_loci_pval_summary, gene_ld_loci_table) {
  
  # First join the gene-type to the gene_ld_loci_table:
  ld_loci_gene_type_summary <- gene_ld_loci_table %>%
    dplyr::left_join(
      dplyr::select(gene_type_table,
                    ensembl_id, biotype_other, biotype_lnc_rna,
                    biotype_protein_coding, biotype_antisense),
      by = "ensembl_id") %>%
    dplyr::select(-ensembl_id) %>%
    group_by(ld_loci_id) %>%
    summarize_all(max) %>%
    ungroup()
  
  # Now join to the dataset with the p-values:
  ld_loci_data <- ld_loci_pval_summary %>%
    dplyr::left_join(ld_loci_gene_type_summary, by = "ld_loci_id")
  
  # join brainvar data:
  ld_loci_model_data <- ld_loci_data %>%
    dplyr::left_join({
      brainvar_eqtl_data %>%
        filter(ensembl_id %in% gene_ld_loci_table$ensembl_id) %>%
        dplyr::left_join(gene_ld_loci_table,
                         by = "ensembl_id") %>%
        create_brainvar_eqtl_cov
    }, by = "ld_loci_id") %>%
    dplyr::left_join(create_wgcna_indicators(gene_ld_loci_table),
                     by = "ld_loci_id") %>%
    # join gtex:
    dplyr::left_join({
      gtex_eqtl_summary %>%
        filter(ensembl_id %in% gene_ld_loci_table$ensembl_id) %>%
        dplyr::left_join(gene_ld_loci_table,
                         by = "ensembl_id") %>%
        create_gtex_eqtl_cov
    }, by = "ld_loci_id")  %>%
    dplyr::left_join(create_gtex_wgcna_indicators(gene_ld_loci_table),
                     by = "ld_loci_id") %>%
    dplyr::left_join(make_loeuf_summary(gene_ld_loci_table),
                     by = "ld_loci_id") %>%
    update_eqtl_vars() %>%
    add_z_stat_vars()
  
  # Add fake p-values
  ld_loci_model_data <- ld_loci_model_data %>%
    arrange(scz_quad_pval) %>%
    mutate(scz_pval_rank = 1:n()) %>%
    left_join(ld_loci_model_data %>% make_fake_scz_pval,
              by = c("scz_pval_rank" = "pval_rank"))
  
  # Make a table of the BrainVar and GTEx eQTL results together:
  eqtl_gene_esnp_counts <- dplyr::select(brainvar_eqtl_data,
                                         ensembl_id, hg19_id) %>%
    bind_rows(gtex_eqtl_summary %>%
                dplyr::select(ensembl_id, hg19_id)) %>%
    distinct() %>%
    dplyr::left_join(gene_ld_loci_table,
                     by = "ensembl_id") %>%
    group_by(ld_loci_id) %>%
    summarize(n_esnps = length(unique(hg19_id)),
              n_eqtl_genes = length(unique(ensembl_id)))
  
  # Join to LD loci then to the model data:
  ld_loci_model_data <- ld_loci_model_data %>%
    dplyr::left_join(eqtl_gene_esnp_counts, by = "ld_loci_id") %>%
    mutate(n_esnps = ifelse(is.na(n_esnps), 0, n_esnps),
           n_brainvar_esnps = ifelse(is.na(n_brainvar_esnps), 0, n_brainvar_esnps),
           n_gtex_esnps = ifelse(is.na(n_gtex_esnps), 0, n_gtex_esnps))
  
  # More gene counts:
  ld_loci_gene_count <- gene_ld_loci_table %>%
    group_by(ld_loci_id) %>%
    summarize(n_genes = n())
  
  ld_loci_model_data <- ld_loci_model_data %>%
    dplyr::left_join(ld_loci_gene_count, by = "ld_loci_id") %>%
    # Normalize SNPs by genes:
    mutate(snps_to_genes = n_snps / n_genes,
           # Fix missing eQTL genes:
           n_brainvar_genes = ifelse(is.na(n_brainvar_genes), 0, n_brainvar_genes),
           n_gtex_genes = ifelse(is.na(n_gtex_genes), 0, n_gtex_genes),
           n_eqtl_genes = ifelse(is.na(n_eqtl_genes), 0, n_eqtl_genes),
           # Ratio of eSNPs and eQTL genes:
           prop_esnps = n_esnps / n_snps,
           prop_eqtl_genes = n_eqtl_genes / n_genes)
  
  return(ld_loci_model_data)
  
}


# Create the AdaPT model datasets -----------------------------------------

# Positional + eSNPs
pos_esnps_model_data_list <-
  map(1:length(pos_esnps_ld_loci_pval_list),
      function(rsquared_i) {
        create_adapt_model_data(pos_esnps_ld_loci_pval_list[[rsquared_i]],
                                pos_esnp_gene_ld_loci_list[[rsquared_i]])
      })

# Positional
pos_model_data_list <-
  map(1:length(pos_ld_loci_pval_list),
      function(rsquared_i) {
        create_adapt_model_data(pos_ld_loci_pval_list[[rsquared_i]],
                                pos_gene_ld_loci_list[[rsquared_i]])
      })



# Save the datasets -------------------------------------------------------

walk(1:length(pos_esnps_model_data_list),
     function(rsquared_i) {
       write_csv(pos_esnps_model_data_list[[rsquared_i]],
                 paste0(
                   "data/adapt_model_data/positional_esnps/rsquared",
                   rsquared_thresholds[rsquared_i] * 100, "_model_data.csv"
                 ))
     })

walk(1:length(pos_model_data_list),
     function(rsquared_i) {
       write_csv(pos_model_data_list[[rsquared_i]],
                 paste0(
                   "data/adapt_model_data/positional/rsquared",
                   rsquared_thresholds[rsquared_i] * 100, "_model_data.csv"
                 ))
     })



