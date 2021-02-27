# PURPOSE: Initialize the table of possible genes using GENCODE v21, matching the
#          gene model used by Werling et al for BrainVar

library(tidyverse)
library(data.table)


# Load the GENCODE v21 data -----------------------------------------------

gencode_v21_gene_info <-
  rtracklayer::readGFF("data/gencode/gencode.v21.annotation.gtf.gz") %>%
  # Make the ensembl_id column:
  as.data.table() %>%
  .[, ensembl_id := str_sub(gene_id, 1, 15)] %>%
  # Filter to only the Gencode genes:
  .[type %in% "gene",] %>%
  # only keep autosomal chromosomes:
  mutate(gene_chr = as.numeric(str_remove(seqid, "chr"))) %>%
  filter(!is.na(gene_chr))


# Create a simplified table -----------------------------------------------

gene_type_table <- gencode_v21_gene_info %>%
  dplyr::select(ensembl_id, gene_chr, start, end, strand, gene_name, gene_type) %>%
  # Create the biotype columns:
  mutate(gene_biotype = case_when(
    gene_type == "protein_coding" ~ "protein_coding",
    gene_type == "antisense" ~ "antisense",
    gene_type %in% c("3prime_overlapping_ncrna", "bidirectional_promoter_lncRNA",
                     "lincRNA", "macro_lncRNA", "non_coding", "processed_transcript",
                     "sense_intronic", "sense_overlapping") ~ "lnc_rna",
    TRUE ~ "other")) %>%
  # Pivot biotype to be indicators:
  mutate(gene_biotype_member = gene_biotype,
         is_member = 1) %>%
  pivot_wider(names_from = gene_biotype_member, values_from = is_member,
              values_fill = list(is_member = 0),
              names_prefix = "biotype_")


# Save this table:
write_csv(gene_type_table,
          "data/gencode/gencode_v21_table.csv")

