# PURPOSE: Initialize the GTEx WGCNA results using only the GENCODE v21 genes

library(tidyverse)
library(data.table)
library(WGCNA)

# Load the gene table -----------------------------------------------------

gene_type_table <-
  read_csv("data/gencode/gencode_v21_table.csv")

# Load reference files ----------------------------------------------------

# See data/eqtls/gtex/wgcna/input/GTEx_v7_Annotations_SampleAttributesDS.txt for all files and see
# data/eqtls/gtex/wgcna/input/GTEx_Analysis_v7_Annotations_SampleAttributesDD.xlsx for an
# appropriate glossary
ref_file <- fread("data/eqtls/gtex/wgcna/input/GTEx_v7_Annotations_SampleAttributesDS.txt",
                  verbose = FALSE)
ref_file_glossary <- readxl::read_excel(
  "data/eqtls/gtex/wgcna/input/GTEx_Analysis_v7_Annotations_SampleAttributesDD.xlsx")


# Load cortical samples ---------------------------------------------------

# This file initializes the expression datasets with just the specific
# tissue types for the WGCNA analysis using the reference file containing
# the tissue type information for each of the samples measured from

# Cortical tissue samples: either Brain - Anterior cingulate cortex (BA24) or Brain - Frontal Cortex (BA9)
cortical_samples <- fread("data/eqtls/gtex/wgcna/input/GTEx_v7_Annotations_SampleAttributesDS.txt",
                          verbose = FALSE) %>%
  filter(str_detect(SMTSD, "(BA24)|(BA9)")) %>%
  pull(SAMPID)


# Init TPM files ----------------------------------------------------------

# Next need to load the TPM file, only selecting columns matching these
# samples and saving the files:
tpm_counts <- fread(input = paste0("zcat < ",
                                   "data/eqtls/gtex/wgcna/input/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz"),
                    verbose = FALSE)
tpm_colnames <- colnames(tpm_counts)
match_cortical_samples <- cortical_samples[which(cortical_samples %in% tpm_colnames)]
tpm_counts[, c("Name", "Description", match_cortical_samples), with = FALSE] %>%
  write_csv("data/eqtls/gtex/wgcna/side_output/tpm_cortical_samples.csv")
rm(tpm_counts)

# Enable parallel processing for WGCNA with ten threads
# NOTE enableWGCNAThreads does not work - but this is a contradiction with the
# documentation unfortunately...
allowWGCNAThreads(nThreads = 10)

# Load cortical samples:
cortical_tpm_data <- readr::read_csv("data/eqtls/gtex/wgcna/side_output/tpm_cortical_samples.csv")

# Since the gene IDs provided in the raw GTEx data are in the GENCODE form,
# need to remove the .version part of them
cortical_tpm_data <- cortical_tpm_data %>%
  mutate(ensembl_id = str_sub(Name, 1, 15)) %>%
  # And then filter to only be genes in the GENCODE v21 data:
  filter(ensembl_id %in% gene_type_table$ensembl_id)

# Let's remove all genes with 0 expression for over 50% of the provided samples
zero_half_cortical_expression <- which(apply(dplyr::select(cortical_tpm_data,
                                                           -c(Name, ensembl_id,
                                                              Description)), 1,
                                             function(x) length(which(x == 0)) / length(x) >= .5))

# Now remove these genes
expressed_cortical_tpm_data <- cortical_tpm_data[-zero_half_cortical_expression,]
dim(expressed_cortical_tpm_data)
# [1] 26513   253
# similar to BrainVar based count from the coverage

# Grab only the samples columns and transpose so each row is a sample, but
# take the log base 2 transformation of the expressions (TPMs) + 1:
cortical_log_tpm <- apply(dplyr::select(expressed_cortical_tpm_data,
                                        -c(Name, ensembl_id,
                                           Description)),
                          1, function(x) log(x + 1, base = 2))

# Set the column names to be the gene ids:
colnames(cortical_log_tpm) <- expressed_cortical_tpm_data$ensembl_id

# Now need to choose which power to use based on
# the scale-free topology criterion
powers <- c(2:20)

# We're going to use the unsigned
cortical_power_search <- pickSoftThreshold(cortical_log_tpm,
                                           dataIsExpr = TRUE,
                                           powerVector = powers,
                                           corFnc = cor,
                                           corOptions = list(use = 'p'),
                                           networkType = "signed")

# Create two plots showing the R-squared and connectivity:
cortical_power_search$fitIndices %>%
  ggplot(aes(x = Power, y = SFT.R.sq)) +
  geom_hline(yintercept = 0.8, color = "darkorange") +
  geom_text(aes(label = Power), color = "darkblue",
            size = 5) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Soft threshold power",
       y = "Scale free topology model fit R^2",
       title = "Scale independence chart for cortical tissue samples") +
  theme_bw()

cortical_power_search$fitIndices %>%
  as_tibble() %>%
  ggplot(aes(x = Power, y = `mean.k.`)) +
  #geom_hline(yintercept = 0.8, color = "darkorange") +
  geom_text(aes(label = Power), color = "darkblue",
            size = 5) +
  #scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Soft threshold power",
       y = "Mean connectivity",
       title = "Scale independence chart for cortical tissue samples") +
  theme_bw()

# Will use a power of 8

# Now generate the WGCNA results using the blockwise implementation to run
# locally with the selected powers from above:

# Set the proper directory to the wcgna_results folder
#setwd("data/eqtls/gtex/wgcna/side_output")
cortical_bwnet <- blockwiseModules(cortical_log_tpm,
                                   maxBlockSize = 30000,
                                   power = 8, TOMType = "signed",
                                   minModuleSize = 200,
                                   detectCutHeight = 0.999,
                                   saveTOMs = TRUE,
                                   mergeCutHeight = 0.15,
                                   saveTOMFileBase = "cortical_TOM_blockwise",
                                   verbose = 3)
# Return back to the project directory

# Convert labels to colors for plotting
cortical_mergedColors <- labels2colors(cortical_bwnet$colors)

# Create a dataframes with the gene names and the labels:
cortical_wgcna_data <- tibble("ensembl_id" =  colnames(cortical_log_tpm),
                              "wgcna_label" = cortical_mergedColors)

# Save them to the wcgna_results folder
write_csv(cortical_wgcna_data,
          "data/eqtls/gtex/wgcna/output/cortical_wgcna_results.csv")


