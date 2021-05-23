# PURPOSE: Examine the number of genes and merged genes from the algorithm
library(tidyverse)

# Load the new updated counts ---------------------------------------------

# Want to load the results for each type of positional assignment and rsquared
# threshold, can do so via a wrapper function to make it easy to load with a
# table to iterate over:
table_types <- expand.grid("rsquared" = c(0.25, 0.5, 0.75),
                           "assignment" = c("positional", "positional_esnps")) %>%
  as_tibble()

# Iterate over these row values to get the number of merged sets, along with
# including a row for all genes (rsquared 1 indicating no merging):
gene_merge_counts <-
  map_dfr(1:nrow(table_types),
          function(i) {
            # Load the type of results:
            result_data <-
              read_csv(paste0("data/tidy_gene_locus/",
                              table_types$assignment[i],
                              "/agglom_rsquared",
                              table_types$rsquared[i] * 100,
                              "_ld_loci.csv"))

            # Get the number of single genes versus loci:
            locus_type_summary <- result_data %>%
              group_by(locus_type) %>%
              count() %>%
              pivot_wider(names_from = locus_type, values_from = n)

            locus_size_summary <- result_data %>%
              filter(locus_type == "locus") %>%
              group_by(ld_locus_id) %>%
              count()
            print(summary(locus_size_summary$n))


            # Return a table with the number of merged gene sets along with the
            # original number (distinct will just remove duplicates later)
            table_types[i,] %>%
              mutate(n_loci = length(unique(result_data$ld_locus_id)),
                     n_single_genes = locus_type_summary$gene,
                     n_multi_gene_loci = n_loci - n_single_genes,
                     n_genes_in_loci = locus_type_summary$locus) %>%
              bind_rows(
                tibble("rsquared" = 1, "assignment" = table_types$assignment[i],
                       "n_loci" = nrow(result_data))
              ) %>%
              return

          }) %>%
  distinct()

gene_merge_counts
# # A tibble: 8 x 6
#         rsquared assignment       n_loci n_single_genes n_multi_gene_loci n_genes_in_loci
#         <dbl>  <fct>             <int>          <int>             <int>           <int>
#   1     0.25   positional        27941          23956              3985           16625
#   2     1      positional        40581             NA                NA              NA
#   3     0.5    positional        33277          30410              2867           10171
#   4     0.75   positional        37114          35454              1660            5127
#   5     0.25   positional_esnps  27522          23386              4136           17915
#   6     1      positional_esnps  41301             NA                NA              NA
#   7     0.5    positional_esnps  33147          30103              3044           11198
#   8     0.75   positional_esnps  37542          35751              1791            5550

# Load the Positional + eSNPs results and see which is the largest locus:
pos_esnps_results <-
  read_csv(paste0("data/tidy_gene_locus/",
                  table_types$assignment[4],
                  "/agglom_rsquared",
                  table_types$rsquared[4] * 100,
                  "_ld_loci.csv"))

locus_size_summary <- pos_esnps_results %>%
  filter(locus_type == "locus") %>%
  group_by(ld_locus_id) %>%
  count() %>%
  arrange(desc(n))
locus_size_summary
# # A tibble: 4,136 x 2
# # Groups:   ld_locus_id [4,136]
# 1 chr11_7        67
# 2 chr16_9        43
# 3 chr16_18       38
# 4 chr5_5         38
# 5 chr8_4         38
# 6 chr3_11        37
# 7 chr6_20        36
# 8 chr10_1        33
# 9 chr19_49       32
# 10 chr15_36       31
# # … with 4,126 more rows

summary(locus_size_summary$n)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 2.000   2.000   3.000   4.331   5.000  67.000

