# PURPOSE: Examine gene counts by type and compute their ratios, and do this for
#          the different rsquared thresholds

library(tidyverse)

read_adapt_gene_locus_table <- function(assign_type, phenotype, rsquared = 25) {
  read_csv(paste0("data/adapt_results/",
                  assign_type, "/", phenotype,
                  "/snp_gene_tables/rsquared", rsquared, "_gene_locus_table.csv"))
}


# Get the results for each phenotype --------------------------------------

asd_pos_esnps_gene_locus_25 <- read_adapt_gene_locus_table("positional_esnps", "asd")
scz_pos_esnps_gene_locus_25 <- read_adapt_gene_locus_table("positional_esnps", "scz")
ea_pos_esnps_gene_locus_25 <- read_adapt_gene_locus_table("positional_esnps", "ea")

asd_pos_gene_locus_25 <- read_adapt_gene_locus_table("positional", "asd")
scz_pos_gene_locus_25 <- read_adapt_gene_locus_table("positional", "scz")
ea_pos_gene_locus_25 <- read_adapt_gene_locus_table("positional", "ea")

# Ratio of clustered to unclustered ---------------------------------------

compute_ratio <- . %>%
  group_by(locus_type) %>%
  count() %>%
  ungroup() %>%
  pivot_wider(names_from = locus_type, values_from = n) %>%
  mutate(ratio = gene / locus, n_genes = gene + locus,
         perc_unclustered = gene / n_genes)

# First for Positional + eSNPs
asd_pos_esnps_gene_locus_25 %>% compute_ratio
# # A tibble: 1 x 4
#     gene locus ratio n_genes perc_unclustered
#     <int> <int> <dbl>   <int>            <dbl>
#   1    78   405 0.193     483            0.161
asd_pos_esnps_gene_locus_25 %>%
  group_by(locus_type) %>%
  summarize(n_loci = length(unique(ld_locus_id)))
# # A tibble: 2 x 2
#     locus_type n_loci
#     <chr>       <int>
#   1 gene           78
#   2 locus          47

scz_pos_esnps_gene_locus_25 %>% compute_ratio
#       gene locus ratio n_genes perc_unclustered
#     <int> <int> <dbl>   <int>            <dbl>
#   1  1566  3779 0.414    5345            0.293
scz_pos_esnps_gene_locus_25 %>%
  group_by(locus_type) %>%
  summarize(n_loci = length(unique(ld_locus_id)))
# # A tibble: 2 x 2
#     locus_type n_loci
#     <chr>       <int>
#   1 gene         1566
#   2 locus         711

ea_pos_esnps_gene_locus_25 %>% compute_ratio
# # A tibble: 1 x 5
#     gene locus ratio n_genes perc_unclustered
#     <int> <int> <dbl>   <int>            <dbl>
#   1  5091  7471 0.681   12562            0.405
ea_pos_esnps_gene_locus_25 %>%
  group_by(locus_type) %>%
  summarize(n_loci = length(unique(ld_locus_id)))
# # A tibble: 2 x 2
#     locus_type n_loci
#     <chr>       <int>
#   1 gene         5091
#   2 locus        1507

# Next for Positional
asd_pos_gene_locus_25 %>% compute_ratio
#      gene locus ratio n_genes perc_unclustered
#     <int> <int> <dbl>   <int>            <dbl>
#   1    77   370 0.208     447            0.172
asd_pos_gene_locus_25 %>%
  group_by(locus_type) %>%
  summarize(n_loci = length(unique(ld_locus_id)))
# # A tibble: 2 x 2
#     locus_type n_loci
#     <chr>       <int>
#   1 gene           77
#   2 locus          54

scz_pos_gene_locus_25 %>% compute_ratio
# # A tibble: 1 x 5
#     gene locus ratio n_genes perc_unclustered
#     <int> <int> <dbl>   <int>            <dbl>
#   1  1439  3319 0.434    4758            0.302
scz_pos_gene_locus_25 %>%
  group_by(locus_type) %>%
  summarize(n_loci = length(unique(ld_locus_id)))
# # A tibble: 2 x 2
#     locus_type n_loci
#     <chr>       <int>
#   1 gene         1439
#   2 locus         644

ea_pos_gene_locus_25 %>% compute_ratio
# # A tibble: 1 x 5
#     gene locus ratio n_genes perc_unclustered
#     <int> <int> <dbl>   <int>            <dbl>
#   1  4991  6775 0.737   11766            0.424
ea_pos_gene_locus_25 %>%
  group_by(locus_type) %>%
  summarize(n_loci = length(unique(ld_locus_id)))
# # A tibble: 2 x 2
#     locus_type n_loci
#     <chr>       <int>
#   1 gene         4991
#   2 locus        1431


# Sum of unique signals ---------------------------------------------------

length(unique(asd_pos_esnps_gene_locus_25$ld_locus_id))
# [1] 125

length(unique(scz_pos_esnps_gene_locus_25$ld_locus_id))
# [1] 2277

length(unique(ea_pos_esnps_gene_locus_25$ld_locus_id))
# [1] 6598



# Repeat for higher thresholds --------------------------------------------


# r^2 = 0.50 --------------------------------------------------------------

asd_pos_esnps_gene_locus_50 <- read_adapt_gene_locus_table("positional_esnps", "asd",
                                                           rsquared = 50)
scz_pos_esnps_gene_locus_50 <- read_adapt_gene_locus_table("positional_esnps", "scz",
                                                           rsquared = 50)
ea_pos_esnps_gene_locus_50 <- read_adapt_gene_locus_table("positional_esnps", "ea",
                                                          rsquared = 50)

asd_pos_gene_locus_50 <- read_adapt_gene_locus_table("positional", "asd",
                                                     rsquared = 50)
scz_pos_gene_locus_50 <- read_adapt_gene_locus_table("positional", "scz",
                                                     rsquared = 50)
ea_pos_gene_locus_50 <- read_adapt_gene_locus_table("positional", "ea",
                                                    rsquared = 50)

# First for Positional + eSNPs
asd_pos_esnps_gene_locus_50 %>% compute_ratio
# # A tibble: 1 x 5
#     gene locus ratio n_genes perc_unclustered
#     <int> <int> <dbl>   <int>            <dbl>
#   1   112   275 0.407     387            0.289
asd_pos_esnps_gene_locus_50 %>%
  group_by(locus_type) %>%
  summarize(n_loci = length(unique(ld_locus_id)))
# # A tibble: 2 x 2
#     locus_type n_loci
#     <chr>       <int>
#   1 gene          112
#   2 locus          38

scz_pos_esnps_gene_locus_50 %>% compute_ratio
# # A tibble: 1 x 5
#     gene locus ratio n_genes perc_unclustered
#     <int> <int> <dbl>   <int>            <dbl>
#   1  2306  2250  1.02    4556            0.506
scz_pos_esnps_gene_locus_50 %>%
  group_by(locus_type) %>%
  summarize(n_loci = length(unique(ld_locus_id)))
# # A tibble: 2 x 2
#     locus_type n_loci
#     <chr>       <int>
#   1 gene         2306
#   2 locus         499

ea_pos_esnps_gene_locus_50 %>% compute_ratio
# # A tibble: 1 x 5
#      gene locus ratio n_genes perc_unclustered
#     <int> <int> <dbl>   <int>            <dbl>
#   1  7054  4484  1.57   11538            0.611
ea_pos_esnps_gene_locus_50 %>%
  group_by(locus_type) %>%
  summarize(n_loci = length(unique(ld_locus_id)))
# # A tibble: 2 x 2
# locus_type n_loci
# <chr>       <int>
#   1 gene         7054
# 2 locus        1098

# Next for Positional
asd_pos_gene_locus_50 %>% compute_ratio
# # A tibble: 1 x 5
#     gene locus ratio n_genes perc_unclustered
#     <int> <int> <dbl>   <int>            <dbl>
#   1    98   212 0.462     310            0.316
asd_pos_gene_locus_50 %>%
  group_by(locus_type) %>%
  summarize(n_loci = length(unique(ld_locus_id)))
# # A tibble: 2 x 2
# locus_type n_loci
# <chr>       <int>
#   1 gene           98
# 2 locus          34

scz_pos_gene_locus_50 %>% compute_ratio
# # A tibble: 1 x 5
#     gene locus ratio n_genes perc_unclustered
#     <int> <int> <dbl>   <int>            <dbl>
#   1  2318  1929  1.20    4247            0.546
scz_pos_gene_locus_50 %>%
  group_by(locus_type) %>%
  summarize(n_loci = length(unique(ld_locus_id)))
# # A tibble: 2 x 2
# locus_type n_loci
# <chr>       <int>
#   1 gene         2318
# 2 locus         453

ea_pos_gene_locus_50 %>% compute_ratio
# # A tibble: 1 x 5
#     gene locus ratio n_genes perc_unclustered
#     <int> <int> <dbl>   <int>            <dbl>
#   1  6882  4009  1.72   10891            0.632
ea_pos_gene_locus_50 %>%
  group_by(locus_type) %>%
  summarize(n_loci = length(unique(ld_locus_id)))
# # A tibble: 2 x 2
#     locus_type n_loci
#     <chr>       <int>
#   1 gene         6882
#   2 locus        1025


# Repeat for r^2 = 0.75 ---------------------------------------------------


asd_pos_esnps_gene_locus_75 <- read_adapt_gene_locus_table("positional_esnps", "asd",
                                                           rsquared = 75)
scz_pos_esnps_gene_locus_75 <- read_adapt_gene_locus_table("positional_esnps", "scz",
                                                           rsquared = 75)
ea_pos_esnps_gene_locus_75 <- read_adapt_gene_locus_table("positional_esnps", "ea",
                                                          rsquared = 75)

asd_pos_gene_locus_75 <- read_adapt_gene_locus_table("positional", "asd",
                                                     rsquared = 75)
scz_pos_gene_locus_75 <- read_adapt_gene_locus_table("positional", "scz",
                                                     rsquared = 75)
ea_pos_gene_locus_75 <- read_adapt_gene_locus_table("positional", "ea",
                                                    rsquared = 75)

# First for Positional + eSNPs
asd_pos_esnps_gene_locus_75 %>% compute_ratio
# # A tibble: 1 x 5
#      gene locus ratio n_genes perc_unclustered
#      <int> <int> <dbl>   <int>            <dbl>
#   1   179   119  1.50     298            0.601
asd_pos_esnps_gene_locus_75 %>%
  group_by(locus_type) %>%
  summarize(n_loci = length(unique(ld_locus_id)))
# # A tibble: 2 x 2
# locus_type n_loci
# <chr>       <int>
#   1 gene          179
# 2 locus          21

scz_pos_esnps_gene_locus_75 %>% compute_ratio
# # A tibble: 1 x 5
# gene locus ratio n_genes perc_unclustered
# <int> <int> <dbl>   <int>            <dbl>
#   1  3297  1127  2.93    4424            0.745
scz_pos_esnps_gene_locus_75 %>%
  group_by(locus_type) %>%
  summarize(n_loci = length(unique(ld_locus_id)))
# # A tibble: 2 x 2
#     locus_type n_loci
#     <chr>       <int>
#   1 gene         3297
#   2 locus         323

ea_pos_esnps_gene_locus_75 %>% compute_ratio
# # A tibble: 1 x 5
#     gene locus ratio n_genes perc_unclustered
#     <int> <int> <dbl>   <int>            <dbl>
#   1  9071  2159  4.20   11230            0.808
ea_pos_esnps_gene_locus_75 %>%
  group_by(locus_type) %>%
  summarize(n_loci = length(unique(ld_locus_id)))
# # A tibble: 2 x 2
# locus_type n_loci
# <chr>       <int>
#   1 gene         9071
# 2 locus         652

# Next for Positional
asd_pos_gene_locus_75 %>% compute_ratio
# # A tibble: 1 x 5
# gene locus ratio n_genes perc_unclustered
# <int> <int> <dbl>   <int>            <dbl>
#   1   150    91  1.65     241            0.622
asd_pos_gene_locus_75 %>%
  group_by(locus_type) %>%
  summarize(n_loci = length(unique(ld_locus_id)))
# # A tibble: 2 x 2
# locus_type n_loci
# <chr>       <int>
#   1 gene          150
# 2 locus          22

scz_pos_gene_locus_75 %>% compute_ratio
# # A tibble: 1 x 5
# gene locus ratio n_genes perc_unclustered
# <int> <int> <dbl>   <int>            <dbl>
#   1  3098   952  3.25    4050            0.765
scz_pos_gene_locus_75 %>%
  group_by(locus_type) %>%
  summarize(n_loci = length(unique(ld_locus_id)))
# # A tibble: 2 x 2
# locus_type n_loci
# <chr>       <int>
#   1 gene         3098
# 2 locus         274

ea_pos_gene_locus_75 %>% compute_ratio
  # # A tibble: 1 x 5
  #      gene locus ratio n_genes perc_unclustered
  #     <int> <int> <dbl>   <int>            <dbl>
  #   1  8683  1918  4.53   10601            0.819
ea_pos_gene_locus_75 %>%
  group_by(locus_type) %>%
  summarize(n_loci = length(unique(ld_locus_id)))
# # A tibble: 2 x 2
#     locus_type n_loci
#     <chr>       <int>
#   1 gene         6882
#   2 locus        1025






