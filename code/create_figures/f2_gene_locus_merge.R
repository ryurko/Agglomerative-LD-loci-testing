# PURPOSE: Create figure comparing the number of genes/loci as a function of
#          the LD r^2 threshold

library(tidyverse)
library(cowplot)

# Load the tidy gene to merged set tables  --------------------------------

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

            # Return a table with the number of merged gene sets along with the
            # original number (distinct will just remove duplicates later)
            table_types[i,] %>%
              mutate(n_loci = length(unique(result_data$ld_locus_id))) %>%
              bind_rows(
                tibble("rsquared" = 1, "assignment" = table_types$assignment[i],
                       "n_loci" = nrow(result_data))
              ) %>%
              return

          }) %>%
  distinct()



# Create the plot with the change in number of gene sets ------------------
library(latex2exp)

gene_merge_plot <- gene_merge_counts %>%
  ggplot(aes(x = rsquared, y = n_loci, color = assignment)) +
  geom_line(size = .5, alpha = 0.8) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("darkred", "darkblue"),
                     labels = c("Positional", "Positional + eSNPs"),
                     guide = guide_legend(reverse = TRUE)) +
  scale_x_continuous(limits = c(.25, 1.1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.50", "0.75",
                                "1\n(w/o merging)")) +
  scale_y_continuous(limits = c(25000, max(gene_merge_counts$n_loci + 100))) +
  theme_bw() +
  labs(x = TeX('LD $r^2$ threshold'),
       y = "Number of genes/loci",
       color = "Assignment") +
  theme(legend.title = element_blank(),
        legend.position = c(.65, .2),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

save_plot("figures/main/f_gene_merge.jpg",
          gene_merge_plot, ncol = 1, nrow = 1, base_asp = 1.2)
save_plot("figures/main/f_gene_merge.pdf",
          gene_merge_plot, ncol = 1, nrow = 1, base_asp = 1.2)
