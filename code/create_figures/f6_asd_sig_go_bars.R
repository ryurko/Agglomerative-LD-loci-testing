# PURPOSE: Create GO enrichment plots to display in manuscript and supplement
#          based on the FUMA output. This is only for the ASD results as the
#          SCZ results are displayed using the REVIGO treemap for ease.

library(tidyverse)
library(cowplot)

# Load the FUMA output of interest ----------------------------------------

# Write a helper function to load the LD loci GO enrichment results

load_fuma_go_results <- function(list_type, assign_type, phenotype,
                                 incl_fun = TRUE) {
  
  read_tsv(paste0("data/adapt_results/fuma/output/",
                  list_type, "/", assign_type, "/", phenotype,
                  ifelse(incl_fun, "_", "_pos_only_"),
                  "rsquared25",
                  "/GS.txt"))
  
}

# First the ASD Positional + eSNPs signals results:
asd_pos_esnps_signals_go_table <- load_fuma_go_results("signals", "positional_esnps", "asd") %>%
  filter(str_detect(Category, "GO_")) %>%
  mutate(list_type = "signals", assign_type = "positional_esnps", phenotype = "asd",
         incl_fun = TRUE)
# Next with all genes (no filtering)
asd_pos_esnps_all_go_table <- load_fuma_go_results("all", "positional_esnps", "asd") %>%
  filter(str_detect(Category, "GO_")) %>%
  mutate(list_type = "all", assign_type = "positional_esnps", phenotype = "asd",
         incl_fun = TRUE)

# Next the positional results for ASD with both signals and then
asd_pos_signals_go_table <- load_fuma_go_results("signals", "positional", "asd",
                                                 incl_fun = FALSE) %>%
  filter(str_detect(Category, "GO_")) %>%
  mutate(list_type = "signals", assign_type = "positional", phenotype = "asd",
         incl_fun = FALSE)
asd_pos_all_go_table <- load_fuma_go_results("all", "positional", "asd",
                                             incl_fun = FALSE) %>%
  filter(str_detect(Category, "GO_")) %>%
  mutate(list_type = "all", assign_type = "positional", phenotype = "asd",
         incl_fun = FALSE)

# Create GO barcharts for main results ------------------------------------

# Create barcharts of GO enrichment for BP and CC results with the main set of
# results corresponding to ASD Positional + eSNPs based on signals

asd_pos_esnps_sig_go_bars <- asd_pos_esnps_signals_go_table %>%
  filter(Category %in% c("GO_bp", "GO_cc")) %>%
  mutate(Category = fct_relevel(Category, "GO_bp", "GO_cc"),
         Category = fct_recode(Category, `Biological processes` = "GO_bp",
                               `Cellular components` = "GO_cc"),
         neg_log_pval = -log10(adjP)) %>%
  mutate(GeneSet = str_remove(GeneSet, "^GO_") %>% str_replace_all("_", " ") %>%
           tolower() %>% str_replace("regulation of ", "regulation of\n"),
         # Custom editing
         GeneSet = ifelse(GeneSet == "plasma membrane bounded cell projection cytoplasm",
                          "plasma membrane bounded\ncell projection cytoplasm",
                          GeneSet),
         GeneSet = ifelse(GeneSet == "glycosyl compound biosynthetic process",
                          "glycosyl compound\nbiosynthetic process",
                          GeneSet),
         # Reorder based on evidence
         GeneSet = tidytext::reorder_within(GeneSet, neg_log_pval, list_type)) %>%
  ggplot(aes(x = GeneSet, y = neg_log_pval)) +
  geom_bar(stat = "identity", width = 0.5, fill = "darkblue") +
  labs(y = "-log10 adjusted p-value", x = "GO term") +
  coord_flip() +
  tidytext::scale_x_reordered() +
  facet_wrap(~Category, scales = "free", ncol = 2) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14),
        plot.title = element_text(size = 16))

# Save this plot:
save_plot("figures/main/f6_asd_sig_go_bars.jpg",
          asd_pos_esnps_sig_go_bars,
          ncol = 2, nrow = 1)
save_plot("figures/main/f6_asd_sig_go_bars.pdf",
          asd_pos_esnps_sig_go_bars,
          ncol = 2, nrow = 1)


# Create the supplementary figures ----------------------------------------


# ASD Positional + eSNPs GO results with all genes ------------------------

asd_pos_esnps_all_go_bars <- asd_pos_esnps_all_go_table %>%
  filter(Category %in% c("GO_bp", "GO_cc")) %>%
  mutate(Category = fct_relevel(Category, "GO_bp", "GO_cc"),
         Category = fct_recode(Category, `Biological processes` = "GO_bp",
                               `Cellular components` = "GO_cc"),
         neg_log_pval = -log10(adjP)) %>%
  mutate(GeneSet = str_remove(GeneSet, "^GO_") %>% str_replace_all("_", " ") %>%
           tolower() %>% str_replace("regulation of ", "regulation of\n"),
         # Custom editing
         GeneSet = ifelse(GeneSet == "plasma membrane bounded cell projection cytoplasm",
                          "plasma membrane bounded\ncell projection cytoplasm",
                          GeneSet),
         # Reorder based on evidence
         GeneSet = tidytext::reorder_within(GeneSet, neg_log_pval, list_type)) %>%
  ggplot(aes(x = GeneSet, y = neg_log_pval)) +
  geom_bar(stat = "identity", width = 0.5, fill = "darkblue") +
  #geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  labs(y = "-log10 adjusted p-value", x = "GO term") +
  coord_flip() +
  tidytext::scale_x_reordered() +
  facet_wrap(~Category, scales = "free", ncol = 2) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14),
        plot.title = element_text(size = 16))

# Save this plot:
save_plot("figures/suppl/f_asd_all_go_bars.jpg",
          asd_pos_esnps_all_go_bars,
          ncol = 2, nrow = 1)
save_plot("figures/suppl/f_asd_all_go_bars.pdf",
          asd_pos_esnps_all_go_bars,
          ncol = 2, nrow = 1)


# Create the ASD Positional results figures  ------------------------------

# Just stack the two different sources of CC results for Positional assignment
# together since neither display any BP terms
asd_pos_go_cc_bars <- asd_pos_signals_go_table %>%
  bind_rows(asd_pos_all_go_table) %>%
  filter(Category == c("GO_cc")) %>%
  mutate(list_type = fct_relevel(list_type, "signals", "all"),
         list_type = fct_recode(list_type, `Signal genes only` = "signals",
                                `All genes` = "all"),
         neg_log_pval = -log10(adjP)) %>%
  mutate(GeneSet = str_remove(GeneSet, "^GO_") %>% str_replace_all("_", " ") %>%
           tolower() %>% str_replace("regulation of ", "regulation of\n"),
         # Reorder based on evidence
         GeneSet = tidytext::reorder_within(GeneSet, neg_log_pval, list_type)) %>%
  ggplot(aes(x = GeneSet, y = neg_log_pval)) +
  geom_bar(stat = "identity", width = 0.5, fill = "darkblue") +
  labs(y = "-log10 adjusted p-value", x = "Cellular components GO term") +
  coord_flip() +
  tidytext::scale_x_reordered() +
  facet_wrap(~list_type, scales = "free", ncol = 2) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14),
        plot.title = element_text(size = 16))

# Save this plot:
save_plot("figures/suppl/f_asd_pos_go_cc_bars.jpg",
          asd_pos_go_cc_bars,
          ncol = 2, nrow = 1)
save_plot("figures/suppl/f_asd_pos_go_cc_bars.pdf",
          asd_pos_go_cc_bars,
          ncol = 2, nrow = 1)
