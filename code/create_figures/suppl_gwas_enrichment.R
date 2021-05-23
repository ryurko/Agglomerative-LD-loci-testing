# PURPOSE: Generate figure comparing enrichment between ASD, SCZ, and EA SNPs

library(tidyverse)
library(data.table)
library(latex2exp)
library(cowplot)

# Load GWAS data ----------------------------------------------------------

scz_asd_ea_gwas_data <- fread("data/gwas/scz_asd_ea_gwas_results.csv") %>%
  .[hg38_id != "",]

# Create qq-plot comparing enrichment -------------------------------------

set.seed(1000)
gwas_snp_enrichment_plot <- scz_asd_ea_gwas_data %>%
  dplyr::select(hg19_id, asd_p, scz_p, ea_p) %>%
  pivot_longer(asd_p:ea_p,
               names_to = "phenotype",
               values_to = "pval") %>%
  mutate(phenotype = toupper(str_remove(phenotype, "_p")),
         phenotype = fct_relevel(phenotype, "EA", "SCZ", "ASD")) %>%
  group_by(phenotype) %>%
  arrange(pval) %>%
  mutate(neglog_p = -log10(pval),
         exp_neglog_p = -log10(ppoints(n())),
         subset_sample = ifelse(ppoints(n()) > .02,
                                rbernoulli(n(), p = 0.2),
                                1)) %>%
  ungroup() %>%
  filter(subset_sample == 1) %>%
  ggplot(aes(x = exp_neglog_p, y = neglog_p,
             color = phenotype)) +
  ggsci::scale_color_npg() +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.5,
              linetype = "dotted", color = "black") +
  theme_bw() +
  theme(legend.direction = "vertical",
        legend.position = c(.25, .65),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) +
  labs(x = TeX('Expected -log$_{10}$(SNP p-value)'),
       y = TeX('Observed -log$_{10}$(SNP p-value)'),
       color = "Phenotype") +
  guides(color = guide_legend(override.aes = list(size = 5)))


# Save figure -------------------------------------------------------------

save_plot("figures/data/f1_snp_gwas_qqplot.jpg",
          gwas_snp_enrichment_plot)




