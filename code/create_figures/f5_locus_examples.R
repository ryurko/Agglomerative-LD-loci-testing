# PURPOSE: Create clean combined version of locus zoom plots

library(tidyverse)
library(cowplot)

# Source helper functions -------------------------------------------------

source("ld_locus_zoom/kernel_smoothing_fns.R")
source("ld_locus_zoom/get_gene_chart_tiers.R")
source("ld_locus_zoom/get_zoom_x_axis_limits.R")

color_function <- colorRampPalette(c("#440154FF", "darkblue", "#00846b", "darkorange"))

# Load the shiny app data -------------------------------------------------

# File path:
data_file_path <- "ld_locus_zoom/data/positional_esnps"

# Load the gene info:
gene_info_table <- read_csv(paste0(data_file_path, "/gene_info.csv")) %>%
  mutate(Gene = paste0(gene_name, " (", ensembl_id, ")")) %>%
  # Create a size column and sort by largest:
  mutate(gene_size = end - start)  %>%
  # Sort by gene size
  arrange(desc(gene_size))

# Load the LD loci level data ----
all_ld_loci_level_data <-
  read_csv(paste0(data_file_path, "/positional_ld_loci_level_smoothing.csv"))

# Load the gene-level data ----
all_gene_level_data <-
  read_csv(paste0(data_file_path, "/positional_gene_level_smoothing.csv")) %>%
  dplyr::left_join(dplyr::select(gene_info_table, ensembl_id,
                                 Gene),
                   by = "ensembl_id")

# Functional SNPs ----
all_functional_snps_data <-
  read_csv(paste0(data_file_path, "/functional_snp_gene_ld_loci_data.csv")) %>%
  mutate(Gene = paste0(gene_name, " (", ensembl_id, ")"))


# Create plot 5a ----------------------------------------------------------

# Initialize the locus gene data:
chr17_35_gene_table <- gene_info_table %>%
  filter(ld_locus_id == "chr17_35") %>%
  dplyr::select(ensembl_id, ld_locus_id, start, end, strand,
                gene_name, Gene)
# Add the chart tiers:
chr17_35_gene_table <- chr17_35_gene_table %>%
  mutate(gene_chart_tier = get_gene_chart_tiers(chr17_35_gene_table$ensembl_id,
                                                chr17_35_gene_table$start,
                                                chr17_35_gene_table$end),
         # make chart start and end depending on the strand:
         chart_start = ifelse(strand == "+", start, end),
         chart_end = ifelse(strand == "+", end, start))

# Vector of target genes for locus to highlight
chr17_35_target_genes <- c("RAI1", "SREBF1", "TOM1L2", "LRRC48", "ATPAF2", "GID4", "DRG2")

# Indicate in table if gene is one of the targets to highlight:
chr17_35_gene_table <- chr17_35_gene_table %>%
  mutate(is_target = as.numeric(gene_name %in% chr17_35_target_genes))

# Determine the plot axes:
chr17_35_plot_axes <-
  get_zoom_x_axis_limits(
    unique(c(pull(filter(all_ld_loci_level_data, ld_locus_id == "chr17_35"), bp),
             pull(filter(all_functional_snps_data, ld_locus_id == "chr17_35"), bp))),
                         chr17_35_gene_table$start, chr17_35_gene_table$end)

chr17_35_target_gene_names <- chr17_35_target_genes[order(chr17_35_target_genes)]

# List of gene-color vectors to use:
chr17_35_gene_color_list <- list("gene_colors" =
                                    color_function(length(unique(chr17_35_target_gene_names))),
                                  "gene_names" = chr17_35_target_gene_names)
chr17_35_fun_gene_names <- all_functional_snps_data %>%
  filter(ld_locus_id == "chr17_35", gene_name %in% chr17_35_target_genes) %>%
  pull(gene_name) %>%
  unique()
chr17_35_fun_gene_names <- chr17_35_fun_gene_names[order(chr17_35_fun_gene_names)]
chr17_35_gene_color_list[["fun_gene_colors"]] <-
  chr17_35_gene_color_list$gene_colors[
    which(chr17_35_target_gene_names %in% chr17_35_fun_gene_names)
  ]

# Make the plot starting with the smooth signal first
chr17_35_smooth_plot <- all_ld_loci_level_data %>%
  filter(ld_locus_id == "chr17_35") %>%
  ggplot(aes(x = bp)) +
  # Start with rugs of the positional SNPs
  geom_rug(data = filter(all_ld_loci_level_data,
                         ld_locus_id == "chr17_35", is_fake_bp == 0),
           alpha = 0.25, sides = "b") +
  labs(x = "BP", y = "ASD squared z-stats", color = "Gene") +
  # Add eSNPs
  geom_rug(data = {
    filter(all_functional_snps_data, ld_locus_id == "chr17_35",
           gene_name %in% chr17_35_target_genes)
  },
  aes(color = gene_name),
  alpha = 0.25, sides = "t") +
  geom_segment( data = {
    filter(all_functional_snps_data, ld_locus_id == "chr17_35",
           gene_name %in% chr17_35_target_genes)
  },
  aes(xend = bp, y = 0, yend = asd_z_squared, color = gene_name),
  size = 0.5, alpha = 0.25) +
  scale_color_manual(values = chr17_35_gene_color_list$gene_colors,
                     breaks = chr17_35_gene_color_list$gene_names,
                     drop = FALSE) +
  # Add SCZ
  geom_line(data = filter(all_ld_loci_level_data,
                          ld_locus_id == "chr17_35"),
            aes(y = smooth_scz_background,
                group = intra_locus_cluster),
            color = "darkred", size = 1) +
  # Add EA
  geom_line(data = filter(all_ld_loci_level_data,
                          ld_locus_id == "chr17_35"),
            aes(y = smooth_ea_background,
                group = intra_locus_cluster),
            color = "lightblue", size = 1) +
  geom_line(data = filter(all_ld_loci_level_data,
                          ld_locus_id == "chr17_35"),
            aes(y = null_95_percent,
                group = intra_locus_cluster),
            color = "gray", linetype = "dotted",
            size = 1) +
  geom_line(data = filter(all_ld_loci_level_data,
                          ld_locus_id == "chr17_35"),
            aes(y = smooth_asd_zsquared,
                group = intra_locus_cluster),
            color = "black", size = 1) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_x_continuous(limits =
                       # Manually set this for a defined range of interest
                       c(17500000, max(chr17_35_gene_table$end) + 1000)) +
  scale_y_continuous(breaks = seq(0, 10, by = 2))


# Make the gene plot:
chr17_35_gene_loc_plot <- chr17_35_gene_table %>%
  filter(is_target == 0) %>%
  # First the genes that
  ggplot() +
  geom_rect(aes(xmin = start, xmax = end,
                ymin = gene_chart_tier - 0.45,
                ymax = gene_chart_tier + 0.45),
            alpha = .75, fill = "gray", color = "gray") +
  geom_rect(data = filter(chr17_35_gene_table, is_target == 1),
            aes(xmin = start, xmax = end,
                ymin = gene_chart_tier - 0.45,
                ymax = gene_chart_tier + 0.45,
                fill = gene_name, color = gene_name),
            alpha = 1) +
  # Add text labels on top to two genes in particular:
  geom_text(data = filter(chr17_35_gene_table,
                          gene_name %in% c("RAI1", "TOM1L2")) %>%
              dplyr::mutate(gene_id_direction =
                              paste0(gene_name, "\n(", strand, " strand)")),
            aes(x = (start + end) / 2, y = gene_chart_tier,
                label = gene_id_direction), color = "white",
            size = 6) +
  scale_fill_manual(values = chr17_35_gene_color_list$gene_colors,
                    breaks = chr17_35_gene_color_list$gene_names,
                    drop = FALSE) +
  scale_color_manual(values = chr17_35_gene_color_list$gene_colors,
                     breaks = chr17_35_gene_color_list$gene_names,
                     drop = FALSE) +
  theme_minimal() +
  scale_y_reverse() +
  labs(x = "BP", fill = "Gene", color = "Gene") +
  theme(legend.position = c(.15, 0.6),
        legend.direction = "horizontal",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        #legend.title = element_blank(),
        #panel.grid.major.y = element_blank(),
        #panel.grid.minor.y = element_blank(),
        panel.grid = element_blank()) +
  scale_x_continuous(limits =
                       # Manually set this for a defined range of interest
                       c(17500000, max(chr17_35_gene_table$end) + 1000))


# Arrange the figures together:
chr17_35_plot <-
  plot_grid(chr17_35_smooth_plot + theme(legend.position = "none",
                                          plot.margin = unit(c(0,0,0,0), "cm")),
            chr17_35_gene_loc_plot, align = "hv", ncol = 1)


# Create 5b plot of just one gene -----------------------------------------

# Initialize the locus gene data:
chr3_651_gene_table <- gene_info_table %>%
  filter(ld_locus_id == "chr3_651") %>%
  dplyr::select(ensembl_id, ld_locus_id, start, end, strand,
                gene_name, Gene)
# Add the chart tiers:
chr3_651_gene_table <- chr3_651_gene_table %>%
  mutate(gene_chart_tier = get_gene_chart_tiers(chr3_651_gene_table$ensembl_id,
                                                chr3_651_gene_table$start,
                                                chr3_651_gene_table$end),
         # make chart start and end depending on the strand:
         chart_start = ifelse(strand == "+", start, end),
         chart_end = ifelse(strand == "+", end, start))

# Vector of target genes for locus to highlight
chr3_651_target_genes <- c("FOXP1")

# Indicate in table if gene is one of the targets to highlight:
chr3_651_gene_table <- chr3_651_gene_table %>%
  mutate(is_target = as.numeric(gene_name %in% chr3_651_target_genes))

# Determine the plot axes:
chr3_651_plot_axes <-
  get_zoom_x_axis_limits(pull(filter(all_ld_loci_level_data,
                                     ld_locus_id == "chr3_651"),
                              bp),
                         chr3_651_gene_table$start, chr3_651_gene_table$end)

chr3_651_target_gene_names <- chr3_651_target_genes[order(chr3_651_target_genes)]

# List of gene-color vectors to use:
chr3_651_gene_color_list <- list("gene_colors" =
                                   color_function(length(unique(chr3_651_target_gene_names))),
                                 "gene_names" = chr3_651_target_gene_names)
chr3_651_fun_gene_names <- all_functional_snps_data %>%
  filter(ld_locus_id == "chr3_651", gene_name %in% chr3_651_target_genes) %>%
  pull(gene_name) %>%
  unique()
chr3_651_fun_gene_names <- chr3_651_fun_gene_names[order(chr3_651_fun_gene_names)]
chr3_651_gene_color_list[["fun_gene_colors"]] <-
  chr3_651_gene_color_list$gene_colors[
    which(chr3_651_target_gene_names %in% chr3_651_fun_gene_names)
  ]

# Make the plot starting with the smooth signal first
chr3_651_smooth_plot <- all_ld_loci_level_data %>%
  filter(ld_locus_id == "chr3_651") %>%
  ggplot(aes(x = bp)) +
  # Start with rugs of the positional SNPs
  geom_rug(data = filter(all_ld_loci_level_data,
                         ld_locus_id == "chr3_651", is_fake_bp == 0),
           alpha = 0.25, sides = "b") +
  labs(x = "BP", y = "ASD squared z-stats", color = "Gene") +
  # Add eSNPs
  geom_rug(data = {
    filter(all_functional_snps_data, ld_locus_id == "chr3_651",
           gene_name %in% chr3_651_target_genes)
  },
  aes(color = gene_name),
  alpha = 0.25, sides = "t") +
  geom_segment( data = {
    filter(all_functional_snps_data, ld_locus_id == "chr3_651",
           gene_name %in% chr3_651_target_genes)
  },
  aes(xend = bp, y = 0, yend = asd_z_squared, color = gene_name),
  size = 0.5, alpha = 0.25) +
  scale_color_manual(values = chr3_651_gene_color_list$gene_colors,
                     breaks = chr3_651_gene_color_list$gene_names,
                     drop = FALSE) +
  # Add SCZ
  geom_line(data = filter(all_ld_loci_level_data,
                          ld_locus_id == "chr3_651"),
            aes(y = smooth_scz_background),
            color = "darkred", size = 1) +
  # Add EA
  geom_line(data = filter(all_ld_loci_level_data,
                          ld_locus_id == "chr3_651"),
            aes(y = smooth_ea_background),
            color = "lightblue", size = 1) +
  geom_line(data = filter(all_ld_loci_level_data,
                          ld_locus_id == "chr3_651"),
            aes(y = null_95_percent),
            color = "gray", linetype = "dotted",
            size = 1) +
  geom_line(data = filter(all_ld_loci_level_data,
                          ld_locus_id == "chr3_651"),
            aes(y = smooth_asd_zsquared),
            color = "black", size = 1) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_x_continuous(limits =
                       c(chr3_651_plot_axes$chart_x_min - 1000,
                         chr3_651_plot_axes$chart_x_max + 1000)) +
  scale_y_continuous(breaks = seq(0, 10, by = 2))


# Make the gene plot:
chr3_651_gene_loc_plot <- chr3_651_gene_table %>%
  dplyr::mutate(gene_id_direction =
                  paste0(gene_name, "\n(", strand, " strand)")) %>%
  ggplot() +
  geom_rect(aes(xmin = start, xmax = end,
                ymin = gene_chart_tier - 0.45,
                ymax = gene_chart_tier + 0.45,
                fill = gene_name, color = gene_name),
            alpha = 1) +
  # Add text labels on top to two genes in particular:
  geom_text(aes(x = (start + end) / 2, y = gene_chart_tier,
                label = gene_id_direction), color = "white",
            size = 9) +
  scale_fill_manual(values = chr3_651_gene_color_list$gene_colors,
                    breaks = chr3_651_gene_color_list$gene_names,
                    drop = FALSE) +
  scale_color_manual(values = chr3_651_gene_color_list$gene_colors,
                     breaks = chr3_651_gene_color_list$gene_names,
                     drop = FALSE) +
  theme_minimal() +
  scale_y_reverse() +
  labs(x = "BP", fill = "Gene", color = "Gene") +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        #legend.title = element_blank(),
        #panel.grid.major.y = element_blank(),
        #panel.grid.minor.y = element_blank(),
        panel.grid = element_blank()) +
  scale_x_continuous(limits = c(chr3_651_plot_axes$chart_x_min - 1000,
                                chr3_651_plot_axes$chart_x_max + 1000))


# Arrange the figures together:
chr3_651_plot <-
  plot_grid(chr3_651_smooth_plot + theme(plot.margin = unit(c(0,0,0,0), "cm")),
            chr3_651_gene_loc_plot, align = "hv", ncol = 1)



# Create focused version of 5c --------------------------------------------

# Initialize the locus gene data:
chr17_1_gene_table <- gene_info_table %>%
  filter(ld_locus_id == "chr17_1") %>%
  dplyr::select(ensembl_id, ld_locus_id, start, end, strand,
                gene_name, Gene)
# Add the chart tiers:
chr17_1_gene_table <- chr17_1_gene_table %>%
  mutate(gene_chart_tier = get_gene_chart_tiers(chr17_1_gene_table$ensembl_id,
                                                chr17_1_gene_table$start,
                                                chr17_1_gene_table$end),
         # make chart start and end depending on the strand:
         chart_start = ifelse(strand == "+", start, end),
         chart_end = ifelse(strand == "+", end, start))

# Vector of target genes for locus to highlight
chr17_1_target_genes <- c("MAPT", "KANSL1", "LRRC37A", "ARL17A", "NSF", "WNT3")

# Indicate in table if gene is one of the targets to highlight:
chr17_1_gene_table <- chr17_1_gene_table %>%
  mutate(is_target = as.numeric(gene_name %in% chr17_1_target_genes))

# Determine the plot axes:
chr17_1_plot_axes <-
  get_zoom_x_axis_limits(pull(filter(all_ld_loci_level_data,
                                     ld_locus_id == "chr17_1"),
                              bp),
                         chr17_1_gene_table$start, chr17_1_gene_table$end)

chr17_1_target_gene_names <- chr17_1_target_genes[order(chr17_1_target_genes)]

# List of gene-color vectors to use:
chr17_1_gene_color_list <- list("gene_colors" =
                                    color_function(length(unique(chr17_1_target_gene_names))),
                                  "gene_names" = chr17_1_target_gene_names)
chr17_1_fun_gene_names <- all_functional_snps_data %>%
  filter(ld_locus_id == "chr17_1", gene_name %in% chr17_1_target_genes) %>%
  pull(gene_name) %>%
  unique()
chr17_1_fun_gene_names <- chr17_1_fun_gene_names[order(chr17_1_fun_gene_names)]
chr17_1_gene_color_list[["fun_gene_colors"]] <-
  chr17_1_gene_color_list$gene_colors[
    which(chr17_1_target_gene_names %in% chr17_1_fun_gene_names)
  ]

# Make the plot starting with the smooth signal first
chr17_1_smooth_plot <- all_ld_loci_level_data %>%
  filter(ld_locus_id == "chr17_1") %>%
  ggplot(aes(x = bp)) +
  # Start with rugs of the positional SNPs
  geom_rug(data = filter(all_ld_loci_level_data,
                         ld_locus_id == "chr17_1", is_fake_bp == 0),
           alpha = 0.25, sides = "b") +
  labs(x = "BP", y = "ASD squared z-stats", color = "Gene") +
  # Add eSNPs
  geom_rug(data = {
    filter(all_functional_snps_data, ld_locus_id == "chr17_1",
           gene_name %in% chr17_1_target_genes)
  },
  aes(color = gene_name),
  alpha = 0.25, sides = "t") +
  geom_segment( data = {
    filter(all_functional_snps_data, ld_locus_id == "chr17_1",
           gene_name %in% chr17_1_target_genes)
  },
  aes(xend = bp, y = 0, yend = asd_z_squared, color = gene_name),
  size = 0.5, alpha = 0.25) +
  scale_color_manual(values = chr17_1_gene_color_list$gene_colors,
                     breaks = chr17_1_gene_color_list$gene_names,
                     drop = FALSE) +
  # Add SCZ
  geom_line(data = filter(all_ld_loci_level_data,
                          ld_locus_id == "chr17_1"),
            aes(y = smooth_scz_background,
                group = intra_locus_cluster),
            color = "darkred", size = 1) +
  # Add EA
  geom_line(data = filter(all_ld_loci_level_data,
                          ld_locus_id == "chr17_1"),
            aes(y = smooth_ea_background,
                group = intra_locus_cluster),
            color = "lightblue", size = 1) +
  geom_line(data = filter(all_ld_loci_level_data,
                          ld_locus_id == "chr17_1"),
            aes(y = null_95_percent,
                group = intra_locus_cluster),
            color = "gray", linetype = "dotted",
            size = 1) +
  geom_line(data = filter(all_ld_loci_level_data,
                          ld_locus_id == "chr17_1"),
            aes(y = smooth_asd_zsquared,
                group = intra_locus_cluster),
            color = "black", size = 1) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_x_continuous(limits =
                       c(min(chr17_1_gene_table$start),
                         max(chr17_1_gene_table$end)))


# Make the gene plot:
chr17_1_gene_loc_plot <- chr17_1_gene_table %>%
  filter(is_target == 0) %>%
  # First the genes that
  ggplot() +
  geom_rect(aes(xmin = start, xmax = end,
                ymin = gene_chart_tier - 0.45,
                ymax = gene_chart_tier + 0.45),
            alpha = .75, fill = "gray", color = "gray") +
  geom_rect(data = filter(chr17_1_gene_table, is_target == 1),
            aes(xmin = start, xmax = end,
                ymin = gene_chart_tier - 0.45,
                ymax = gene_chart_tier + 0.45,
                fill = gene_name, color = gene_name),
            alpha = 1) +
  # Add text labels on top to two genes in particular:
  geom_text(data = filter(chr17_1_gene_table,
                          gene_name %in% c("MAPT", "KANSL1")) %>%
              dplyr::mutate(gene_id_direction =
                              paste0(gene_name, "\n(", strand, " strand)")),
            aes(x = (start + end) / 2, y = gene_chart_tier,
                label = gene_id_direction), color = "white",
            size = 6) +
  scale_fill_manual(values = chr17_1_gene_color_list$gene_colors,
                    breaks = chr17_1_gene_color_list$gene_names,
                    drop = FALSE) +
  scale_color_manual(values = chr17_1_gene_color_list$gene_colors,
                     breaks = chr17_1_gene_color_list$gene_names,
                     drop = FALSE) +
  theme_minimal() +
  scale_y_reverse() +
  labs(x = "BP", fill = "Gene", color = "Gene") +
  theme(legend.position = c(.2, .2),
        legend.direction = "horizontal",
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        #legend.title = element_blank(),
        #panel.grid.major.y = element_blank(),
        #panel.grid.minor.y = element_blank(),
        panel.grid = element_blank()) +
  scale_x_continuous(limits =
                       c(min(chr17_1_gene_table$start),
                         max(chr17_1_gene_table$end)))

# Arrange the figures together:
chr17_1_plot <-
  plot_grid(chr17_1_smooth_plot + theme(legend.position = "none",
                                          plot.margin = unit(c(0,0,0,0), "cm")),
            chr17_1_gene_loc_plot, align = "hv", ncol = 1)


# Arrange all three together ----------------------------------------------


fig_five_grid <-
  plot_grid(NULL, chr17_35_plot, NULL, chr3_651_plot, NULL, chr17_1_plot,
            labels = c("", "A", "", "B", "", "C"), label_fontface = "plain",
            rel_heights = c(0.1, 1, 0.1, 1, 0.1, 1),
            hjust = 0, vjust = -0.5, label_size = 36, align = "hv",
            nrow = 6)
save_plot("figures/main/f5_zoom_grid_ex.jpg",
          fig_five_grid, ncol = 1, nrow = 6, base_asp = 4)
save_plot("figures/main/f5_zoom_grid_ex.pdf",
          fig_five_grid, ncol = 1, nrow = 6, base_asp = 4)


