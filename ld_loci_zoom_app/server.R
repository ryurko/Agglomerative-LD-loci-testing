# PURPOSE: Create the display for the LD loci zoom app, separated by the type
#          of content: (1) our paper's ASD results and (2) custom results
#          uploaded by the user. There are similar steps between each set with
#          the key difference being the need to generate the kernel smoothing
#          results for the user uploaded files.

library(shiny)
library(shinyjs)
library(shinyalert)
options(shiny.maxRequestSize = 30*1024^2)
library(tidyverse)
library(data.table)
library(plotly)

# TODO: Source the helper functions for simplifying the steps in this code
source("kernel_smoothing_fns.R")
source("get_gene_chart_tiers.R")
source("get_zoom_x_axis_limits.R")

# Initialize the color function:
color_function <- colorRampPalette(c("#440154FF", "darkblue", "#00846b", "darkorange"))

# Define server logic required to draw a histogram
shinyServer(function(input, output) {

# Construct the output for the paper's ASD results ------------------------

    # Load the datasets given the SNP to gene assignment type, by first initializing
    # the file path for the data
    data_file_path <- reactive({
        paste0("data/",
               janitor::make_clean_names(str_remove(tolower(input$assign_type), "-")))
    })

    # Load the gene info:
    gene_info_table <- reactive({
        req(data_file_path())
        read_csv(paste0(data_file_path(), "/gene_info.csv")) %>%
            mutate(Gene = paste0(gene_name, " (", ensembl_id, ")")) %>%
            # Create a size column and sort by largest:
            mutate(gene_size = end - start)  %>%
            # Sort by gene size
            arrange(desc(gene_size))
    })

    # Load the LD loci level data ----
    all_ld_loci_level_data <- reactive({
        req(data_file_path())
        read_csv(paste0(data_file_path(), "/positional_ld_loci_level_smoothing.csv"))
    })

    # Load the gene-level data ----
    all_gene_level_data <- reactive({
        read_csv(paste0(data_file_path(), "/positional_gene_level_smoothing.csv")) %>%
            dplyr::left_join(dplyr::select(gene_info_table(), ensembl_id,
                                           Gene),
                             by = "ensembl_id")
    })

    # Functional SNPs ----
    all_functional_snps_data <- reactive({
        if (input$assign_type == "Positional") {
            tibble(ld_loci_id = NULL)
        } else {
            read_csv(paste0(data_file_path(), "/functional_snp_gene_ld_loci_data.csv")) %>%
                mutate(Gene = paste0(gene_name, " (", ensembl_id, ")"))
        }
    })

    # Select CHR ----
    output$chr_selection <-
        renderUI({
            req(all_ld_loci_level_data())
            selectizeInput(inputId = "chr",
                           label = "Select CHR:",
                           choices = sort(unique(all_ld_loci_level_data()$chr)),
                           selected = sort(unique(all_ld_loci_level_data()$chr))[1])
        })

    # Render input for selecting LD loci id in selected chromosome
    output$ld_loci_selection <-
        renderUI({
            req(all_ld_loci_level_data(), input$chr)
            candidate_ld_loci_ids <- all_ld_loci_level_data() %>%
                    filter(chr == input$chr) %>%
                    pull(ld_loci_id)  %>%
                    unique()
            selectizeInput(inputId = "ld_loci",
                           label = "Select LD loci ID:",
                           choices = c("choose" = "", candidate_ld_loci_ids),
                           selected = candidate_ld_loci_ids[1])
        })

    # Filter the respective LD loci, gene, and functional datasets:
    ld_loci_level_smooth_data <- reactive({
        req(input$ld_loci, all_ld_loci_level_data())
        all_ld_loci_level_data() %>%
            filter(ld_loci_id == input$ld_loci) %>%
            dplyr::rename(`Smooth SCZ z-squared` = smooth_scz_background,
                          `Smooth EA z-squared` = smooth_ea_background,
                          `Smooth ASD z-squared` = smooth_asd_zsquared,
                          `Null 95th percentile` = null_95_percent,
                          BP = bp)
    })
    gene_level_smooth_data <- reactive({
        req(input$ld_loci, all_gene_level_data())
        all_gene_level_data() %>%
            filter(ld_loci_id == input$ld_loci) %>%
            dplyr::rename(`Smooth ASD z-squared` = smooth_asd_zsquared,
                          BP = bp)
    })
    # Get the functional data:
    functional_snps_data <- reactive({
        req(input$ld_loci, all_functional_snps_data())
        if (input$assign_type == "Positional") {
            all_functional_snps_data()
        } else {
            all_functional_snps_data() %>%
                filter(ld_loci_id == input$ld_loci) %>%
                dplyr::rename(BP = bp, `ASD z-squared` = asd_z_squared)
        }
    })

    # Load the SNP-gene data:
    snp_gene_pairs <- reactive({
        req(input$ld_loci, data_file_path())
        read_csv(paste0(data_file_path(), "/snp_gene_table.csv")) %>%
            filter(ld_loci_id == input$ld_loci)
    })


    # Create the appropriate gene-location table as well:
    gene_loc_table_display <- reactive({
        req(gene_info_table(), input$ld_loci)

        gene_loc_table <- gene_info_table() %>%
            filter(ld_loci_id == input$ld_loci) %>%
            dplyr::select(ensembl_id, ld_loci_id, start, end, strand,
                          gene_name, Gene)

        gene_loc_table %>%
            # Add the chart tier:
            mutate(gene_chart_tier = get_gene_chart_tiers(gene_loc_table$ensembl_id,
                                                          gene_loc_table$start,
                                                          gene_loc_table$end),
                   # make chart start and end depending on the strand:
                   chart_start = ifelse(strand == "+", start, end),
                   chart_end = ifelse(strand == "+", end, start))

    })

    # List of gene-color vectors to use:
    gene_color_list <- reactive({
        req(functional_snps_data(), gene_loc_table_display())

        gene_vec <- unique(gene_loc_table_display()$Gene)
        ordered_gene_vec <- gene_vec[order(gene_vec)]

        gene_colors <- color_function(length(unique(ordered_gene_vec)))
        colors_result <- list("gene_colors" = gene_colors,
                              "gene_names" = ordered_gene_vec)

        # Now if there are any eSNPs - get their gene colors:
        if (nrow(functional_snps_data()) > 0) {

            fun_gene_vec <- functional_snps_data() %>%
                dplyr::pull(Gene) %>%
                unique()
            ordered_fun_gene_vec <- fun_gene_vec[order(fun_gene_vec)]
            colors_result[["fun_gene_colors"]] <-
                gene_colors[which(ordered_gene_vec %in% ordered_fun_gene_vec)]

        }

        colors_result

    })

    # Determine the plot axes:
    plot_axes <- reactive({
        req(gene_loc_table_display(), ld_loci_level_smooth_data())
        get_zoom_x_axis_limits(ld_loci_level_smooth_data()$BP,
                               gene_loc_table_display()$start,
                               gene_loc_table_display()$end)
    })

    output$gene_table <- DT::renderDataTable(server = FALSE, {
        req(gene_loc_table_display(), input$ld_loci)
        gene_loc_table_display() %>%
            dplyr::select(ensembl_id, gene_name, start, end, strand, ld_loci_id) %>%
            mutate(gwas_catalog = paste0("https://www.ebi.ac.uk/gwas/genes/",
                                         gene_name),
                   gwas_catalog = paste0("<a href='", gwas_catalog, "'>",
                                         gwas_catalog, "</a>")) %>%
            dplyr::select(-ld_loci_id) %>%
            DT::datatable(escape = FALSE,
                          extensions = c("Buttons"),
                          options = list(
                              dom = 'Bfrtip',
                              buttons = list(
                                  list(extend = 'collection',
                                       buttons = list(
                                           list(extend = 'csv',
                                                filename = paste0(input$ld_loci, "_genes")),
                                           list(extend = 'excel',
                                                filename = paste0(input$ld_loci, "_genes"))),
                                       text = "Download gene table"
                                  )
                              )
                          ),
                          rownames = FALSE)
    })

    output$snp_gene_table <- DT::renderDataTable(server = FALSE, {
        req(snp_gene_pairs(), input$ld_loci)
        snp_gene_pairs() %>%
            mutate(gwas_catalog = paste0("<a href='", gwas_catalog, "'>",
                                         gwas_catalog, "</a>")) %>%
            dplyr::select(-ld_loci_id) %>%
            DT::datatable(escape = FALSE,
                          extensions = c("Buttons", "KeyTable", "FixedColumns"),
                          options = list(
                              dom = 'Bfrtip',
                              scrollX = TRUE,
                              fixedColumns = list(leftColumns = 3),
                              keys = TRUE,
                              buttons = list(
                                  list(extend = 'collection',
                                       buttons = list(
                                           list(extend = 'csv',
                                                filename = paste0(input$ld_loci, "_snp_gene_pairs")),
                                           list(extend = 'excel',
                                                filename = paste0(input$ld_loci, "_snp_gene_pairs"))),
                                       text = "Download SNP-gene table"
                                  )
                              )
                          ),
                          rownames = FALSE)
    })

    output$zoomPlot <- renderPlotly({
        #browser()
        req(gene_loc_table_display(), plot_axes())

        # While the main dataset for the plot is ld_loci_level_smooth_data() -
        # the input determines what is added to the display, but we can start
        # with the base layer regardless of the input:
        base_snp_plot <- ld_loci_level_smooth_data() %>%
            ggplot(aes(x = BP)) +
            # Start with rugs of the positional SNPs
            geom_rug(data = filter(ld_loci_level_smooth_data(), is_fake_bp == 0),
                     alpha = 0.25, sides = "b") +
            labs(x = "BP", y = "ASD squared z-stats",
                 color = "Gene") +
            theme_bw()

        # Add the eQTLs to the plot separately?
        if (input$display_esnps) {

            # First check that there are eSNPs:
            if (nrow(functional_snps_data()) > 0) {

                base_snp_plot <- base_snp_plot +
                    geom_rug(data = {
                        functional_snps_data()
                    },
                    aes(color = Gene),
                    alpha = 0.25, sides = "t") +
                    geom_segment(
                        data = {
                            functional_snps_data()
                        }, aes(xend = BP, y = 0, yend = `ASD z-squared`,
                               color = Gene),
                        size = 0.5, alpha = 0.25)
            }

        }

        # Next determine if there is background:
        if ("scz" %in% input$background_pheno) {

            base_snp_plot <- base_snp_plot +
                geom_line(data = ld_loci_level_smooth_data(),
                          aes(y = `Smooth SCZ z-squared`),
                          color = "darkred", alpha = 0.8)
        }

        if ("ea" %in% input$background_pheno) {

            base_snp_plot <- base_snp_plot +
                geom_line(data = ld_loci_level_smooth_data(),
                          aes(y = `Smooth EA z-squared`),
                          color = "lightblue", alpha = 0.8)
        }

        # Now display the null results:
        base_snp_plot <- base_snp_plot +
            geom_line(data = ld_loci_level_smooth_data(),
                      aes(y = `Null 95th percentile`), alpha = .8,
                      color = "gray", linetype = "dotted",
                      size = 1)

        # If LD loci level is one of the choices then display it:
        if ("ld_loci" %in% input$smooth_level) {
            base_snp_plot <- base_snp_plot +
                geom_line(data = ld_loci_level_smooth_data(),
                          aes(y = `Smooth ASD z-squared`), alpha = .8,
                          color = "black", size = 1)
        }

        # Add gene-level smoothing?
        if ("gene" %in% input$smooth_level) {
            #browser()
            base_snp_plot <- base_snp_plot +
                geom_line(data = gene_level_smooth_data(), alpha = 0.8,
                          aes(y = `Smooth ASD z-squared`,
                              group = Gene, color = Gene))
        }

        # Finally add the color scale that should be used regardless of any settings
        base_snp_plot <- base_snp_plot +
            scale_color_manual(values = gene_color_list()$gene_colors,
                               breaks = gene_color_list()$gene_names,
                               drop = FALSE)

        # Next create the gene display below:
        gene_loc_plot <- gene_loc_table_display() %>%
            dplyr::mutate(gene_id_direction =
                              paste0(Gene, "\nstrand: ", strand,
                                     ifelse(strand == "+", " (left to right)",
                                            " (right to left)"))) %>%
            ggplot(aes(xmin = start, xmax = end,
                       ymin = gene_chart_tier - 0.45,
                       ymax = gene_chart_tier + 0.45,
                       fill = Gene, color = Gene,
                       text = gene_id_direction)) +
            geom_rect(alpha = 1) +
            scale_fill_manual(values = gene_color_list()$gene_colors,
                              breaks = gene_color_list()$gene_names,
                              drop = FALSE) +
            scale_color_manual(values = gene_color_list()$gene_colors,
                               breaks = gene_color_list()$gene_names,
                               drop = FALSE) +
            theme_minimal() +
            scale_y_reverse() +
            theme(legend.position = "bottom",
                  axis.text.y = element_blank()) +
            labs(x = "BP", fill = "Gene", color = "Gene")

        # First create the ggplotly object for the SNP plot, then modify its
        # legend layers before joining with the gene table in a subplot display:
        base_snp_plotly <- ggplotly(
            base_snp_plot +
                theme(legend.title = element_blank()) +
                scale_x_continuous(limits =
                                       c(plot_axes()$chart_x_min - 1000,
                                         plot_axes()$chart_x_max + 1000)),
            tooltip = c("BP", "Smooth ASD z-squared", "Smooth SCZ z-squared",
                        "Smooth EA z-squared", "Gene", "ASD z-squared",
                        "Null 95th percentile")
        )

        # If there are any functional SNPs  displayed separately then update the legend
        # to handle these genes:
        if (input$display_esnps) {
            # First check that there are functional SNPs:
            if (nrow(functional_snps_data()) > 0) {
                base_snp_plotly <- base_snp_plotly %>%
                    style(showlegend = FALSE,
                          traces = c(2:((length(gene_color_list()$fun_gene_colors) * 2) + 1)))
            }
        }

        # Next if gene-level then remove the final layer of colors from the legend:
        if ("gene" %in% input$smooth_level) {

            # These are the final layer so how many layers are there currently?
            n_base_snp_traces <- length(base_snp_plotly$x$data)
            # How many genes are displayed:
            n_gene_lines <- gene_level_smooth_data() %>%
                pull(Gene) %>%
                unique() %>%
                length()

            # Now update the trace
            start_gene_trace <- n_base_snp_traces - n_gene_lines + 1
            base_snp_plotly <- base_snp_plotly %>%
                style(showlegend = FALSE,
                      traces = c(start_gene_trace:n_base_snp_traces))
        }

        # Make the gene_loc_plotly object:
        gene_loc_plotly <- ggplotly(gene_loc_plot +
                                        theme(legend.title = element_blank(),
                                              panel.grid.major.y = element_blank(),
                                              panel.grid.minor.y = element_blank()) +
                                        scale_x_continuous(limits = c(plot_axes()$chart_x_min - 1000,
                                                                      plot_axes()$chart_x_max + 1000)),
                                    tooltip = "text")


        # Create the plotly subplot display
        subplot_display <-
            subplot(base_snp_plotly,
                    gene_loc_plotly,
                    shareX = TRUE, shareY = FALSE, titleY = TRUE, nrows = 2)


        subplot_display %>%
            layout(legend = list(orientation = "v", x = 100, y = 0.5)) %>%
            config(
                toImageButtonOptions = list(
                    format = "svg",
                    filename = paste0(input$ld_loci, "_zoom_plot"),
                    width = 1200,
                    height = 900
                )
            )


    })



# Construct the output for the custom results -----------------------------


# First make the table of the upload genes --------------------------------

        # Load the gene info:
    custom_gene_info_table <- reactive({
        req(input$custom_genes)

        req_gene_table_cols <- c("gene_id", "start", "end", "chr", "loci_id")
        custom_gene_data <- read_csv(input$custom_genes$datapath)

        if (all(req_gene_table_cols %in% colnames(custom_gene_data))) {
            custom_gene_data %>%
                # Create a size column and sort by largest:
                mutate(gene_size = end - start) %>%
                return()
        } else {
            shinyalert("Column Error",
                       paste0("Your gene table is missing the following columns: ",
                              paste0(dplyr::setdiff(req_gene_table_cols,
                                                    colnames(custom_gene_data)),
                                     collapse = ", ")),
                       type = "error")
            returnValue()
        }

    })

# Generate the smoothing results with the upload data ---------------------

    # First load the custom SNP-gene dataset
    custom_snp_gene_table <- reactive({
        req(input$custom_snp_gene, input$custom_genes)

        req_snp_gene_table_cols <- c("snp_id", "gene_id", "bp", "snp_signal")
        custom_snp_gene_data <- read_csv(input$custom_snp_gene$datapath)

        if (all(req_snp_gene_table_cols %in% colnames(custom_snp_gene_data))) {
            custom_snp_gene_data  %>%
                # Join the loci ID from the gene table:
                dplyr::left_join(dplyr::select(custom_gene_info_table(),
                                               gene_id, loci_id),
                                 by = "gene_id") %>%
                return()
        } else {
            shinyalert("Column Error",
                       paste0("Your SNP-gene table is missing the following columns: ",
                              paste0(dplyr::setdiff(req_snp_gene_table_cols,
                                                    colnames(custom_snp_gene_data)),
                                     collapse = ", ")),
                       type = "error")
            returnValue()
        }

    })

    # Now create the SNP-loci dataset for the smoothing
    custom_snp_loci_table <- reactive({
        custom_snp_gene_table() %>%
            dplyr::select(-gene_id) %>%
            dplyr::distinct()
    })

    # Generate the loci-level GCV bandwidths
    custom_loci_gcv <- reactive({
        get_gcv_bw(custom_snp_loci_table(), min_n_snps = input$min_n_snps)
    })

    # Next generate the loci level kernel smoothing results
    custom_snp_loci_smoothing <- reactive({
        get_loci_level_smoothing(custom_snp_loci_table(),
                                   custom_loci_gcv(),
                                   n_points = input$loci_interp)
    })

    # Generate the gene level kernel smoothing results
    custom_snp_gene_smoothing <- reactive({
        get_gene_level_smoothing(custom_snp_gene_table(),
                                 custom_loci_gcv(),
                                 n_points = input$gene_interp,
                                 min_n_snps = input$min_n_snps)
    })


# Load functional and reference data if provided --------------------------

    custom_fun_snp_gene_table <- reactive({
        req(input$custom_fun_snp_gene, input$custom_genes)

        req_snp_gene_table_cols <- c("snp_id", "gene_id", "bp", "snp_signal")
        custom_fun_snp_gene_data <- read_csv(input$custom_fun_snp_gene$datapath)

        if (all(req_snp_gene_table_cols %in% colnames(custom_fun_snp_gene_data))) {
            custom_fun_snp_gene_data  %>%
                # Join the loci ID from the gene table:
                dplyr::left_join(dplyr::select(custom_gene_info_table(),
                                               gene_id, loci_id),
                                 by = "gene_id") %>%
                # Add indicator column to denote to these are functional
                mutate(is_functional = 1, is_positional = 0) %>%
                return()
        } else {
            shinyalert("Column Error",
                       paste0("Your functional SNP-gene table is missing the following columns: ",
                              paste0(dplyr::setdiff(req_snp_gene_table_cols,
                                                    colnames(custom_fun_snp_gene_data)),
                                     collapse = ", ")),
                       type = "error")
            returnValue()
        }

    })

    custom_reference_table <- reactive({
        req(input$custom_ref_points)

        req_ref_table_cols <- c("bp", "reference_signal", "loci_id")
        custom_ref_data <- read_csv(input$custom_ref_points$datapath)

        if (all(req_ref_table_cols %in% colnames(custom_ref_data))) {
            custom_ref_data %>%
                return()
        } else {
            shinyalert("Column Error",
                       paste0("Your reference signal dataset is missing the following columns: ",
                              paste0(dplyr::setdiff(req_ref_table_cols,
                                                    colnames(custom_ref_data)),
                                     collapse = ", ")),
                       type = "error")
            returnValue()
        }


    })

# Select which custom loci to display -----------------------------------

    # Select CHR ----
    output$custom_chr_selection <-
        renderUI({
            #req(input$custom_snp_gene, input$custom_genes)
            if(!is.null(input$custom_snp_gene) & !is.null(input$custom_genes)) {
                chr_options <- sort(unique(custom_gene_info_table()$chr))
            } else {
                chr_options <- ""
            }

            selectizeInput(inputId = "custom_chr",
                           label = "Select CHR:",
                           choices = chr_options,
                           selected = chr_options[1])
        })

    # Render input for selecting LD loci id in selected chromosome
    output$custom_loci_selection <-
        renderUI({

            if(!is.null(input$custom_snp_gene) & !is.null(input$custom_genes)) {
                loci_options <- custom_gene_info_table() %>%
                    filter(chr == input$custom_chr) %>%
                    pull(loci_id)  %>%
                    unique()
            } else {
                loci_options <- ""
            }

            selectizeInput(inputId = "custom_loci",
                           label = "Select loci ID:",
                           choices = loci_options,
                           selected = loci_options[1])

        })

# Make reactive expressions to display warnings input ---------------------

    custom_fun_status <- reactive({
        req_snp_gene_table_cols <- c("snp_id", "gene_id", "bp", "snp_signal")

        if (input$custom_display_fun_snps) {

            if (!is.null(input$custom_fun_snp_gene)) {

                if (all(req_snp_gene_table_cols %in% colnames(custom_fun_snp_gene_table()))) {
                    TRUE

                } else {
                    shinyalert("Column Error",
                               paste0("Your functional SNP-gene table is missing the following columns: ",
                                      paste0(dplyr::setdiff(req_snp_gene_table_cols,
                                                            colnames(custom_fun_snp_gene_table())),
                                             collapse = ", ")),
                               type = "error")
                    FALSE
                }


            } else {
                shinyalert("Missing functional SNP-gene data",
                           "You need to upload a dataset of functional SNP-gene pairs to display",
                           type = "warning")
                FALSE
            }
        } else {
            FALSE
        }
    })

    custom_ref_status <- reactive({
        req_ref_table_cols <- c("bp", "reference_signal", "loci_id")
        if (input$custom_display_reference) {

            if (!is.null(input$custom_ref_points)) {

                if (all(req_ref_table_cols %in% colnames(custom_reference_table()))) {
                    TRUE
                } else {
                    shinyalert("Column Error",
                               paste0("Your reference signal dataset is missing the following columns: ",
                                      paste0(dplyr::setdiff(req_ref_table_cols,
                                                            colnames(custom_reference_table())),
                                             collapse = ", ")),
                               type = "error")
                    FALSE
                }

            } else {
                shinyalert("Missing reference signal data",
                           "You need to upload a reference signal dataset to display",
                           type = "warning")
                FALSE
            }
        } else {
            FALSE
        }
    })


# Filter to the selected data for display ---------------------------------

    # Filter the respective LD loci, gene, and functional datasets:
    custom_loci_smooth_data <- reactive({
        custom_snp_loci_smoothing() %>%
            filter(loci_id == input$custom_loci) %>%
            dplyr::rename(`Smooth signal` = smooth_signal,
                          BP = bp)
    })
    custom_gene_smooth_data <- reactive({
        custom_snp_gene_smoothing() %>%
            filter(loci_id == input$custom_loci) %>%
            dplyr::rename(`Smooth signal` = smooth_signal,
                          BP = bp)
    })

    # Get the functional data:
    custom_loci_fun_snp_data <- reactive({
        req_snp_gene_table_cols <- c("snp_id", "gene_id", "bp", "snp_signal")

        if (all(req_snp_gene_table_cols %in% colnames(custom_fun_snp_gene_table()))) {
            custom_fun_snp_gene_table() %>%
                filter(loci_id == input$custom_loci) %>%
                dplyr::rename(BP = bp,
                              `Functional SNP signal` = snp_signal) %>%
                return()

        } else {
            shinyalert("Column Error",
                       paste0("Your functional SNP-gene table is missing the following columns: ",
                              paste0(dplyr::setdiff(req_snp_gene_table_cols,
                                                    colnames(custom_fun_snp_gene_table())),
                                     collapse = ", ")),
                       type = "error")
            returnValue()
        }

        # custom_fun_snp_gene_table() %>%
        #     filter(loci_id == input$custom_loci) %>%
        #     dplyr::rename(BP = bp,
        #                   `Functional SNP signal` = snp_signal)

    })

    # Get the reference data:
    custom_loci_ref_data <- reactive({
        custom_reference_table() %>%
            filter(loci_id == input$custom_loci) %>%
            dplyr::rename(BP = bp,
                          `Reference signal` = reference_signal)
    })

    custom_loci_gene_data <- reactive({
        result_data <- custom_gene_info_table() %>%
            filter(loci_id == input$custom_loci)  %>%
            # Sort by gene size
            arrange(desc(gene_size))
        result_data %>%
            # Add the chart tier:
            mutate(gene_chart_tier = get_gene_chart_tiers(result_data$gene_id,
                                                          result_data$start,
                                                          result_data$end))
    })

    # List of gene-color vectors to use:
    custom_gene_color_list <- reactive({
        custom_gene_vec <- unique(custom_loci_gene_data()$gene_id)
        custom_ordered_gene_vec <- custom_gene_vec[order(custom_gene_vec)]

        custom_gene_colors <- color_function(length(unique(custom_ordered_gene_vec)))
        custom_colors_result <- list("gene_colors" = custom_gene_colors,
                                     "gene_names" = custom_ordered_gene_vec)

        # Now if there are any eSNPs - get their gene colors:
        if (custom_fun_status()) {

            if (nrow(custom_loci_fun_snp_data()) > 0) {

                custom_fun_gene_vec <- unique(custom_loci_fun_snp_data()$gene_id)
                custom_ordered_fun_gene_vec <- custom_fun_gene_vec[order(custom_fun_gene_vec)]
                custom_colors_result[["fun_gene_colors"]] <-
                    custom_gene_colors[which(custom_ordered_gene_vec %in% custom_ordered_fun_gene_vec)]

            }

        }

        custom_colors_result

    })



# Create the output tables to display for the loci ----------------------

    output$custom_gene_table <- DT::renderDataTable(server = FALSE, {
        req(input$custom_genes)

        if (is.null(input$custom_snp_gene)) {
            custom_gene_table_display <- custom_gene_info_table() %>%
                dplyr::select(-gene_size)

        } else {
            custom_gene_table_display <- custom_loci_gene_data() %>%
                # Drop the gene chart tier and size columns:
                dplyr::select(-gene_chart_tier, -gene_size)
        }

        custom_gene_table_display %>%
            DT::datatable(escape = FALSE,
                          extensions = c("Buttons"),
                          options = list(
                              dom = 'Bfrtip',
                              buttons = list(
                                  list(extend = 'collection',
                                       buttons = list(
                                           list(extend = 'csv',
                                                filename = paste0(input$custom_loci, "_genes")),
                                           list(extend = 'excel',
                                                filename = paste0(input$custom_loci, "_genes"))),
                                       text = "Download gene table"
                                  )
                              )
                          ),
                          rownames = FALSE)
    })

    output$custom_snp_gene_table <- DT::renderDataTable(server = FALSE, {
        req(input$custom_snp_gene, input$custom_genes)

        # First decide how to construct the output table of custom_snp_gene_pairs
        # depending on the inclusion of functional SNPs. If there there was
        # no file uploaded to include functional SNPs than just use the required
        # positional SNPs:
        custom_snp_gene_pairs <- custom_snp_gene_table() %>%
            filter(loci_id == input$custom_loci) %>%
            mutate(is_positional = 1, is_functional = 0)

        # Now if there are functional SNPs loaded, join them but keep the
        # SNP-gene pairs distinct to allow for pairs that are both types
        if (!is.null(input$custom_fun_snp_gene)) {

            if (!is.null(custom_loci_fun_snp_data())) {

            # Only need to do this if there are any functional SNPs in the
            # selected loci:
            if (nrow(custom_loci_fun_snp_data()) > 0) {

                # Change the column names back:
                custom_fun_snp_gene_pairs <- custom_loci_fun_snp_data() %>%
                    dplyr::rename(bp = BP,
                                  snp_signal = `Functional SNP signal`)

                # Stack together with the positional pairs, but drop the indicators
                # for the type and just take the distinct rows:
                all_custom_snp_gene_pairs <- dplyr::select(custom_snp_gene_pairs,
                                                           -is_positional,
                                                           -is_functional) %>%
                    bind_rows(dplyr::select(custom_fun_snp_gene_pairs,
                                            -is_functional, -is_positional)) %>%
                    dplyr::distinct()

                # Now summarize the SNP-gene pair types to join back over:
                snp_gene_pair_types <- custom_snp_gene_pairs %>%
                    bind_rows(custom_fun_snp_gene_pairs) %>%
                    # Group by the snp_id, gene, and loci:
                    group_by(snp_id, gene_id, loci_id) %>%
                    summarize(is_positional = max(is_positional, na.rm = TRUE),
                              is_functional = max(is_functional, na.rm = TRUE)) %>%
                    ungroup()

                # Join these back over:
                custom_snp_gene_pairs <- all_custom_snp_gene_pairs %>%
                    dplyr::left_join(snp_gene_pair_types,
                                     by = c("snp_id", "gene_id", "loci_id"))

                }

            }
        }


        custom_snp_gene_pairs %>%
            DT::datatable(escape = FALSE,
                          extensions = c("Buttons", "KeyTable", "FixedColumns"),
                          options = list(
                              dom = 'Bfrtip',
                              scrollX = TRUE,
                              fixedColumns = list(leftColumns = 3),
                              keys = TRUE,
                              buttons = list(
                                  list(extend = 'collection',
                                       buttons = list(
                                           list(extend = 'csv',
                                                filename = paste0(input$custom_loci, "_snp_gene_pairs")),
                                           list(extend = 'excel',
                                                filename = paste0(input$custom_loci, "_snp_gene_pairs"))),
                                       text = "Download SNP-gene table"
                                  )
                              )
                          ),
                          rownames = FALSE)
    })



# Create the custom data plotly displays ----------------------------------

    # First determine the custom plot axes
    custom_plot_axes <- reactive({
        get_zoom_x_axis_limits(custom_loci_smooth_data()$BP,
                               custom_loci_gene_data()$start,
                               custom_loci_gene_data()$end)
    })


    output$custom_zoom_plot <- renderPlotly({

        # Start with the base-layer displaying rugs across the whole loci
        # regardless of the gene-level smoothing
        custom_base_snp_plot <- custom_loci_smooth_data() %>%
            ggplot(aes(x = BP)) +
            # Start with rugs of the positional SNPs
            geom_rug(data = filter(custom_loci_smooth_data(), is_fake_bp == 0),
                     alpha = 0.25, sides = "b") +
            labs(x = "BP", y = "Phenotype signal",
                 color = "Gene") +
            theme_bw()

        # Add the functional SNPs to the plot separately?
        if (custom_fun_status()) {
            # First check that there are eSNPs:
            if (nrow(custom_loci_fun_snp_data()) > 0) {

                custom_base_snp_plot <- custom_base_snp_plot +
                    geom_rug(data = {
                        custom_loci_fun_snp_data()
                    },
                    aes(color = gene_id),
                    alpha = 0.25, sides = "t") +
                    geom_segment(
                        data = {
                            custom_loci_fun_snp_data()
                        }, aes(xend = BP, y = 0, yend = `Functional SNP signal`,
                               color = gene_id),
                        size = 0.5, alpha = 0.25)
            }
        }

        # if (!is.null(input$custom_fun_snp_gene)) {
        #
        #     if (input$custom_display_fun_snps) {
        #
        #         # First check that there are eSNPs:
        #         if (nrow(custom_loci_fun_snp_data()) > 0) {
        #
        #             custom_base_snp_plot <- custom_base_snp_plot +
        #                 geom_rug(data = {
        #                     custom_loci_fun_snp_data()
        #                 },
        #                 aes(color = gene_id),
        #                 alpha = 0.25, sides = "t") +
        #                 geom_segment(
        #                     data = {
        #                         custom_loci_fun_snp_data()
        #                     }, aes(xend = BP, y = 0, yend = `Functional SNP signal`,
        #                            color = gene_id),
        #                     size = 0.5, alpha = 0.25)
        #         }
        #
        #     }
        #
        # }

        # Add reference results?
        if (custom_ref_status()) {
            # First check that there are reference signal data
            if (nrow(custom_loci_ref_data()) > 0) {

                custom_base_snp_plot <- custom_base_snp_plot +
                    geom_line(data = custom_loci_ref_data(),
                              aes(y = `Reference signal`), alpha = 0.8,
                              color = "gray", linetype = "dotted", size = 1)
            }
        }

        # if (!is.null(input$custom_ref_points)) {
        #     if (input$custom_display_reference) {
        #
        #         # First check that there are eSNPs:
        #         if (nrow(custom_loci_ref_data()) > 0) {
        #
        #             custom_base_snp_plot <- custom_base_snp_plot +
        #                 geom_line(data = custom_loci_ref_data(),
        #                           aes(y = `Reference signal`), alpha = 0.8,
        #                           color = "gray", linetype = "dotted", size = 1)
        #         }
        #     }
        # }


        # If LD loci level is one of the choices then display it:
        if ("ld_loci" %in% input$custom_smooth_level) {
            custom_base_snp_plot <- custom_base_snp_plot +
                geom_line(data = custom_loci_smooth_data(),
                          aes(y = `Smooth signal`), alpha = .8,
                          color = "black", size = 1)
        }

        # Add gene-level smoothing?
        if ("gene" %in% input$custom_smooth_level) {
            custom_base_snp_plot <- custom_base_snp_plot +
                geom_line(data = custom_gene_smooth_data(), alpha = 0.8,
                          aes(y = `Smooth signal`,
                              group = gene_id, color = gene_id))
        }

        # Finally add the color scale that should be used regardless of any settings
        custom_base_snp_plot <- custom_base_snp_plot +
            scale_color_manual(values = custom_gene_color_list()$gene_colors,
                               breaks = custom_gene_color_list()$gene_names,
                               drop = FALSE)

        # Next create the gene display below:
        custom_gene_plot <- custom_loci_gene_data() %>%
            ggplot(aes(xmin = start, xmax = end,
                       ymin = gene_chart_tier - 0.45,
                       ymax = gene_chart_tier + 0.45,
                       fill = gene_id, color = gene_id,
                       text = gene_id)) +
            geom_rect(alpha = 1) +
            scale_fill_manual(values = custom_gene_color_list()$gene_colors,
                              breaks = custom_gene_color_list()$gene_names,
                              drop = FALSE) +
            scale_color_manual(values = custom_gene_color_list()$gene_colors,
                               breaks = custom_gene_color_list()$gene_names,
                               drop = FALSE) +
            theme_minimal() +
            scale_y_reverse() +
            theme(legend.position = "bottom",
                  axis.text.y = element_blank()) +
            labs(x = "BP", fill = "gene_id", color = "gene_id")

        # First create the ggplotly object for the SNP plot, then modify its
        # legend layers before joining with the gene table in a subplot display:
        custom_base_snp_plotly <- ggplotly(
            custom_base_snp_plot +
                theme(legend.title = element_blank()) +
                scale_x_continuous(limits =
                                       c(custom_plot_axes()$chart_x_min - 1000,
                                         custom_plot_axes()$chart_x_max + 1000)),
            tooltip = c("BP", "Smooth signal", "gene_id", "Reference signal",
                        "Functional SNP signal")
        )

        # If there are any functional SNPs  displayed separately then update the legend
        # to handle these genes:
        if (!is.null(input$custom_fun_snp_gene)) {

            if (custom_fun_status()) {
                # First check that there are functional SNPs:
                if (nrow(custom_loci_fun_snp_data()) > 0) {
                    custom_base_snp_plotly <- custom_base_snp_plotly %>%
                        style(showlegend = FALSE,
                              traces = c(2:((length(custom_gene_color_list()$fun_gene_colors) * 2) + 1)))
                }
            }
        }


        # Next if gene-level then remove the final layer of colors from the legend:
        if ("gene" %in% input$custom_smooth_level) {

            # These are the final layer so how many layers are there currently?
            custom_n_base_snp_traces <- length(custom_base_snp_plotly$x$data)
            # How many genes are displayed:
            custom_n_gene_lines <- custom_gene_smooth_data() %>%
                pull(gene_id) %>%
                unique() %>%
                length()

            # Now update the trace
            custom_start_gene_trace <- custom_n_base_snp_traces -
                custom_n_gene_lines + 1
            custom_base_snp_plotly <- custom_base_snp_plotly %>%
                style(showlegend = FALSE,
                      traces = c(custom_start_gene_trace:custom_n_base_snp_traces))
        }

        # Make the plotly object:
        custom_gene_plotly <-
            ggplotly(custom_gene_plot +
                         theme(legend.title = element_blank(),
                               panel.grid.major.y = element_blank(),
                               panel.grid.minor.y = element_blank()) +
                         scale_x_continuous(limits = c(custom_plot_axes()$chart_x_min - 1000,
                                                       custom_plot_axes()$chart_x_max + 1000)),
                                    tooltip = "text")


        # Create the plotly subplot display
        custom_subplot_display <-
            subplot(custom_base_snp_plotly,
                    custom_gene_plotly,
                    shareX = TRUE, shareY = FALSE, titleY = TRUE, nrows = 2)


        custom_subplot_display %>%
            layout(legend = list(orientation = "v", x = 100, y = 0.5)) %>%
            config(
                toImageButtonOptions = list(
                    format = "svg",
                    filename = paste0(input$custom_loci, "_custom_zoom_plot"),
                    width = 1200,
                    height = 900
                )
            )


    })



})
