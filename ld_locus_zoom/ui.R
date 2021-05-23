# PURPOSE: Construct the UI for the LD loci zoom app using navbar pages to
#          separate between the display of our paper's ASD results, versus
#          the option for the user to upload their own results. There will also
#          be a separate page for the documentation of the app, providing info
#          with how to use the application.

library(shiny)
library(shinycssloaders)
library(shinycustomloader)
library(shinythemes)
library(shinyjs)
library(shinyalert)
library(markdown)
library(plotly)


# Define UI for application that draws a histogram
shinyUI(fluidPage(

    # Navbarpage structure for UI with cosmo theme
    navbarPage("LD locus zoom", theme = shinytheme("cosmo"),

               # Now display the tab for the ASD results (our paper's results)
               tabPanel("ASD results", fluid = TRUE, icon = icon("chart-area"),
                        # Sidebar layout with input and output definitions ----
                        sidebarLayout(

                            # Sidebar panel for inputs ----
                            sidebarPanel(
                                #titlePanel("Display options"),
                                # Select SNP-to-gene assignment type for display ----
                                selectizeInput(inputId = "assign_type",
                                               label = "Select SNP-to-gene assignment type:",
                                               choices = c("Positional", "Positional + eSNPs"),
                                               selected = "Positional + eSNPs"),

                                # Select CHR and loci ----
                                splitLayout(
                                    tags$head(tags$style(HTML(".shiny-split-layout > div {overflow: visible;}"))),
                                    uiOutput("chr_selection"),
                                    uiOutput("ld_loci_selection"),
                                    cellWidths = c("0%","50%", "50%")
                                ),

                                # Display smoothing at gene or LD loci level?
                                checkboxGroupInput("smooth_level",
                                                   label = "Level of smoothing displayed:",
                                                   choices = list("LD locus (black)" = "ld_loci",
                                                                  "Gene (colored by gene)" = "gene"),
                                                   selected = "ld_loci", inline = TRUE),

                                # Option for displaying eSNPs as peaks:
                                radioButtons(inputId = "display_esnps",
                                             label = "Display eSNPs?",
                                             choices = list("Yes" = TRUE, "No" = FALSE),
                                             selected = FALSE, inline = TRUE),

                                # Options for displaying SCZ and EA as background:
                                checkboxGroupInput(inputId = "background_pheno",
                                                   label = "Background phenotypes?",
                                                   choices = list("SCZ (dark red)" = "scz", "EA (light blue)" = "ea"),
                                                   selected = NULL, inline = TRUE)

                            ),

                            # Main panel for displaying outputs ----
                            mainPanel(
                                # Output: Visualization, Data, Description
                                tabsetPanel(type = "tabs",
                                            tabPanel("Visualization", withLoader(plotlyOutput(outputId = "zoomPlot"),
                                                                                 type = "html", loader = "dnaspin")),
                                            tabPanel("Genes", withSpinner(DT::dataTableOutput("gene_table"))),
                                            tabPanel("SNPs", withSpinner(DT::dataTableOutput("snp_gene_table"))))
                            )
                        ),
               ),

               tabPanel("Upload results", fluid = TRUE, icon = icon("file-upload"),
                        useShinyalert(),
                        sidebarLayout(

                            # Sidebar panel for inputs ----
                            sidebarPanel(

                                fileInput("custom_genes", "(Required) Upload table of gene-loci pairs",
                                          accept = c("text/csv",
                                                     "text/comma-separated-values,
                                                            text/plain", ".csv")),
                                fileInput("custom_snp_gene", "(Required) Upload table of SNP-gene pairs",
                                          accept = c("text/csv",
                                                     "text/comma-separated-values,
                                                            text/plain", ".csv")),
                                fileInput("custom_fun_snp_gene", "(Optional) Upload table of functional SNP-gene pairs",
                                          accept = c("text/csv",
                                                     "text/comma-separated-values,
                                                            text/plain", ".csv")),
                                fileInput("custom_ref_points", "(Optional) Upload reference signal data",
                                          accept = c("text/csv",
                                                     "text/comma-separated-values,
                                                            text/plain", ".csv")),

                                # Select CHR and loci ----
                                splitLayout(
                                    tags$head(tags$style(HTML(".shiny-split-layout > div {overflow: visible;}"))),
                                    uiOutput("custom_chr_selection"),
                                    uiOutput("custom_loci_selection"),
                                    cellWidths = c("0%","50%", "50%")
                                ),

                                # Display smoothing at gene or LD loci level?
                                checkboxGroupInput("custom_smooth_level",
                                                   label = "Level of smoothing displayed:",
                                                   choices = list("Locus (black)" = "ld_loci",
                                                                  "Gene (colored by gene)" = "gene"),
                                                   selected = "ld_loci", inline = TRUE),

                                # Option for displaying eSNPs as peaks:
                                #withSpinner(uiOutput("custom_add_fun_snps"), type = 0),
                                radioButtons(inputId = "custom_display_fun_snps",
                                             label = "Display functional SNPs?",
                                             choices = list("Yes" = TRUE, "No" = FALSE),
                                             selected = FALSE, inline = TRUE),


                                # Option for displaying null reference:
                                #withSpinner(uiOutput("custom_add_reference"), type = 0)
                                radioButtons(inputId = "custom_display_reference",
                                             label = "Display reference signal?",
                                             choices = list("Yes" = TRUE, "No" = FALSE),
                                             selected = FALSE, inline = TRUE)

                            ),

                            # Main panel for displaying outputs ----
                            mainPanel(
                                #h4(p("About the Project"))
                            #
                                # Output: Visualization, Data, Description
                                tabsetPanel(type = "tabs",
                                            tabPanel("Visualization",
                                                     withLoader(plotlyOutput(outputId = "custom_zoom_plot"),
                                                                type = "html", loader = "dnaspin")),
                                            tabPanel("Genes", withSpinner(DT::dataTableOutput("custom_gene_table"))),
                                            tabPanel("SNPs", withSpinner(DT::dataTableOutput("custom_snp_gene_table"))),
                                            tabPanel("Smoothing settings",
                                                     h4(p("You can change three settings for the displayed kernel smoothing interpolation:")),
                                                     numericInput("min_n_snps", "Minimum number of SNPs required for kernel smoothing:", 3, min = 3),
                                                     numericInput("loci_interp", "Minimum number of interpolation points for displayed locus-level smoothing:", 1000, min = 100),
                                                     numericInput("gene_interp", "Minimum number of interpolation points for displayed gene-level smoothing:", 100, min = 25)
                                                     ))
                            )
                        ),
               ),

               # Now display the tab for the ASD results (our paper's results)
               tabPanel("About", fluid = TRUE, icon = icon("info-circle"),
                        includeMarkdown("description.md"),
                        withMathJax())

               )
    )
)
