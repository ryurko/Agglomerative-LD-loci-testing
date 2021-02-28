# About the LD loci zoom application

The [Shiny](https://shiny.rstudio.com/) LD loci zoom application is divided between 
a display of the results for autism spectrum disorder (ASD) in the _ASD results_ tab, 
and the option to upload data to display custom results using the _Upload results_ tab. 
This documentation first provides an overview of the display in the _ASD results_ tab 
followed by instructions for uploading custom results. 

Use the following links to jump ahead to different sections:

1. [Visualization overview for ASD results](#viz-overview)

2. [Export gene and SNP tables](#gene-snps)

3. [Directions for uploading results](#upload)

4. [Contact information](#contact)


## Visualization overview for ASD results <a name="viz-overview"></a>

### LD loci selection 

There are two sets of ASD results to choose from based on the SNP-to-gene assignment type for initializing gene $g$'s set of SNPs $S$:

1. _Positional_: $S_g = S_g^{pos}$ is the set of SNPs with genomic positions $BP$ located within gene $g$'s start and end positions (without padding).

2. _Positional + eSNPs_: $S_g = S_g^{pos} \cup S_g^{eQTL}$, where $S_g^{eQTL}$ are the set of eSNPs associated with gene $g$.

The _Positional + eSNPs_ results are the default selection. Then you choose to examine a LD loci associated with ASD by first picking the chromosome (CHR), and then an LD loci ID located in the selected CHR. (_Note: LD loci IDs have no meaning beyond ensuring unique identifiers for the LD locis within each CHR_.) 

### Kernel smoothing localization

Given the set of SNPs $S$ for a selected LD loci $m^*$, we smooth over the squared $z$ statistics of the LD loci's positional SNPs $S^{Pos}_{m^*} \subseteq S_{m^*}$ using the Nadaraya–Watson estimator:

$$
\hat{z}_{m^*}^2(BP) = \sum_{i \in S_{m^*}^{Pos}} z_i^2 \frac{K_h(BP_i, BP)}{\sum_{j \in S_{m^*}^{Pos}} K_h(BP_j, BP)},
$$

where $K_h$ is a one-dimensional Gaussian kernel with bandwidth $h$ selected separately for each LD loci using generalized cross-validation (as implemented in the [`np` package](https://github.com/JeffreyRacine/R-Package-np)). The default option is to display the smooth signal on the y-axis with interpolation between the positional SNPs at the LD loci level (denoted by the black line). The LD loci's SNPs are represented by the rugs along the x-axis. For reference, the gray dotted line denotes a point-wise $95^{th}$ percentile of 1,000 simulations for squared null Gaussian random variables, $z_{m^*}^2$ where $z_{m^*} \sim \text{Normal}(\boldsymbol{0}_{m^*}, \boldsymbol{\Sigma}_{m^*})$, given the LD structure $\boldsymbol{\Sigma}_{m^*}$ of the selected LD loci $m^*$. 

Additionally, all of the genes in the selected LD loci are displayed below the x-axis with rectangles denoting their start and end positions, with different colors for each gene. (_Note: gene arrangement has no meaning beyond ensuring no overlap between each gene's start and end positions._) You can choose to display the gene-level smooth signal (lines colored by gene) which uses the Nadaraya-Watson estimator above, but apply separately for each gene.

### Display eSNPs and background phenotypes

You can display any subset of eSNPs $S^{eQTL} \subseteq S$ separately as bars with their heights indicating individual SNP-level signal (colors matching their associated genes). This separation is due to the presence of intergenic eSNPs, however any eSNPs that are also positionally assigned to genes, $S^{Pos} \cap S^{eQTL}$, are included in the positional smoothing above.

For additional reference, you can choose to display the kernel smoothing results for schizophrenia (SCZ, dark red) and educational attainment (EA, light blue) in the background. Both are normalized to appear on the same signal scale as ASD to simplify comparisons.

### Plotly features

The visualization is built with the [plotly](https://plotly.com/r/) library enabling interactive features such as hover text describing the signal strength or indicating the gene. You can click on the gene names in the legend to the right for highlighting subsets of genes, which will also select the subset of eSNPs or gene-level smoothing lines. Additionally, the plotly features enable you to zoom in and out, pan the figure, as well as download the plot (the camera icon) as an SVG file.

## Export gene and SNP tables <a name="gene-snps"></a>

You can view and export tables of the genes and SNP-gene pairs in the selected LD loci via the __Genes__ and __SNPs__ tabs respectively. Both tabs include download buttons above the tables enabling the option to export all of the genes or SNP-gene pairs in either csv or Excel format. __All genomic positions in the visualization and tables (gene's start and end, SNPs' BP) are based on genome assembly version GRCh38__.

The genes table includes columns for Ensembl ID, gene name, [GENCODE v21](https://www.gencodegenes.org/human/release_21.html) start and end positions, strand  (+ or -) indicating direction, as well as a link to the gene's [GWAS catalog](https://www.ebi.ac.uk/gwas/home) page. _Note: the GWAS catalog link is constructed using the gene's name and does not guarantee that the gene is in the GWAS catalog._

The SNP-gene table includes a row for every SNP-gene pair in the selected LD loci, meaning that if a SNP is assigned to multiple genes then it will have a row for each gene. The SNP-gene table includes columns for:

* SNP's BP, alleles, rs ID, and GWAS summary statistics for ASD, SCZ, and EA
* Ensembl ID for the gene it is assigned to and how: `is_positional = 1` if it is positionally assigned, `is_esnp = 1` if it is an associated eSNP (SNPs can be both be positional and eSNPs for a gene)
* [GWAS catalog](https://www.ebi.ac.uk/gwas/home) link. _Note: the GWAS catalog link is constructed using the SNP's rs ID and does not guarantee that the SNP is in the GWAS catalog._

## Directions for uploading results <a name="upload"></a>

From the _Upload results_ tab, you can choose to upload custom results for displaying 
localized signal in a similar manner to visualization in the _ASD results_ tab. There
are two required datasets (__in CSV format__) for you to upload in order for the visualization to be
displayed:

+ __Gene-loci table__ - each row corresponds to unique gene-loci pair with five required columns (_and value types_):

  1. `gene_id` - unique gene ID, e.g., Ensembl ID (_character or numeric_)
  
  2. `start` - gene's start position (_numeric_)
  
  3. `end` - gene's end position (_numeric_)
  
  4. `chr` - gene's CHR (_character or numeric_)
  
  5. `loci_id` - unique loci ID, e.g., LD loci IDs in _ASD results_ tab, which each gene is assigned to that must be unique within a CHR (_character or numeric_)


  
+ __SNP-gene table__ - each row corresponds to unique SNP-gene pair with four required columns (_and value types_):

  1. `snp_id` - unique SNP ID, e.g., rs ID or CHR:BP (_character or numeric_)
  
  2. `gene_id` - unique gene ID designating which gene the SNP is assigned to and corresponds to the same set of gene IDs in the gene-loci table above (_character or numeric_)
  
  3. `bp` - SNP's genomic position (_numeric_)
  
  4. `snp_signal` - SNP's signal to apply smoothing over, e.g., $-\log_{10}$($p$-value) or squared $z$ statistics (_numeric_)
  
  
Using these two provided tables, the loci IDs from the gene-loci table are joined to the SNP-gene table and then the kernel smoothing localization from above is applied to user' `snp_signal` column. As before for the ASD results, the bandwidth is selected separately for each loci using generalized cross-validation. _Note: This process may take a few minutes to generate the smoothing results for all provided locis depending on the size of the uploaded datasets, but once complete you can switch between the displays for selected locis without delay._ As before with the _ASD results_ tab, the localized signal will be displayed in the _Visualization_ tab, with options to export tables of genes and SNPs using the _Genes_ and _SNPs_ tabs respectively. _Note: Any additional columns provided in the gene-loci and SNP-gene tables will be included in the exported tables_.

### Display functional SNPs and reference signals

Additionally, you can optionally upload two additional datasets (__in CSV format__) for display:

+ __Functional SNP-gene table__ - each row corresponds to a unique functional SNP-gene pair with four required columns (_and value types_):

  1. `snp_id` - unique SNP ID, e.g., rs ID or CHR:BP (_character or numeric_)
  
  2. `gene_id` - unique gene ID designating which gene the SNP is functionally assigned to and corresponds to the same set of gene IDs in the gene-loci table above (_character or numeric_)
  
  3. `bp` - SNP's genomic position (_numeric_)
  
  4. `snp_signal` - SNP's signal to display, e.g., $-\log_{10}$($p$-value) or squared $z$ statistics, that is assumed to be on the same scale as the SNP-gene table's `snp_signal` (_numeric_)


+ __Reference signal data__ - each row corresponds to a unique genomic position in a loci with three required columns (_and value types_):

  1. `loci_id` - unique loci ID corresponding to the same set of loci IDs in the gene-loci table above (_character or numeric_)
  
  2. `bp` - unique genomic position with a loci (_numeric_)
  
  3. `reference_signal` - reference signal value to display in comparison to the kernel smoothing signals (_numeric_)


The functional SNP-gene pairs can then be optionally displayed in the same manner as the associated eSNPs in the _ASD results_ tab, with vertical bars denoting the functional SNP's signal (colors matching their associated genes). The reference signal can then be optionally displayed as dashed, grey reference line in the background of the smoothing signal display, similarly to the pointwise null reference line in the _ASD results_ tab. (See here for example code demonstrating how to generate pointwise null values: __TODO INSERT GITHUB LINK__).

### Change kernel smoothing interpolation 

Finally, you can change three particular settings for the displayed kernel smoothing interpolation in the _Smoothing settings_ tab:

1. Minimum number of SNPs required for kernel smoothing in a loci or gene (default and minimum value is 3)

2. Minimum number of interpolation points for displayed loci-level smoothing (default is 1000, minimum value is 100)

3. Minimum number of interpolation points for displayed loci-level smoothing (default is 100, minimum value is 25)

The default options correspond to the settings used in the _ASD results_ tab, and if any of the three are changed then the kernel smoothing results must be updated for all locis.


## Contact <a name="contact"></a>

Email questions and comments to Ron Yurko at <ryurko@andrew.cmu.edu>.

