# Overview of manuscript repository

This repository contains the code for reproducing the results and figures in _An agglomerative approach for LD loci testing with metadata and interactive localization_.

The folders are organized in the following manner:

- [code](https://github.com/ryurko/Agglomerative-LD-loci-testing/blob/master/code) - all scripts for initializing datasets, generating results and figures in manuscript

- [figures](https://github.com/ryurko/Agglomerative-LD-loci-testing/blob/master/figures) - files for final figures displayed in the manuscript (including supplementary materials)

- [data](https://github.com/ryurko/Agglomerative-LD-loci-testing/blob/master/data) - folder containing files necessary for the analysis (_NOTE: public GitHub repository version of folder only contains README file explaining the organization and datasets of final results_)

- [ld_loci_zoom_app](https://github.com/ryurko/Agglomerative-LD-loci-testing/blob/master/ld_loci_zoom_app) - LD loci zoom Shiny application folder, containing the necessary app files, scripts, and datasets

Due to their respective terms and conditions, we are unable to publicly post the initial datasets used in the analysis. However, all datasets considered in the manuscript are publicly available in the following locations:

- __Autism spectrum disorder (ASD)__ (Grove et al., 2019) and __Schizophrenia (SCZ)__ (Ruderfer et al., 2018) GWAS results were accessed from the [Psychiatric Genomics Consortium](https://www.med.unc.edu/pgc/download-results/).

- __Educational attainment (EA)__ (Lee et al., 2018) GWAS results were accessed from the [Social Science Genetic Association Consortium](https://www.thessgac.org/data).

- 1,000 Genomes reference data were accessed from the [MAGMA software page](https://ctg.cncr.nl/software/magma).

- __BrainVar__ eQTLs and WGCNA results were accessed from the [Supplemental Information section of the publication (Werling et al., 2020)](https://www.sciencedirect.com/science/article/pii/S2211124720303673#app2).

- [__Genotype-Tissue Expression (GTEx)__ V7 project datasets](https://gtexportal.org/home/datasets) (GTEx Consortium, 2015) were accessed for identifying eQTLs and generating gene co-expression modules. [The Tissue-Specific All SNP Gene Associations](https://gtexportal.org/home/datasets#filesetFilesDiv651) files for *Brain_Frontal_Cortex_BA9* and *Brain_Anterior_cingulate_cortex_BA24* were accessed for identifying GTEx eQTLs. The [RNA-Seq Data](https://gtexportal.org/home/datasets#filesetFilesDiv54) gene TPM counts were used for creating gene co-expression covariates.

- __Loss-of-function observed / expected upper fraction (LOEUF)__ values were accessed from the [gnomaAD browser (v2.1.1)](https://gnomad.broadinstitute.org/downloads)

## Contact

Ron Yurko: [ryurko@andrew.cmu.edu](mailto:ryurko@andrew.cmu.edu)

## References

- J Grove, et al., Identification of common genetic risk variants for autism spectrum disorder. _Nature genetics_ 51(3):431–444 (2019).

- DM Ruderfer, et al., Genomic dissection of bipolar disorder and schizophrenia, including 28 subphenotypes. _Cell_ 173(7):1705–1715 (2018).

- JJ Lee, et al., Gene discovery and polygenic prediction from a 1.1-million-person GWAS of educational attainment. _Nature genetics_ 50(8):1112 (2018).

- 1000 Genomes Project Consortium and others. An integrated map of genetic variation from 1,092 human genomes. _Nature_ 491:56 (2012).

- DM Werling, et al., Whole-genome and RNA sequencing reveal variation and transcriptomic coordination in the developing human prefrontal cortex. _Cell Reports_ 31(1):107489 (2020).

- GTEx Consortium, The genotype-tissue expression (gtex) pilot analysis: Multitissue gene regulation in humans. _Science_ 348:648–660 (2015).

- KJ Karczewski, et al., The mutational constraint spectrum quantified from variation in 141,456 humans. _Nature_ 581(7809):434-443 (2020).


