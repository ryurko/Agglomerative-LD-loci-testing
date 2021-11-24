# Organization of code

This folder contains the scripts for all steps of the analysis with each folder prefix, `sXX`, indicating the step number `XX`, e.g., `s01_prep_input_data` contains the scripts for the first step of the analysis in the manuscript. Additionally, all code for generating the manuscript and supplement figures is located in the `create_figures` folder.

All code is presented uder the assumption of working in a R Project version of this repository.

## Order of scripts

The scripts within the folders for the first seven steps should be run in the following order:

1. [s01_prep_input_data](https://github.com/ryurko/Agglomerative-LD-loci-testing/blob/master/code/s01_prep_input_data)

+ `init_gencode_v21_table.R`

+ `init_asd_scz_ea_gwas.R`

2. [s02_find_eqtl_pairs](https://github.com/ryurko/Agglomerative-LD-loci-testing/blob/master/code/s02_find_eqtl_pairs)

+ `get_brainvar_eqtls.R` 

+ `get_gtex_brain_eqtls.R`

3. [s03_assign_snps_to_genes](https://github.com/ryurko/Agglomerative-LD-loci-testing/blob/master/code/s03_assign_snps_to_genes)

+ `assign_positional_and_esnps.R`

+ `init_server_version.R`

4. [s04_create_ld_induced_genes_loci](https://github.com/ryurko/Agglomerative-LD-loci-testing/blob/master/code/s04_create_ld_induced_genes_loci)

+ Helper functions required in same directory: `agglom_algo_functions.R`

+ `update_tidy_snp_gene_data.R`

+ `init_ld_reference_cor_by_chr.R`

+ `merge_positional_esnps_genes.R`

+ `merge_positional_genes.R`

+ `create_tidy_loci_tables.R`

+ `examine_algo_counts.R`

5. [s05_create_gene_locus_level_data](https://github.com/ryurko/Agglomerative-LD-loci-testing/blob/master/code/s05_create_gene_locus_level_data)

+ Helper functions required in same directory using `Rcpp`: `sim_quad_test_stat_fast.cpp`

+ `get_positional_esnps_locus_mc_pvalues.R`

+ `get_positional_locus_mc_pvalues.R`

+ `construct_adapt_model_metadata.R`

6. [s06_perform_adapt](https://github.com/ryurko/Agglomerative-LD-loci-testing/blob/master/code/s06_perform_adapt)

+ `tune_positional_esnps.R`

+ `tune_positional.R`

+  Used to decide settings for AdaPT gradient boosted trees: `examine_tune_results.R`

+ `get_adapt_int_only.R`

+ `get_adapt_cv_positional_esnps.R`

+ `get_adapt_cv_positional.R`

+ `save_snp_gene_result_tables.R`

+ `save_gene_lists.R`

7. [s07_prep_app_smoothing](https://github.com/ryurko/Agglomerative-LD-loci-testing/blob/master/code/s07_prep_app_smoothing)

+ Helper functions required in same directory: `kernel_smoothing_fns.R`

+ `init_ABC_snp_gene_kernel_input.R` (where `ABC` is either `asd`, `scz`, or `ea`)

+ `get_ABC_locus_gcv_bandwidth.R`

+ `generate_ABC_kernel_smoothing_results.R`

Note that `s08_examine_results` contains scripts for exploring results in various ways.

