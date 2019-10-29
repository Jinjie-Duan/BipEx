00_merged_phenotypic_information.r and 00_merge_phenotypic_subtype_information.r check the phenotype information from BSP and read in the latest subphenotype information from the collaborators. We carefully go through and ensure things match correctly, and that not much has changed as a sanity check. 

This then generates BIP_phenotype_information_cleaned_new_subtype_information_added.tsv

From this file, we then noticed a few discrepancies between the coarse and fine phenotype definitions - these are manually curated in docs/qc.Rmd and written to BIP_phenotype_information_cleaned_new_subtype_information_added_and_harmonised.tsv.

Finally we checked the names of phenotypes in the vcf and compared them to those in this final phenotype file. 

To do this we just run 00_get_vcf_names.py to read in the initial matrix table and export the sample IDs

We then read in the sample IDs and compare them to the cleaned phenotype IDs to identify any differences and obvious name changes.

The resultant cleaned names are then written to BIP_phenotype_information_cleaned_new_subtype_information_added_and_harmonised_final.tsv in 00_combine_pheno_file_and_vcf_names.r