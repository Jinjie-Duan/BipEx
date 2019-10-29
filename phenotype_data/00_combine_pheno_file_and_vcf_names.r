library(data.table)

vcf_names <- fread("gsutil cat gs://raw_data_bipolar_dalio_w1_w2_hail_02/dalio_bipolar_w1_w2/Dalio_W1_W2_GRCh38_exomes_sample_IDs.tsv")
pheno_file <- fread("BIP_phenotype_information_cleaned_new_subtype_information_added_and_harmonised.tsv")

# Differences between the names in each of the files
setdiff(vcf_names$s, pheno_file$SAMPLE_ALIAS)
setdiff(pheno_file$SAMPLE_ALIAS, vcf_names$s)

# The vast majority of these differences are because the pheno file has names that contain '.', ',', or ' '. 
# Let's replace those

pheno_file$SAMPLE_ALIAS <- gsub("[\\,,\\.,\\ ]", "_", pheno_file$SAMPLE_ALIAS)

fwrite(pheno_file, file="BIP_phenotype_information_cleaned_new_subtype_information_added_and_harmonised_final.tsv", row.names=FALSE, sep='\t')