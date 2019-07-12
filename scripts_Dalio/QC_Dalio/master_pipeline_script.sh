# Run read_in.py to obtain a vds and list of sample IDs from the joint called .vcf.

# This time around I obtained a phenotype information file from DSP and BSP, merged them and cleaned them using 
# 00_merge_phenotypic_information.r and 00_clean_phenotypic_information.r.
# Move this data to the cloud.

gsutil cp ~/Repositories/BipEx/phenotype_data/BIP_phenotype_information_cleaned_new_subtype_information_added.tsv gs://dalio_bipolar_w1_w2_hail_02/data/samples/Dalio_phenotypes.tsv

# Start the cluster
hailctl dataproc start mycluster

hailctl dataproc submit mycluster 00_create_phenotype_keytable.py
hailctl dataproc submit mycluster 01_load_vcf_filterGT.py

# Ensure target intervals and low complexity regions for hg38 are on the cloud.
gsutil cp ~/Repositories/BipEx/variants_Dalio/LCR-hs38.bed gs://raw_data_bipolar_dalio_w1_w2_hail_02/inputs/LCR-hs38.bed
gsutil cp ~/Repositories/BipEx/variants_Dalio/ice_coding_v1_targets.interval_list gs://raw_data_bipolar_dalio_w1_w2_hail_02/inputs/ice_coding_v1_targets.interval_list
gsutil cp ~/Repositories/BipEx/variants_Dalio/ice_coding_v1_padded_targets.interval_list gs://raw_data_bipolar_dalio_w1_w2_hail_02/inputs/ice_coding_v1_padded_targets.interval_list

hailctl dataproc submit mycluster 02_prefilter_variants.py
hailctl dataproc submit mycluster 03_initial_sample_qc.py

# Plot the QC metrics locally using 03_initial_sample_qc_plot.r
# Choose thresholds and save these in 03_initial_sample_qc_filter.r and run.

# Move the resulting list of samples to the cloud.
gsutil cp ~/Repositories/BipEx/samples_Dalio/03_initial_qc.keep.sample_list gs://dalio_bipolar_w1_w2_hail_02/data/samples/03_initial_qc.keep.sample_list

# Move the interval list of high LD regions.
# First need to create it by lifting over an interval list or bed file.
# I lifted over https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)
# to build 38.

# Ended up using CrossMap, as hail resulted in a lot of NAs for intervals.

# The website is http://crossmap.sourceforge.net/
# I downloaded the chain file from sourceforge: GRCh37_to_GRCh38.chain.gz (found on the above website).
# insalled via pip3 and ran
# CrossMap.py bed GRCh37_to_GRCh38.chain.gz ../../variants_Dalio/b37_high_ld.bed ../../variants_Dalio/b38_high_ld.bed

gsutil cp ~/Repositories/BipEx/variants_Dalio/b38_high_ld.bed gs://raw_data_bipolar_dalio_w1_w2/inputs/b38_high_ld.bed

# Run 04_export_plink.py
hailctl dataproc submit mycluster 04_export_plink.py

# Copy the plink files to the broad cluster.
ssh dpalmer@ag
mkdir /stanley/genetics/analysis/bipolar_dalio/plink_b38
gsutil cp gs://dalio_bipolar_w1_w2_hail_02/data/plink/filterGT* /stanley/genetics/analysis/bipolar_dalio/plink_b38/

# Run the pruning algorithms and merge the results together.
# Just open an interactive node and use these scripts.

# 04_prune.sh

# This should move everything to the cloud...if not, run the following.
gsutil cp /stanley/genetics/analysis/bipolar_dalio/plink_b38/pruned/04_prune.keep.variant_list gs://dalio_bipolar_w1_w2_hail_02/data/variants/
gsutil cp /stanley/genetics/analysis/bipolar_dalio/plink_b38/pruned/04_chrX.prune.in gs://dalio_bipolar_w1_w2_hail_02/data/variants/

# Run 05_impute_sex.py, 
hailctl dataproc submit mycluster 05_impute_sex.py
hailctl dataproc submit mycluster 06_ibd.py

# Plot, and filter: Sex check and IBD check. Send the resultant list to the cloud.
# Plot using 05_impute_sex_plot.r and 06_ibd_plot.r.

# Run 06_ibd_filter.r 
gsutil cp ../../samples_Dalio/05_sexcheck.remove.sample_list gs://dalio_bipolar_w1_w2_hail_02/data/samples/05_sexcheck.remove.sample_list
gsutil cp ../../samples_Dalio/06_ibd.remove.sample_list gs://dalio_bipolar_w1_w2_hail_02/data/samples/06_ibd.remove.sample_list

# Move the nonpsych gnomad variant list across to the latest bucket.
gsutil cp gs://dalio_bipolar/data/variants/gnomad.r2.0.1.nonpsych.variant_list gs://dalio_bipolar_w1_w2_hail_02/data/variants/gnomad.r2.0.1.nonpsych.variant_list

# Move vep configuration file to the cluster.
gsutil cp vep_configuration.json gs://raw_data_bipolar/inputs/vep_configuration.json

# Start up a new cluster specifically for vep annotation.
hailctl dataproc start --vep mycluster 
# Perform variant annotation on the cloud.
hailctl dataproc submit mycluster 07_annotate_variants_vep_38.py
hailctl dataproc submit mycluster 07_liftovers_for_annotation.py
hailctl dataproc submit mycluster 07_annotate_variants.py
# The following just chounts the number of variants in a bunch of the classes - does this make sense?
hailctl dataproc submit mycluster 07_check_annotate_variants.py

# Count the number of ultra rare variants and perform PCA.
hailctl dataproc submit mycluster 08_ultra_rare_counts.py
hailctl dataproc submit mycluster 09_pca.py

# Make sure that the 1000 Genomes matrix table is present
# gs://hail-datasets-hail-data/1000_Genomes_autosomes.phase_3.GRCh38.mt is the full 1KG dataset. Create the hardcalls from this if required.
hailctl dataproc submit mycluster 10_pca_1kg.py

# Plot and perform filtering.
# 08_ultra_rare_counts_plot.r, 09_10_pca_plot.r.
# Move the resultant list of samples to be filtered out to the cloud.
gsutil cp ../../samples_Dalio/10_european.strict.sample_list gs://dalio_bipolar_w1_w2_hail_02/data/samples/10_european.strict.sample_list
gsutil cp ../../samples_Dalio/10_european.loose.sample_list gs://dalio_bipolar_w1_w2_hail_02/data/samples/10_european.loose.sample_list

# Run further PCA - this uses the strict definition of the Europeans.
# Potentially remove outliers here.
hailctl dataproc submit mycluster 11_pca_EUR.py
# Run 11_pca_EUR_plot.r to determine if there are any outliers.

hailctl dataproc submit mycluster 12_pca_aj.py

# Run 12_pca_us_aj_plot.r to create a list of AJs and combine it with our Europeans, 
# after removing outliers from the output of 11_pca_EUR_plot.r and 12_pca_us_aj_plot.r.
gsutil cp ../../samples_Dalio/12_european_and_aj.sample_list gs://dalio_bipolar_w1_w2_hail_02/data/samples/12_european_and_aj.sample_list
gsutil cp ../../samples_Dalio/12_mainland_european.sample_list gs://dalio_bipolar_w1_w2_hail_02/data/samples/12_mainland_european.sample_list
gsutil cp ../../samples_Dalio/12_swedes.sample_list gs://dalio_bipolar_w1_w2_hail_02/data/samples/12_swedes.sample_list
gsutil cp ../../samples_Dalio/12_european.sample_list gs://dalio_bipolar_w1_w2_hail_02/data/samples/12_european.sample_list

# Move the output to the cloud.
# There are different parts to run PCA on.
# 1) Both the Swedes and the mainland , together with 1kg and AJs.
# 2) All the strictly defined 'Europeans' and 1kg.

hailctl dataproc submit mycluster 13_pca_aj_1kg.py
hailctl dataproc submit mycluster 13_pca_EUR_1kg.py

# Plot the PCA results - do things look clean?
# 13_pca_aj_1kg_plot.r.
# Determine whether further samples should be removed.
gsutil cp ../../samples_Dalio/13_aj_classify.sample_list gs://dalio_bipolar_w1_w2_hail_02/data/samples/13_aj_classify.sample_list
gsutil cp ../../samples_Dalio/12_aj.sample_list gs://dalio_bipolar_w1_w2_hail_02/data/samples/12_aj.sample_list


hailctl dataproc submit mycluster 14_final_variant_qc.py

# Plot, filter, and move the subset of variants to the cloud.
# Plot using 15_final_variant_qc_plot.r, filter using 15_final_variant_qc_filter.r.
gsutil cp ../../variants_Dalio/14_final_qc.keep.variant_list gs://dalio_bipolar_w1_w2_hail_02/data/variants/14_final_qc.keep.variant_list

# Perform sample QC.
hailctl dataproc submit mycluster 15_final_sample_qc.py

# Create plots to visualise the QC'ed data - does it look good?
# Plot using 15_final_sample_qc_plot.r, filter using 15_filter_final_sample.r.
gsutil cp ../../samples_Dalio/15_final_qc.keep.sample_list gs://dalio_bipolar_w1_w2_hail_02/data/samples/15_final_qc.keep.sample_list

# Perform final creation of plink files and pruning.
hailctl dataproc submit mycluster 16_export_plink_final.py

# Copy the plink files to the broad cluster.
ssh dpalmer@ag
gsutil cp gs://dalio_bipolar_w1_w2_hail_02/data/plink/final_qc* /stanley/genetics/analysis/bipolar_dalio/plink_b38/

# 16_prune_final.sh

# This should move everything to the cloud...if not, run the following.
gsutil cp 16_prune.final_qc.keep.variant_list gs://dalio_bipolar_w1_w2_hail_02/data/variants/

# Create QC'ed .vds files.
hailctl dataproc submit mycluster 17_create_qc_mt.py
hailctl dataproc submit mycluster 17_check_qc_mt.py

# Create final PCA plots.
17_pca_final_plot.r


