# Run read_in.py to obtain a vds and list of sample IDs from the joint called .vcf.
# Move the resultant sample names to local directory
gsutil cp gs://raw_data_bipolar_dalio_w1_w2/bipolar_WES_Dalio_W1_W2/Dalio_W1_W2_samples.tsv ~/Repositories/bipolar_WES_3/data/Dalio/Dalio_W1_W2_samples.tsv

# Copy the joint called .vcf to the local directory so we can create a clean phenotype file.
scp dpalmer@ag:/seq/dax/BiPolar_CasesControls1_Exomes/Exome/v2/BiPolar_CasesControls1_Exomes.calling_metadata.txt ~/Repositories/bipolar_WES_3/data/Dalio/BiPolar_CasesControls1_Exomes.calling_metadata.txt

# Run matching_to_vcf_samples_Dalio.r, to obtain 
# ~/Repositories/bipolar_WES_3/Dalio_phenotypes.tsv
# Move this data to the cloud.

gsutil cp ~/Repositories/bipolar_WES_3/Dalio_phenotypes.tsv gs://dalio_bipolar_w1_w2/data/samples/Dalio_phenotypes.tsv

# Move the vcf sample names (obtained using hail on the cloud initially) to the cluster as well.
ssh dpalmer@ag
gsutil cp gs://raw_data_bipolar_dalio_w1_w2/bipolar_WES_Dalio_W1_W2/Dalio_W1_W2_samples.tsv ~/Repositories/bipolar_WES_Dalio_W1_W2/Dalio_W1_W2_samples.tsv

# Run 00_get_bam_info.r on the broad cluster, and move the result to to the cloud.
gsutil cp ~/Repositories/bipolar_WES_Dalio_W1_W2/bam_metrics.tsv gs://dalio_bipolar_w1_w2/data/samples/bam_metrics.tsv

# Start the cluster
cluster start mycluster

cluster submit mycluster 00_create_phenotype_keytable.py
cluster submit mycluster 01_load_vcf_filterGT.py

# Similarly for the low complexity regions.
gsutil cp gs://raw_data_bipolar/inputs/low_complexity_regions.interval_list gs://raw_data_bipolar_dalio_w1_w2/inputs/low_complexity_regions.interval_list

cluster submit mycluster 02_prefilter_variants.py
cluster submit mycluster 03_initial_sample_qc.py

# Move the results from the cloud for plotting:
gsutil cp gs://dalio_bipolar_w1_w2_hail_02/data/samples/03_initial_sample_qc.tsv ~/Repositories/bipolar_WES_4_hail_02/samples_Dalio/03_initial_sample_qc.tsv
# Also download the .bam metrics. These are used as sanity check.
gsutil cp gs://dalio_bipolar_w1_w2/data/samples/bam_metrics.tsv ~/Repositories/bipolar_WES_4_hail_02/samples_Dalio/bam_metrics.tsv
# Plot the QC metrics locally using 03_initial_sample_qc_plot.r
# Choose thresholds and save these in 03_initial_sample_qc_filter.r and run.

# Move the resulting list of samples to the cloud.
gsutil cp ~/Repositories/bipolar_WES_4_hail_02/samples_Dalio/03_initial_qc.keep.sample_list gs://dalio_bipolar_w1_w2_hail_02/data/samples/03_initial_qc.keep.sample_list
# Move the interval list of high LD regions.
gsutil cp gs://dalio_bipolar/inputs/high_LD_regions_b37.interval_list gs://raw_data_bipolar_dalio_w1_w2/inputs/high_LD_regions_b37.interval_list
# Run 04_export_plink.py
cluster submit mycluster 04_export_plink.py

# Copy the plink files to the broad cluster.
ssh dpalmer@ag
gsutil cp gs://dalio_bipolar_w1_w2_hail_02/data/plink/filterGT* /stanley/genetics/analysis/bipolar_dalio/plink/

# Run the pruning algorithm, and merge the results together.
scp 04_prune_autosomes.sh dpalmer@ag:Repositories/bipolar_WES_Dalio_W1_W2/prune_autosomes.sh
scp 04_prune_chrX.sh dpalmer@ag:Repositories/bipolar_WES_Dalio_W1_W2/prune_chrX.sh

qsub ~/Repositories/bipolar_WES_Dalio_W1_W2/prune_autosomes.sh
qsub ~/Repositories/bipolar_WES_Dalio_W1_W2/prune_chrX.sh

# Combine the resulting pruned files into one file using 04_combine_pruned_autosomes.sh.
scp 04_combine_pruned_autosome.sh dpalmer@ag:Repositories/bipolar_WES_Dalio_W1_W2/combine_pruned_autosome.sh
# bash combine_pruned_autosome.sh

# Move the pruned files to the cloud.
gsutil cp /stanley/genetics/analysis/bipolar_dalio/plink/pruned/04_prune.keep.variant_list gs://dalio_bipolar_w1_w2_hail_02/data/variants/
gsutil cp /stanley/genetics/analysis/bipolar_dalio/plink/pruned/04_chrX.prune.in gs://dalio_bipolar_w1_w2_hail_02/data/variants/

# Run 05_impute_sex.py, 
cluster submit mycluster 05_impute_sex.py
cluster submit mycluster 06_ibd.py

# Move the output files to local folder.
gsutil cp gs://dalio_bipolar_w1_w2_hail_02/data/samples/05_imputesex.tsv ../../samples_Dalio/05_imputesex.tsv
gsutil cp gs://dalio_bipolar_w1_w2_hail_02/data/samples/05_ycalled.tsv ../../samples_Dalio/05_ycalled.tsv
gsutil cp gs://dalio_bipolar_w1_w2_hail_02/data/samples/06_ibd.tsv ../../samples_Dalio/06_ibd.tsv

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
cluster start --vep mycluster 
# Perform variant annotation on the cloud.
cluster submit mycluster 07_annotate_variants.py
# The following just chounts the number of variants in a bunch of the classes - does this make sense?
cluster submit mycluster 07_check_annotate_variants.py

# Count the number of ultra rare variants and perform PCA.
cluster submit mycluster 08_ultra_rare_counts.py
cluster submit mycluster 09_pca.py

# Make sure that the 1000 Genomes vds is present
# gs://hail-datasets-hail-data/1000_Genomes_autosomes.phase_3.GRCh37.mt is the full 1KG dataset.
gsutil cp -r gs://raw_data_bipolar/data/ALL.1KG.qc.hardcalls.vds gs://raw_data_bipolar_dalio_w1_w2/data/ALL.1KG.qc.hardcalls.vds
cluster submit mycluster 10_pca_1kg.py

# Move the resulting files from the cloud to local folder.
gsutil cp gs://dalio_bipolar_w1_w2_hail_02/data/samples/08_URVs.tsv ../../samples_Dalio/08_URVs.tsv
gsutil cp gs://dalio_bipolar_w1_w2_hail_02/data/samples/09_pca_scores.tsv ../../samples_Dalio/09_pca_scores.tsv
gsutil cp gs://dalio_bipolar_w1_w2_hail_02/data/samples/10_pca_scores_1kg.tsv ../../samples_Dalio/10_pca_scores_1kg.tsv

# Plot and perform filtering.
# 08_ultra_rare_counts_plot.r, 09_10_pca_plot.r.
# Move the resultant list of samples to be filtered out to the cloud.
gsutil cp ../../samples_Dalio/10_european.strict.sample_list gs://dalio_bipolar_w1_w2_hail_02/data/samples/10_european.strict.sample_list
gsutil cp ../../samples_Dalio/10_european.loose.sample_list gs://dalio_bipolar_w1_w2_hail_02/data/samples/10_european.loose.sample_list

# Run further PCA - this uses the strict definition of the Europeans.
# Potentially remove outliers here.
cluster submit mycluster 11_pca_EUR.py
# Move the resulting files from the cloud to local folder.
gsutil cp gs://dalio_bipolar_w1_w2_hail_02/data/samples/11_pca_scores.strict_european.tsv ../../samples_Dalio/11_pca_scores.strict_european.tsv
# Run 11_pca_EUR_plot.r to determine if there are any outliers.

cluster submit mycluster 12_pca_AJ.py
gsutil cp gs://dalio_bipolar_w1_w2_hail_02/data/samples/12_pca_scores_MGH_JH_samples.tsv ../../samples_Dalio/12_pca_scores_MGH_JH_samples.tsv

# Run 12_pca_us_AJ_plot.r to create a list of AJs and combine it with our Europeans, 
# after removing outliers from the output of 11_pca_EUR_plot.r and 12_pca_us_AJ_plot.r.
gsutil cp ../../samples_Dalio/12_european_and_AJ.sample_list gs://dalio_bipolar_w1_w2_hail_02/data/samples/12_european_and_AJ.sample_list
gsutil cp ../../samples_Dalio/12_mainland_european.sample_list gs://dalio_bipolar_w1_w2_hail_02/data/samples/12_mainland_european.sample_list
gsutil cp ../../samples_Dalio/12_swedes.sample_list gs://dalio_bipolar_w1_w2_hail_02/data/samples/12_swedes.sample_list
gsutil cp ../../samples_Dalio/12_european.sample_list gs://dalio_bipolar_w1_w2_hail_02/data/samples/12_european.sample_list

# Move the output to the cloud.
# There are different parts to run PCA on.
# 1) The mainland Europeans, together with EUR 1kg.
# 2) The Swedes, together with EUR 1kg.
# 3) Both the Swedes and the mainland Europeans, together with 1kg and AJs.
# 4) All the strictly defined 'Europeans' and 1kg.

# Output of these can then be used to determine where the Finns are.
cluster submit mycluster 13_pca_AJ_1kg.py
cluster submit mycluster 13_pca_mainland_EUR_1kg.py
cluster submit mycluster 13_pca_SWE_1kg.py
cluster submit mycluster 13_pca_EUR_1kg.py

# Move the results down from the cloud.
gsutil cp gs://dalio_bipolar_w1_w2_hail_02/data/samples/13_pca_scores.european_and_AJ.1kg.tsv ../../samples_Dalio/13_pca_scores.european_and_AJ.1kg.tsv
gsutil cp gs://dalio_bipolar_w1_w2_hail_02/data/samples/13_pca_scores.mainland_european.1kg.tsv ../../samples_Dalio/13_pca_scores.mainland_european.1kg.tsv
gsutil cp gs://dalio_bipolar_w1_w2_hail_02/data/samples/13_pca_scores.swedes.1kg.tsv ../../samples_Dalio/13_pca_scores.swedes.1kg.tsv
gsutil cp gs://dalio_bipolar_w1_w2_hail_02/data/samples/13_pca_scores.european.1kg.tsv ../../samples_Dalio/13_pca_scores.european.1kg.tsv

# Plot the PCA results - do things look clean?
# 13_pca_AJ_1kg_plot.r.
# Determine whether further samples should be removed.
gsutil cp ../../samples_Dalio/13_AJ_classify.sample_list gs://dalio_bipolar_w1_w2_hail_02/data/samples/13_AJ_classify.sample_list

cluster submit mycluster 15_final_variant_qc.py

# Move the results to local folder.
gsutil cp gs://dalio_bipolar_w1_w2_hail_02/data/variants/15_final_qc.variants.tsv.bgz ../../variants_Dalio/15_final_qc.variants.tsv.bgz

# Plot, filter, and move the subset of variants to the cloud.
# Plot using 15_final_variant_qc_plot.r, filter using 15_final_variant_qc_filter.r.
gsutil cp ../../variants_Dalio/15_final_qc.keep.variant_list gs://dalio_bipolar_w1_w2_hail_02/data/variants/15_final_qc.keep.variant_list

# Perform sample QC.
cluster submit mycluster 16_final_sample_qc.py

# Move the resulting files from the cloud to local folder.
gsutil cp gs://dalio_bipolar_w1_w2_hail_02/data/samples/16_final_qc.before.samples.tsv ../../samples_Dalio/16_final_qc.before.samples.tsv
gsutil cp gs://dalio_bipolar_w1_w2_hail_02/data/samples/16_final_qc.after.samples.tsv ../../samples_Dalio/16_final_qc.after.samples.tsv

# Create plots to visualise the QC'ed data - does it look good?
# Plot using 16_final_sample_qc_plot.r, filter using 16_filter_final_sample.r.
gsutil cp ../../samples_Dalio/16_final_qc.keep.sample_list gs://dalio_bipolar_w1_w2_hail_02/data/samples/16_final_qc.keep.sample_list
# gsutil cp ../../samples_Dalio/16_final_qc_remove_singleton_outliers.keep.sample_list gs://dalio_bipolar_w1_w2/data/samples/16_final_qc_remove_singleton_outliers.keep.sample_list

# Perform final creation of plink files and pruning.
cluster submit mycluster 17_export_plink_final.py

# Copy the plink files to the broad cluster.
ssh dpalmer@ag
gsutil cp gs://dalio_bipolar_w1_w2_hail_02/data/plink/final_qc* /psych/genetics_data/dpalmer/Dalio/plink/

# Run the pruning algorithm, and merge the results together.
scp 17_prune_autosomes_final.sh dpalmer@ag:Repositories/bipolar_WES_Dalio_W1_W2/final_prune_autosomes.sh
# scp 17_prune_chrX_final.sh dpalmer@ag:Repositories/bipolar_WES_Dalio_W1_W2/final_prune_chrX.sh

# Log on to the Broad cluster and submit.
qsub final_prune_autosomes.sh
# qsub final_prune_chrX.sh

# Combine the resulting pruned files into one file using 17_combine_pruned_final_autosomes.sh.

scp 17_combine_pruned_final_autosome.sh dpalmer@ag:Repositories/bipolar_WES_Dalio_W1_W2/final_combine_pruned_autosome.sh

# Move the resulting combined file to the cloud
gsutil cp 17_prune.final_qc.keep.variant_list gs://dalio_bipolar_w1_w2_hail_02/data/variants/17_prune.final_qc.keep.variant_list

# Create QC'ed .vds files.
cluster submit mycluster 18_create_qc_mt.py
cluster submit mycluster 18_check_qc_mt.py

gsutil cp gs://dalio_bipolar_w1_w2_hail_02/data/samples/18_final_qc.samples.tsv.bgz ../../samples_Dalio/18_final_qc.samples.tsv.bgz
gsutil cp gs://dalio_bipolar_w1_w2_hail_02/data/variants/18_final_qc.variants.tsv.bgz ../../samples_Dalio/18_final_qc.variants.tsv.bgz

# Create final PCA plots.
# 18_pca_final_plot.r


