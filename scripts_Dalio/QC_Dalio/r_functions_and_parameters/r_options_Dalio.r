# 00_get_bam_info.r

VCF_SAMPLE_IDS <- '~/Repositories/bipolar_WES_Dalio_W1_W2/Dalio_W1_W2_samples.tsv'
VCF_META <- '/seq/dax/BiPolar_CasesControls1_Exomes/Exome/v2/BiPolar_CasesControls1_Exomes.calling_metadata.txt'
BAM_METRICS <- '~/Repositories/bipolar_WES_Dalio_W1_W2/bam_metrics.tsv'

# 03_initial_sample_qc_filter.r
QC_FILE <- '../../samples_Dalio/03_initial_sample_qc.tsv'
SAMPLE_LIST_INITIAL_QC <- '../../samples_Dalio/03_initial_qc.keep.sample_list'

# 03_initial_sample_qc_plot.r
PLOTS <- '../../QC_plots/sample_plots/'
# Define some thresholds 
# T_sample_callRate <- 0.9
T_sample_callRate <- 0.94
# T_pct_contaminination <- 0.004
T_pct_contaminination <- 0.002
T_pct_chimeras <- 0.0015
T_dpMean <- 30
# T_gqMean <- 82.5
T_gqMean <- 85

# 05_impute_sex_plot.r
IMPUTESEX_FILE <- '../../samples_Dalio/05_imputesex.tsv'
SEXCHECK_LIST <- '../../samples_Dalio/05_sexcheck.remove.sample_list'
# PLOTS <- '../../QC_plots/sample_plots/'
Y_NCALLED_FILE <- '../../samples_Dalio/05_ycalled.tsv'

# 06_ibd_plot.r
IBD_FILE <- '../../samples_Dalio/06_ibd.tsv'
IBD_THRESHOLD <- 0.2

# 06_ibd_filtered.r
SAMPLE_LIST_IBD <- '../../samples_Dalio/06_ibd.remove.sample_list'

# 08_ultra_rare_counts_plot.r
URV_FILE <- '../../samples_Dalio/08_URVs.tsv'

# 09_10_pca_plot.r
PCA_SCORES <- '../../samples_Dalio/09_pca_scores.tsv'
PCA_1KG_SCORES <- '../../samples_Dalio/10_pca_scores_1kg.tsv'
EUROPEAN_SAMPLES_STRICT <- '../../samples_Dalio/10_european.strict.sample_list'
EUROPEAN_SAMPLES_LOOSE <- '../../samples_Dalio/10_european.loose.sample_list'
EUROPEAN_SAMPLES_EXCLUDING_URV_OUTLIERS <- '../../samples_Dalio/10_european.no_URV_outliers.sample_list'

T_nURVSNP <- 180
T_nURVIndel <- 17

# 12_pca_us_PCA_plot.r
PCA_SCORES_USA <- "../../samples_Dalio/12_pca_scores_MGH_JH_samples.tsv"

EUROPEAN_AND_AJ_SAMPLES <- "../../samples_Dalio/12_european_and_AJ.sample_list"
EUROPEANS <- "../../samples_Dalio/12_european.sample_list"
MAINLAND_EUROPEANS <- "../../samples_Dalio/12_mainland_european.sample_list"
SWEDES <- "../../samples_Dalio/12_swedes.sample_list"

# 13_pca_AJ_1kg_plots.r
PCA_EUR_1KG_AJ_SCORES <- '../../samples_Dalio/13_pca_scores.european_and_AJ.1kg.tsv'
PCA_MAINLAND_EUR_1KG_SCORES <- '../../samples_Dalio/13_pca_scores.mainland_european.1kg.tsv'
PCA_SWE_1KG_SCORES <- '../../samples_Dalio/13_pca_scores.swedes.1kg.tsv'
PCA_EUR_1KG_SCORES <- '../../samples_Dalio/13_pca_scores.european.1kg.tsv'

# 15_final_variant_qc_plot.r
VARIANT_QC_FILE <- '../../variants_Dalio/15_final_qc.variants.tsv.bgz'
T_variant_call_rate  <- 0.97
T_absdiff <- 0.02
T_pHWE <- 1e-6

# 15_final_variant_qc_filter.r
VARIANT_LIST <- '../../variants_Dalio/15_final_qc.keep.variant_list'

# 16_final_sample_qc_plot.r
SAMPLE_BEFORE_QC_FILE <- '../../samples_Dalio/16_final_qc.before.samples.tsv'
SAMPLE_AFTER_QC_FILE <- '../../samples_Dalio/16_final_qc.after.samples.tsv'

# 16_final_sample_qc_filter.r
FINAL_SAMPLE_LIST <- '../../samples_Dalio/16_final_qc.keep.sample_list'
FINAL_SAMPLE_LIST_REMOVE_SINGLETON_OUTLIERS <- '../../samples_Dalio/16_final_qc_remove_singleton_outliers.keep.sample_list'

# 18_pca_final_plot.r
FINAL_SAMPLE_QC <- '../../samples_Dalio/18_final_qc.samples.tsv'
FINAL_VARIANT_QC <- '../../variants_Dalio/18_final_qc.variants.tsv'

FINAL_SAMPLE_QC_NO_URV <- '../../samples_Dalio/18_final_qc_no_URV_sample_outliers.samples.tsv'
