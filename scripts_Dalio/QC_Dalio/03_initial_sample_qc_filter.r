library(dplyr)
library(data.table)
source("r_functions_and_parameters/r_options_Dalio.r")

# Run the plotting again to ensure that the thresholds are as in the plots.
source("03_initial_sample_qc_plot.r")

QC_FILE <- "gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/03_initial_sample_qc.tsv"
df <- fread(QC_FILE, sep='\t', header=TRUE, data.table=FALSE)
names(df) <- gsub("phenotype\\.", "", names(df))
names(df) <- gsub("qc\\.", "", names(df))

# Also want to count the number of samples not present in the vcf file.
df_pheno <- fread("../../phenotype_data/BIP_phenotype_information_cleaned_new_subtype_information_added.tsv")
df_pheno <- df_pheno[-which(is.na(df_pheno$PI)),]

# Remove samples with unknown phenotype, those from 1kg that were included, and those with high
# contamination or coverage as defined by the spreadsheet from Laura.
df <- df[-which(is.na(df$PI)),]

df_initial_summary_count <- data.table(
	Filter = c("Initial samples", "Unknown phenotype", "Low coverage or high contamination"),
	Samples = c(nrow(df),
				nrow(filter(df, PHENOTYPE_COARSE=='Unknown')),
				nrow(df[which(df$s %in% low_coverage),])))

df <- filter(df, PHENOTYPE_COARSE!='Unknown')
df <- df[-which(df$s %in% low_coverage),]

df_initial_summary_count <- rbindlist(list(df_initial_summary_count, list("Samples after initial filter", nrow(df))), use.names=FALSE)
fwrite(df_initial_summary_count, file='../../samples_Dalio/03_initial_sample_count.tsv', quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')

names(df) <- gsub("qc_padded_ice\\.", "", names(df))

df_out <- filter(df, call_rate > T_sample_callRate) %>%
	filter(PCT_CONTAMINATION < T_pct_contamination) %>%
	filter(PCT_CHIMERAS < T_pct_chimeras) %>%
	filter(dp_stats.mean > T_dpMean) %>%
	filter(gq_stats.mean > T_gqMean)

df_out <- df_out %>% select(s)
print(dim(df_out))

fwrite(df_out, file=SAMPLE_LIST_INITIAL_QC, quote=FALSE, row.names=FALSE, col.names=FALSE)

# Create the table too
df_summary_count <- data.table(
	Filter = c("Initial samples",
			   paste0("Sample call rate < ", T_sample_callRate),
			   paste0("% FREEMIX contamination > ", T_pct_contamination),
			   paste0("% chimeric reads > ", T_pct_chimeras),
			   paste0("Mean DP < ", T_dpMean),
			   paste0("Mean GQ < ", T_gqMean),
			   "Samples after sample QC filters"),
	Samples = c(nrow(df),
			    nrow(filter(df, call_rate <= T_sample_callRate)),
				nrow(filter(df, PCT_CONTAMINATION >= T_pct_contamination)),
				nrow(filter(df, PCT_CHIMERAS >= T_pct_chimeras)),
				nrow(filter(df, dp_stats.mean <= T_dpMean)),
				nrow(filter(df, gq_stats.mean <= T_gqMean)),
				nrow(df_out)),
	"Bipolar Cases" = c(nrow(df %>% filter(PHENOTYPE_COARSE == "Bipolar Disorder")),
			    nrow(filter(df, call_rate <= T_sample_callRate) %>% filter(PHENOTYPE_COARSE == "Bipolar Disorder")),
				nrow(filter(df, PCT_CONTAMINATION >= T_pct_contamination) %>% filter(PHENOTYPE_COARSE == "Bipolar Disorder")),
				nrow(filter(df, PCT_CHIMERAS >= T_pct_chimeras) %>% filter(PHENOTYPE_COARSE == "Bipolar Disorder")),
				nrow(filter(df, dp_stats.mean <= T_dpMean) %>% filter(PHENOTYPE_COARSE == "Bipolar Disorder")),
				nrow(filter(df, gq_stats.mean <= T_gqMean) %>% filter(PHENOTYPE_COARSE == "Bipolar Disorder")),
				nrow(df_out %>% filter(PHENOTYPE_COARSE == "Bipolar Disorder"))),
	Controls = c(nrow(df %>% filter(PHENOTYPE_COARSE == "Bipolar Disorder")),
			    nrow(filter(df, call_rate <= T_sample_callRate) %>% filter(PHENOTYPE_COARSE == "Control")),
				nrow(filter(df, PCT_CONTAMINATION >= T_pct_contamination) %>% filter(PHENOTYPE_COARSE == "Control")),
				nrow(filter(df, PCT_CHIMERAS >= T_pct_chimeras) %>% filter(PHENOTYPE_COARSE == "Control")),
				nrow(filter(df, dp_stats.mean <= T_dpMean) %>% filter(PHENOTYPE_COARSE == "Control")),
				nrow(filter(df, gq_stats.mean <= T_gqMean) %>% filter(PHENOTYPE_COARSE == "Control")),
				nrow(df_out %>% filter(PHENOTYPE_COARSE == "Control"))))

fwrite(df_summary_count, file='../../samples_Dalio/03_sample_count.tsv', quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')
