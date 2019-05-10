library(dplyr)
source("r_functions_and_parameters/r_options_Dalio.r")

# Run the plotting again to ensure that the thresholds are as in the plots.
source("03_initial_sample_qc_plot.r")

df <- fread(QC_FILE, sep='\t', header=TRUE, data.table=FALSE)
names(df) <- gsub("phenotype\\.", "", names(df))
names(df) <- gsub("qc\\.", "", names(df))

df <- filter(df, Phenotype!='Unknown')

df_out <- filter(df, call_rate > T_sample_callRate) %>%
	filter(pct_contamination < T_pct_contaminination) %>%
	filter(pct_chimeras < T_pct_chimeras) %>%
	filter(dp_stats.mean > T_dpMean) %>%
	filter(gq_stats.mean > T_gqMean)

df_out <- df_out %>% select(s)
print(dim(df_out))

fwrite(df_out, file=SAMPLE_LIST_INITIAL_QC, quote=FALSE, row.names=FALSE, col.names=FALSE)
