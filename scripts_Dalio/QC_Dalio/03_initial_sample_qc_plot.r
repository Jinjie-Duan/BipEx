library(data.table)
# Plotting functions:
source('r_functions_and_parameters/pretty_plotting.r')
# Thresholds and plotting file locations defined in r_options_Dalio.r
source("r_functions_and_parameters/r_options_Dalio.r")

# Watch out for O'Donovan!
df <- fread(QC_FILE, stringsAsFactors=FALSE, sep='\t', header=TRUE, data.table=FALSE)

# Let's rename the PI column so it's not stupid.
# Remove the 1000G samples.
names(df) <- gsub("phenotype\\.", "", names(df))
names(df) <- gsub("qc\\.", "", names(df))

df <- df[-which(is.na(df$PI)),]
df <- df[sample(nrow(df), replace=FALSE),]

fwrite(df, '../../samples_Dalio/03_initial_sample_qc_cleaned_check.tsv', sep='\t')

save_figures <- TRUE
# PDFs, no splitting.
create_pretty_hist(df, aes(x=call_rate), 'Call Rate', T_sample_callRate,
    title='Call Rate', save_figure=save_figures, file=paste0(PLOTS,'03_callRate_hist'))
create_pretty_hist(df, aes(x=pct_contamination), 'Contamination', T_pct_contaminination,
    binwidth=0.0001, xlim=c(0,0.01), title='% Contamination', save_figure=save_figures, file=paste0(PLOTS,'03_contamination_hist'))
create_pretty_hist(df, aes(x=pct_chimeras), 'Chimeric Reads', T_pct_chimeras,
    binwidth=0.00002, xlim=c(0,0.01), title='% Chimeric Reads', save_figure=save_figures, file=paste0(PLOTS,'03_chimeras_hist'))
create_pretty_hist(df, aes(x=dp_stats.mean), 'Mean Depth', T_dpMean,
    binwidth=2, xlim=c(10, 150), title='Mean Depth', save_figure=save_figures, file=paste0(PLOTS,'03_dpMean_hist'))
create_pretty_hist(df, aes(x=gq_stats.mean), 'Mean Genotype Quality', T_gqMean,
    binwidth=0.5, xlim=c(70, 100), title='Mean Genotype Quality', save_figure=save_figures, file=paste0(PLOTS,'03_gqMean_hist'))

# CDFs, no splitting.
create_pretty_cumulative(df, aes(call_rate), 'Call Rate', T_sample_callRate,
    xlim=c(0.75,1), title='Call Rate', save_figure=save_figures, file=paste0(PLOTS,'03_callRate_cdf'))
create_pretty_cumulative(df, aes(pct_contamination), 'Contamination', T_pct_contaminination,
    xlim=c(0,0.01), title='% Contamination', save_figure=save_figures, file=paste0(PLOTS,'03_contamination_cdf'))
create_pretty_cumulative(df, aes(pct_chimeras), 'Chimeric Reads', T_pct_chimeras,
    xlim=c(0,0.004), title='% Chimeric Reads', save_figure=save_figures, file=paste0(PLOTS,'03_chimeras_cdf'))
create_pretty_cumulative(df, aes(dp_stats.mean), 'Mean Depth', T_dpMean,
    xlim=c(10,150), title='Mean Depth', save_figure=save_figures, file=paste0(PLOTS,'03_dpMean_cdf'))
create_pretty_cumulative(df, aes(gq_stats.mean), 'Mean Genotype Quality', T_gqMean,
    xlim=c(70,100), title='Mean Genotype Quality', save_figure=save_figures, file=paste0(PLOTS,'03_gqMean_cdf'))

# Split by Location.
# DEV: need to be careful here - location is the location of the centre, not the location of where the sampling was done.
# For the presentation, make things larger and clearer.
legend_batch <- FALSE
legend_collection <- TRUE
legend_phenotype <- TRUE
save_figures <- TRUE
# y_label='Collection'
y_label_batch <- 'Batch'
y_label_batch <- ''
titles <- c('Call Rate',
    '% Contamination',
    '% Chimeric Reads',
    'Mean Depth (DP)',
    'Mean Genotype Quality (GQ)')
titles <- c('', '', '', '', '')

create_pretty_boxplots(df, aes(x=Location, y=call_rate), aes(color=batch),
    T_sample_callRate, x_label='Call Rate', y_label=y_label_batch, key_label='Batch',
    xlim=c(0.75,1), legend=legend_batch, title=titles[1], save_figure=save_figures,
    file=paste0(PLOTS,'03_callRate_by_collection'), n_ticks=5)
create_pretty_boxplots(df, aes(x=Location, y=pct_contamination), aes(color=batch),
    T_pct_contaminination, x_label='% Contamination', y_label=y_label_batch, key_label='Batch',
    xlim=c(0,0.01), legend=legend_batch, title=titles[2], save_figure=save_figures,
    file=paste0(PLOTS,'03_contaminiation_by_collection'), n_ticks=5)
create_pretty_boxplots(df, aes(x=Location, y=pct_chimeras), aes(color=batch),
    T_pct_chimeras, x_label='% Chimeric Reads', y_label=y_label_batch, key_label='Batch',
    xlim=c(0,0.01), legend=legend_batch, title=titles[3], save_figure=save_figures,
    file=paste0(PLOTS,'03_chimeras_by_collection'), n_ticks=5)
create_pretty_boxplots(df, aes(x=Location, y=dp_stats.mean), aes(color=batch),
    T_dpMean, x_label='Mean Depth', y_label=y_label_batch, key_label='Batch',
    xlim=c(10,150), legend=legend_batch, title=titles[4], save_figure=save_figures,
    file=paste0(PLOTS,'03_dpMean_by_collection'))
create_pretty_boxplots(df, aes(x=Location, y=gq_stats.mean), aes(color=batch),
    T_gqMean, x_label='Mean Genotype Quality', y_label=y_label_batch, key_label='Batch',
    xlim=c(70,100), legend=legend_batch, title=titles[5], save_figure=save_figures,
    file=paste0(PLOTS,'03_gqMean_by_collection'))

create_pretty_boxplots(df, aes(x=Location, y=call_rate), aes(color=Phenotype),
    T_sample_callRate, x_label='Call Rate', y_label='Collection', key_label='Phenotype',
    xlim=c(0.75,1), legend=legend_phenotype, title=titles[1], save_figure=save_figures,
    file=paste0(PLOTS,'03_callRate_by_collection_col_CC'))
create_pretty_boxplots(df, aes(x=Location, y=pct_contamination), aes(color=Phenotype),
    T_pct_contaminination, x_label='% Contamination', y_label='Collection', key_label='Phenotype',
    xlim=c(0,0.01), legend=legend_phenotype, title=titles[2], save_figure=save_figures,
    file=paste0(PLOTS,'03_contaminiation_by_collection_col_CC'))
create_pretty_boxplots(df, aes(x=Location, y=pct_chimeras), aes(color=Phenotype),
    T_pct_chimeras, x_label='% Chimeric Reads', y_label='Collection', key_label='Phenotype',
    xlim=c(0,0.01), legend=legend_phenotype, title=titles[3], save_figure=save_figures,
    file=paste0(PLOTS,'03_chimeras_by_collection_col_CC'))
create_pretty_boxplots(df, aes(x=Location, y=dp_stats.mean), aes(color=Phenotype),
    T_dpMean, x_label='Mean Depth', y_label='Collection', key_label='Phenotype',
    xlim=c(10,150), legend=legend_phenotype, title=titles[4], save_figure=save_figures,
    file=paste0(PLOTS,'03_dpMean_by_collection_col_CC'))
create_pretty_boxplots(df, aes(x=Location, y=gq_stats.mean), aes(color=Phenotype),
    T_gqMean, x_label='Mean Genotype Quality', y_label='Collection', key_label='Batch',
    xlim=c(70,100), legend=legend_phenotype, title=titles[5], save_figure=save_figures,
    file=paste0(PLOTS,'03_gqMean_by_collection_col_CC'))

# Split by case/control
df <- df[-which(df$Phenotype == "Unknown"),]
create_pretty_boxplots(df, aes(x=Phenotype, y=call_rate), aes(color=batch),
    T_sample_callRate, x_label='Call Rate', y_label='Status',
    key_label='Batch', xlim=c(0.75,1), legend=legend_batch, title=titles[1],
    save_figure=save_figures, file=paste0(PLOTS,'03_callRate_by_status'))
create_pretty_boxplots(df, aes(x=Phenotype, y=pct_contamination), aes(color=batch),
    T_pct_contaminination, x_label='% Contamination', y_label='Status',
    key_label='Batch', xlim=c(0,0.01), legend=legend_batch, title=titles[2],
    save_figure=save_figures, file=paste0(PLOTS,'03_contamination_by_status'))
create_pretty_boxplots(df, aes(x=Phenotype, y=pct_chimeras), aes(color=batch),
    T_pct_chimeras, x_label='% Chimeric Reads', y_label='Status',
    key_label='Batch', xlim=c(0,0.01), legend=legend_batch, title=titles[3],
    save_figure=save_figures, file=paste0(PLOTS,'03_chimeras_by_status'))
create_pretty_boxplots(df, aes(x=Phenotype, y=dp_stats.mean), aes(color=batch),
    T_dpMean, x_label='Mean Depth', y_label='Status',
    key_label='Batch', xlim=c(10,150), legend=legend_batch, title=titles[4],
    save_figure=save_figures, file=paste0(PLOTS,'03_dpMean_by_status'))
create_pretty_boxplots(df, aes(x=Phenotype, y=gq_stats.mean), aes(color=batch),
    T_gqMean, x_label='Mean Genotype Quality', y_label='Status',
    key_label='Batch', xlim=c(70,100), legend=legend_batch, title=titles[5],
    save_figure=save_figures, file=paste0(PLOTS,'03_gqMean_by_status'))

# Stratify by batch.
create_pretty_boxplots(df, aes(x=batch, y=call_rate), aes(color=Location),
    T_sample_callRate, y_label='Batch', x_label='Call Rate', key_label='Collection',
    xlim=c(0.75,1), legend=legend_collection, title=titles[1],
    save_figure=save_figures, file=paste0(PLOTS,'03_callRate_by_batch'))
create_pretty_boxplots(df, aes(x=batch, y=pct_contamination), aes(color=Location),
    T_pct_contaminination, y_label='Batch', x_label='% Contamination', key_label='Collection',
    xlim=c(0,0.01), legend=legend_collection, 
    title=titles[2],
    save_figure=save_figures, file=paste0(PLOTS,'03_contamination_by_batch'))
create_pretty_boxplots(df, aes(x=batch, y=pct_chimeras), aes(color=Location),
    T_pct_chimeras, y_label='Batch', x_label='% Chimeric Reads', key_label='Collection',
    xlim=c(0,0.01), legend=legend_collection, title=titles[3],
    save_figure=save_figures, file=paste0(PLOTS,'03_chimeras_by_batch'))
create_pretty_boxplots(df, aes(x=batch, y=dp_stats.mean), aes(color=Location),
    T_dpMean, y_label='Batch', x_label='Mean Depth', key_label='Collection',
    xlim=c(10,150), legend=legend_collection, title=titles[4],
    save_figure=save_figures, file=paste0(PLOTS,'03_dpMean_by_batch'))
create_pretty_boxplots(df, aes(x=batch, y=gq_stats.mean), aes(color=Location),
    T_gqMean, y_label='Batch', x_label='Mean Genotype Quality', key_label='Collection',
    xlim=c(70,100), legend=legend_collection, title=titles[5],
    save_figure=save_figures, file=paste0(PLOTS,'03_gqMean_by_batch'))

# Check that the bimodal region is due to resequencing...
df$breaks <- cut(df$readgroups, breaks=c(0,7,10,12,20))
create_pretty_boxplots(df, aes(x=batch, y=gq_stats.mean), aes(color=breaks),
    T_gqMean, y_label='gqMean', x_label='Mean Genotype Quality', key_label='Resequencing',
    legend=TRUE, jitter_size=1, save_figure=save_figures, file=paste0(PLOTS, '03_reseq_gq_batch'))
create_pretty_boxplots(df, aes(x=batch, y=dp_stats.mean), aes(color=breaks),
    T_dpMean, y_label='dpMean', x_label='Call Rate', key_label='Resequencing',
    legend=TRUE, jitter_size=1, save_figure=save_figures, file=paste0(PLOTS, '03_reseq_dp_batch'))
create_pretty_boxplots(df, aes(x=batch, y=call_rate), aes(color=breaks),
    T_sample_callRate, y_label='Batch', x_label='Call Rate', xlim=c(0.9,1), key_label='Resequencing',
    legend=TRUE, jitter_size=1, save_figure=save_figures, file=paste0(PLOTS, '03_reseq_callRate_batch'))
