rm(list=ls())
library(ggplot2)
library(ggsci)
library(dplyr)
library(data.table)

source("r_functions_and_parameters/pretty_plotting.r")
source("r_functions_and_parameters/r_options_Dalio.r")

save_figure <- TRUE

df_before <- fread(SAMPLE_BEFORE_QC_FILE, sep='\t', stringsAsFactors=FALSE, header=TRUE, data.table=FALSE) %>%
  mutate(phase='Before Variant QC')
df_after <- fread(SAMPLE_AFTER_QC_FILE, sep='\t', stringsAsFactors=FALSE, header=TRUE, data.table=FALSE) %>%
  mutate(phase='After Variant QC')
df <- bind_rows(df_before, df_after) %>%
  mutate(phase=factor(phase, levels=c('Before Variant QC', 'After Variant QC')))

## Split by Location.
# Colour by Case-status.
# Number of singletons.

y_labels <- 'Location'
y_labels <- ''

create_pretty_boxplots(df, aes(y=sample_qc.n_singleton, x=factor(phenotype.Location)),
    aes(color=phenotype.Phenotype), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='Number of Singletons',
    save_figure=save_figure, file=paste0(PLOTS, '16_nSingletonsbyLocationColPheno'), y_label=y_labels)

# rHetHomVar
create_pretty_boxplots(df, aes(y=sample_qc.r_het_hom_var, x=factor(phenotype.Location)),
    aes(color=phenotype.Phenotype), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rHetHomVar',
    save_figure=save_figure, file=paste0(PLOTS, '16_rHetHomVarbyLocationColPheno'), y_label=y_labels)

# rInsertionDeletion
create_pretty_boxplots(df, aes(y=sample_qc.r_insertion_deletion, x=factor(phenotype.Location)),
    aes(color=phenotype.Phenotype), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rInsertionDeletion',
    save_figure=save_figure, file=paste0(PLOTS, '16_rInsertionDeletionbyLocationColPheno'), y_label=y_labels)

# rTiTv
create_pretty_boxplots(df, aes(y=sample_qc.r_ti_tv, x=factor(phenotype.Location)),
    aes(color=phenotype.Phenotype), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rTiTv',
    save_figure=save_figure, file=paste0(PLOTS, '16_rTiTvbyLocationColPheno'), n_ticks=5, y_label=y_labels)

## Colour by Batch.
# Number of singletons.
create_pretty_boxplots(df, aes(y=sample_qc.n_singleton, x=factor(phenotype.Location)),
    aes(color=phenotype.batch), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='Number of Singletons',
    save_figure=save_figure, file=paste0(PLOTS, '16_nSingletonsbyLocationColBatch', y_label=y_labels))

# rHetHomVar
create_pretty_boxplots(df, aes(y=sample_qc.r_het_hom_var, x=factor(phenotype.Location)),
    aes(color=phenotype.batch), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rHetHomVar',
    save_figure=save_figure, file=paste0(PLOTS, '16_rHetHomVarbyLocationColBatch'), y_label=y_labels)

# rInsertionDeletion
create_pretty_boxplots(df, aes(y=sample_qc.r_insertion_deletion, x=factor(phenotype.Location)),
    aes(color=phenotype.batch), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rInsertionDeletion',
    save_figure=save_figure, file=paste0(PLOTS, '16_rInsertionDeletionbyLocationColBatch'), y_label=y_labels)

# rTiTv
create_pretty_boxplots(df, aes(y=sample_qc.r_ti_tv, x=factor(phenotype.Location)),
    aes(color=phenotype.batch), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rTiTv',
    save_figure=save_figure, file=paste0(PLOTS, '16_rTiTvbyLocationColBatch'), n_ticks=5, y_label=y_labels)

### Split by batch.
y_labels <- 'Batch'
y_labels <- ''
## Colour by Case-status.
# Number of singletons.
create_pretty_boxplots(df, aes(y=sample_qc.n_singleton, x=phenotype.batch),
    aes(color=phenotype.Phenotype), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='Number of Singletons',
    legend=TRUE, save_figure=save_figure, file=paste0(PLOTS, '16_nSingletonsbyBatchColPheno'), y_label=y_labels)

# rHetHomVar
create_pretty_boxplots(df, aes(y=sample_qc.r_het_hom_var, x=phenotype.batch),
    aes(color=phenotype.Phenotype), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rHetHomVar',
    legend=TRUE, save_figure=save_figure, file=paste0(PLOTS, '16_rHetHomVarbyBatchColPheno'), y_label=y_labels)

# rInsertionDeletion
create_pretty_boxplots(df, aes(y=sample_qc.r_insertion_deletion, x=phenotype.batch),
    aes(color=phenotype.Phenotype), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rInsertionDeletion',
    legend=TRUE, save_figure=save_figure, file=paste0(PLOTS, '16_rInsertionDeletionbyBatchColPheno'), y_label=y_labels)

# rTiTv
create_pretty_boxplots(df, aes(y=sample_qc.r_ti_tv, x=phenotype.batch),
    aes(color=phenotype.Phenotype), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rTiTv',
    legend=TRUE, save_figure=save_figure, file=paste0(PLOTS, '16_rTiTvbyBatchColPheno'), n_ticks=5, y_label=y_labels)

## Colour by Location.
# Number of singletons
create_pretty_boxplots(df, aes(y=sample_qc.n_singleton, x=phenotype.batch),
    aes(color=phenotype.Location), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='Number of Singletons',
    legend=TRUE, save_figure=save_figure, file=paste0(PLOTS, '16_nSingletonsbyBatchColLocation'), y_label=y_labels)

# rHetHomVar
create_pretty_boxplots(df, aes(y=sample_qc.r_het_hom_var, x=phenotype.batch),
    aes(color=phenotype.Location), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rHetHomVar',
    legend=TRUE, save_figure=save_figure, file=paste0(PLOTS, '16_rHetHomVarbyBatchColLocation'), y_label=y_labels)

# rInsertionDeletion
create_pretty_boxplots(df, aes(y=sample_qc.r_insertion_deletion, x=phenotype.batch),
    aes(color=phenotype.Location), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rInsertionDeletion',
    legend=TRUE, save_figure=save_figure, file=paste0(PLOTS, '16_rInsertionDeletionbyBatchColLocation'), y_label=y_labels)
# rTiTv
create_pretty_boxplots(df, aes(y=sample_qc.r_ti_tv, x=phenotype.batch),
    aes(color=phenotype.Location), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rTiTv',
    legend=TRUE, save_figure=save_figure, file=paste0(PLOTS, '16_rTiTvbyBatchColLocation'), n_ticks=5, y_label=y_labels)

q <- function(x) {
  print(subset(x, x < (mean(x) - 3*sd(x))| x > (mean(x) + 3*sd(x))))
  y <- subset(x, x < (mean(x) - 3*sd(x))| x > (mean(x) + 3*sd(x)))
  if(length(y)==0) y <- NA
  return(y)
}

# Number of singletons
p <- create_pretty_boxplots(df, aes(y=sample_qc.n_singleton, x=phenotype.Location),
    aes(color=phenotype.Location), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rHetHomVar',
    legend=TRUE) + stat_summary(fun.y=q, geom="point", position=position_dodge(1), color='grey')
print(p)

# rHetHomVar
p <- create_pretty_boxplots(df, aes(y=sample_qc.r_het_hom_var, x=phenotype.batch),
    aes(color=phenotype.Location), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rHetHomVar',
    legend=TRUE) + stat_summary(fun.y=q, geom="point", position=position_dodge(1), color='grey')
print(p)

# rInsertionDeletion
p <- create_pretty_boxplots(df, aes(y=sample_qc.r_insertion_deletion, x=phenotype.batch),
    aes(color=phenotype.Location), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rInsertionDeletion',
    legend=TRUE) + stat_summary(fun.y=q, geom="point", position=position_dodge(1), color='grey')
print(p)

# rTiTv
p <- create_pretty_boxplots(df, aes(y=sample_qc.r_ti_tv, x=phenotype.batch),
    aes(color=phenotype.Location), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rTiTv',
    legend=TRUE) + stat_summary(fun.y=q, geom="point", position=position_dodge(1), color='grey')
print(p)

