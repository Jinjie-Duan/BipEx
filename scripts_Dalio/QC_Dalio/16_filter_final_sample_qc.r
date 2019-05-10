library(ggplot2)
library(ggsci)
library(dplyr)
library(data.table)

source("r_functions_and_parameters/r_options_Dalio.r")

df_after <- fread(SAMPLE_AFTER_QC_FILE, sep='\t', stringsAsFactors=FALSE, header=TRUE, data.table=FALSE) %>%
  mutate(phase='After Variant QC')

df_keep <- df_after['s']
print(paste0("Started with: ", nrow(df_keep), " samples"))

# r_ti_tv
df_keep <- group_by(df_after, phenotype.batch) %>%
  summarise(mean=mean(sample_qc.r_ti_tv), sd=sd(sample_qc.r_ti_tv)) %>%
  inner_join(df_after, by='phenotype.batch') %>%
  filter((sample_qc.r_ti_tv >= mean - 3*sd & sample_qc.r_ti_tv <= mean + 3*sd) | is.na(sd)) %>%
  select(s) %>%
  inner_join(df_keep, by='s')
print(paste0("Remove Ti/Tv outliers: ", nrow(df_keep), " samples remain"))

# r_het_hom_var
df_keep <- group_by(df_after, phenotype.batch) %>%
  summarise(mean=mean(sample_qc.r_het_hom_var), sd=sd(sample_qc.r_het_hom_var)) %>%
  inner_join(df_after, by='phenotype.batch') %>%
  filter((sample_qc.r_het_hom_var >= mean - 3*sd & sample_qc.r_het_hom_var <= mean + 3*sd) | is.na(sd)) %>%
  select(s) %>%
  inner_join(df_keep, by='s')
print(paste0("Remove Het/HomVar outliers: ", nrow(df_keep), " samples remain"))

# r_insertion_deletion
df_keep <- group_by(df_after, phenotype.batch) %>%
  summarise(mean=mean(sample_qc.r_insertion_deletion), sd=sd(sample_qc.r_insertion_deletion)) %>%
  inner_join(df_after, by='phenotype.batch') %>%
  filter((sample_qc.r_insertion_deletion >= mean - 3*sd & sample_qc.r_insertion_deletion <= mean + 3*sd) | is.na(sd)) %>%
  select(s) %>%
  inner_join(df_keep, by='s')
print(paste0("Remove Ins/Del outliers: ", nrow(df_keep), " samples remain"))

# n_singletons
df_keep <- group_by(df_after, phenotype.Location) %>%
  summarise(mean=mean(sample_qc.n_singleton), sd=sd(sample_qc.n_singleton)) %>%
  inner_join(df_after, by='phenotype.Location') %>%
  filter((sample_qc.n_singleton >= mean - 3*sd & sample_qc.n_singleton <= mean + 3*sd) | is.na(sd)) %>%
  select(s) %>%
  inner_join(df_keep, by='s')
print(paste0("Remove n_singletons outliers: ", nrow(df_keep), " samples remain"))

# write out
fwrite(df_keep, file=FINAL_SAMPLE_LIST, quote=FALSE, row.names=FALSE, col.names=FALSE)
# fwrite(df_keep_remove_singleton_outliers, file=FINAL_SAMPLE_LIST_REMOVE_SINGLETON_OUTLIERS, quote=FALSE, row.names=FALSE, col.names=FALSE)
