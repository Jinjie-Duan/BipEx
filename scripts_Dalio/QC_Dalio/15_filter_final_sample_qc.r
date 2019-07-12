library(ggplot2)
library(ggsci)
library(dplyr)
library(data.table)

source("r_functions_and_parameters/r_options_Dalio.r")

df_after <- fread(SAMPLE_AFTER_QC_FILE, sep='\t', stringsAsFactors=FALSE, header=TRUE, data.table=FALSE) %>%
  mutate(phase='After Variant QC')

df_keep <- df_after['s']
print(paste0("Started with: ", nrow(df_keep), " samples"))

n_SDs <- 3

# r_ti_tv
df_keep_ti_tv <- group_by(df_after, phenotype.PROJECT_OR_COHORT) %>%
  summarise(mean=mean(sample_qc.r_ti_tv), sd=sd(sample_qc.r_ti_tv)) %>%
  inner_join(df_after, by='phenotype.PROJECT_OR_COHORT') %>%
  filter((sample_qc.r_ti_tv >= mean - n_SDs*sd & sample_qc.r_ti_tv <= mean + n_SDs*sd) | is.na(sd)) %>%
  select(s)

df_keep <- df_keep_ti_tv %>% inner_join(df_keep, by='s')
print(paste0("Remove Ti/Tv outliers: ", nrow(df_keep), " samples remain"))

# r_het_hom_var
df_keep_het_hom_var <- group_by(df_after, phenotype.PROJECT_OR_COHORT) %>%
  summarise(mean=mean(sample_qc.r_het_hom_var), sd=sd(sample_qc.r_het_hom_var)) %>%
  inner_join(df_after, by='phenotype.PROJECT_OR_COHORT') %>%
  filter((sample_qc.r_het_hom_var >= mean - n_SDs*sd & sample_qc.r_het_hom_var <= mean + n_SDs*sd) | is.na(sd)) %>%
  select(s)

df_keep <- df_keep_het_hom_var %>% inner_join(df_keep, by='s')
print(paste0("Remove Het/HomVar outliers: ", nrow(df_keep), " samples remain"))

# r_insertion_deletion
df_keep_insertion_deletion <- group_by(df_after, phenotype.PROJECT_OR_COHORT) %>%
  summarise(mean=mean(sample_qc.r_insertion_deletion), sd=sd(sample_qc.r_insertion_deletion)) %>%
  inner_join(df_after, by='phenotype.PROJECT_OR_COHORT') %>%
  filter((sample_qc.r_insertion_deletion >= mean - n_SDs*sd & sample_qc.r_insertion_deletion <= mean + n_SDs*sd) | is.na(sd)) %>%
  select(s)

df_keep <- df_keep_insertion_deletion %>% inner_join(df_keep, by='s')
print(paste0("Remove Ins/Del outliers: ", nrow(df_keep), " samples remain"))

# n_singletons
df_keep_n_singletons <- group_by(df_after, phenotype.LOCATION) %>%
  summarise(mean=mean(sample_qc.n_singleton), sd=sd(sample_qc.n_singleton)) %>%
  inner_join(df_after, by='phenotype.LOCATION') %>%
  filter((sample_qc.n_singleton >= mean - n_SDs*sd & sample_qc.n_singleton <= mean + n_SDs*sd) | is.na(sd)) %>%
  select(s)

df_keep <- df_keep_n_singletons %>% inner_join(df_keep, by='s')

print(paste0("Remove n_singletons outliers: ", nrow(df_keep), " samples remain"))

df_final_sample_summary <- data.table(Filter = c("Samples after population filters",
                          paste0("Within batch Ti/Tv ratio outside ", n_SDs, " standard deviations"),
                          paste0("Within batch Het/HomVar ratio outside ", n_SDs, " standard deviations"),
                          paste0("Within batch Insertion/Deletion ratio outside ", n_SDs, " standard deviations"),
                          paste0("Within location n singletons outside ", n_SDs, " standard deviations"),
                          "Samples after filters"),
                     Samples = c(nrow(df_after),
                                nrow(df_after) - nrow(df_keep_ti_tv),
                                nrow(df_after) - nrow(df_keep_het_hom_var),
                                nrow(df_after) - nrow(df_keep_insertion_deletion),
                                nrow(df_after) - nrow(df_keep_n_singletons),
                                nrow(df_keep)))

fwrite(df_final_sample_summary, file='../../samples_Dalio/15_sample_count.tsv', quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')

# write out
fwrite(df_keep, file=FINAL_SAMPLE_LIST, quote=FALSE, row.names=FALSE, col.names=FALSE)
