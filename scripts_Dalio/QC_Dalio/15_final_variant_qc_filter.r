library(data.table)
library(dplyr)

source("r_functions_and_parameters/r_options_Dalio.r")

df <- fread(cmd = paste("gzcat", VARIANT_QC_FILE),  data.table=FALSE, header=TRUE, sep='\t')
print(paste0("Initial number of variants: ", nrow(df)))
current_variants <- nrow(df)

df_out <- df %>% filter((qc.AF > 0) & (qc.AF < 1))
print(paste0("After removing invariant sites in the cleaned dataset: ", nrow(df_out)))
print(nrow(df_out) - current_variants)
current_variants <- nrow(df_out)

df_out <- df_out %>% filter(qc.call_rate >= T_variant_call_rate)
print(paste0("After ensuring that the overall call rate is good: ", nrow(df_out)))
print(nrow(df_out) - current_variants)
current_variants <- nrow(df_out)

df_out <- df_out %>% filter(control_qc.qc.call_rate >= T_variant_call_rate)
print(paste0("After ensuring that the control call rate is good: ", nrow(df_out)))
print(nrow(df_out) - current_variants)
current_variants <- nrow(df_out)

df_out <- df_out %>% filter(case_qc.qc.call_rate >= T_variant_call_rate)
print(paste0("After ensuring that the case call rate is good: ", nrow(df_out)))
print(nrow(df_out) - current_variants)
current_variants <- nrow(df_out)

df_out <- df_out %>% filter(abs(case_qc.qc.call_rate - control_qc.qc.call_rate) <= T_absdiff)
print(paste0("After ensuring that the difference in call rates is low: ", nrow(df_out)))
print(nrow(df_out) - current_variants)
current_variants <- nrow(df_out)

df_out <- df_out %>% filter(qc.p_value_hwe > T_pHWE)
print(paste0("After ensuring that sites pass HWE p-value threshold: ", nrow(df_out)))
print(nrow(df_out) - current_variants)
current_variants <- nrow(df_out)

df <- fread(cmd = paste("gzcat", VARIANT_QC_FILE),  data.table=FALSE, header=TRUE, sep='\t')
print("After removing invariant sites in the cleaned dataset: ")
print(nrow(df) - nrow(df %>% filter((qc.AF > 0) & (qc.AF < 1))))

print("After ensuring that the overall call rate is good: ")
print(nrow(df) - nrow(df %>% filter(qc.call_rate >= T_variant_call_rate)))

print("After ensuring that the control call rate is good: ")
print(nrow(df) - nrow(df %>% filter(control_qc.qc.call_rate >= T_variant_call_rate)))

print("After ensuring that the case call rate is good: ")
print(nrow(df) - nrow(df %>% filter(case_qc.qc.call_rate >= T_variant_call_rate)))

print("After ensuring that the difference in call rates is low: ")
print(nrow(df) - nrow(df %>% filter(abs(case_qc.qc.call_rate - control_qc.qc.call_rate) <= T_absdiff)))

print("After ensuring that sites pass HWE p-value threshold: " )
print(nrow(df) - nrow(df %>% filter(qc.p_value_hwe > T_pHWE)))

fwrite(df_out %>% select("locus", "alleles"), file=VARIANT_LIST, quote=FALSE, row.names=FALSE, col.names=TRUE, sep='\t')
