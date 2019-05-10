library(ggplot2)
library(dplyr)
library(ggsci)
library(data.table)

source('r_functions_and_parameters/pretty_plotting.r')
source('r_functions_and_parameters/r_options_Dalio.r')

df <- fread(FINAL_SAMPLE_QC, sep='\t', stringsAsFactors=FALSE, header=TRUE, data.table=FALSE)
df <- df[sample(nrow(df)),]
save_figures <- TRUE

for (i in c(1,3,5)) {
  aes <- aes_string(x=paste0('pca.PC',i), y=paste0('pca.PC',i+1), color='Phenotype')
  create_pretty_scatter(df, aes,
  	save_figure=save_figures,
  	file=paste0(PLOTS,'18_PC',i,'_PC',i+1,'_final_PCs'),
  	x_label=paste0('Principal Component ', i),
  	y_label=paste0('Principal Component ', i+1))
  aes <- aes_string(x=paste0('pca.PC',i), y=paste0('pca.PC',i+1), color='Location')
  create_pretty_scatter(df, aes,
  	save_figure=save_figures,
  	file=paste0(PLOTS,'18_PC',i,'_PC',i+1,'_final_PCs_Location'),
  	x_label=paste0('Principal Component ', i),
  	y_label=paste0('Principal Component ', i+1))
}

# Take a look at the allele frequency in discovEHR.
df_var <- fread(FINAL_VARIANT_QC, sep='\t', stringsAsFactors=FALSE, header=TRUE, data.table=FALSE)

discovEHR_AF <- as.numeric(gsub("^.*\\:([^\\}]*)}", "\\1", df_var$DiscovEHR_AF))
discovEHR_AF <- pmin(1-discovEHR_AF, discovEHR_AF)
hist(discovEHR_AF)
