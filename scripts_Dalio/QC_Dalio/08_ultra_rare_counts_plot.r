library(hexbin)
library(gridExtra)
library(ggplot2)
library(ggExtra)
library(data.table)

# Load plotting functions:
source('r_functions_and_parameters/pretty_plotting.r')

# Get file locations and plotting locations.
source("r_functions_and_parameters/r_options_Dalio.r")

df <- fread(URV_FILE, sep='\t', header=TRUE, data.table=FALSE)
df <- df[sample(nrow(df), replace=FALSE),]

# Scatters of URV-SNPs against URV-Indels.
d <- ggplot(df, aes(x=n_URV_SNP, y=n_URV_indel, colour=phenotype.batch)) + 
geom_point(size=0.5, alpha=0.5) + 
scale_color_d3('category10') + theme_minimal()
d <- ggExtra::ggMarginal(d, type = "density",
  xparams = list(adjust=1), yparams=list(adjust=1.5))
print(d)

save_figures <- TRUE

y_labels <- c('Batch', 'Location')
y_label_batch <- c('', '')
titles <- c('Number of Singletons split by Batch and coloured by Location',
	'Number of Singletons split by Location and coloured by Phenotype')
titles <- c('', '')

create_pretty_boxplots(df, aes(x=phenotype.batch, y=n_URV_SNP), aes(colour=factor(phenotype.Location)), 
	y_label=y_label_batch[1], x_label='Number of Singletons', key_label='Location',
	title=titles[1], legend=TRUE, save_figure=save_figures,  file=paste0(PLOTS,'08_URVs_by_batch'))
create_pretty_boxplots(df, aes(x=phenotype.Location, y=n_URV_SNP), aes(colour=factor(phenotype.Phenotype)),
	y_label=y_label_batch[2], x_label='Number of Singletons', key_label='Phenotype',
	title=titles[2], legend=TRUE, save_figure=save_figures,  file=paste0(PLOTS,'08_URVs_by_location'))
