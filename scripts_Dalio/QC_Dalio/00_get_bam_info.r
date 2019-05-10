library(data.table)
# Take the .bam file information from the phenotype file.

# source('r_options_Dalio.r')
VCF_SAMPLE_IDS <- '~/Repositories/bipolar_WES_Dalio_W1_W2/Dalio_W1_W2_samples.tsv'
VCF_META <- '/seq/dax/BiPolar_CasesControls1_Exomes/Exome/v2/BiPolar_CasesControls1_Exomes.calling_metadata.txt'
BAM_METRICS <- '~/Repositories/bipolar_WES_Dalio_W1_W2/bam_metrics.tsv'

# Read in the .vcf sample names.
vcf_sample_ids <- read.table(VCF_SAMPLE_IDS, colClasses='character', header=1, sep='\t')[,1]

# Compare this to the collection of sample names that supposedly make it into the .vcf. 
vcf_meta <- fread(VCF_META, colClasses='character', header=1, sep='\t', data.table=FALSE)
vcf_meta <- vcf_meta[vcf_meta$included_in_final_vcf=='true',]

# Are they the same?
all(sort(vcf_meta$sample_name_in_vcf) == sort(vcf_sample_ids))
n_samples <- nrow(vcf_meta)

# Create a quality control data.frame.
df <- data.frame(Samples=vcf_meta$sample_name_in_vcf, batch=vcf_meta$project_or_cohort)

for(i in 1:n_samples) {
    # Read in the different files, and pull out the right column.
    tmp_table <- readLines(gsub('.bam', '.hybrid_selection_metrics', vcf_meta$bam[i]), 8)
    ix <- which(strsplit(tmp_table[7], split='\t')[[1]] == 'BAIT_SET')
    df$platform[i] <- strsplit(tmp_table[8], split='\t')[[1]][ix]

    df$pct_contamination[i] <- fread(gsub('.bam', '.selfSM', vcf_meta$bam[i]))$FREEMIX
    df$pct_chimeras[i] <- fread(gsub('.bam', '.alignment_summary_metrics', vcf_meta$bam[i]))$PCT_CHIMERAS[3]
    if (i %% 100 == 0) print(i)
    # Read in the different files, and pull out the right column.
    df$readgroups[i] <- as.integer(strsplit(system(paste0('wc -l ', gsub('/[^/]*.bam', '/', vcf_meta$bam[i]), 'analysis_files.txt'), intern=TRUE), split=' ')[[1]][1])-1
	if(df$readgroups[i]!=8) print(paste0("sample:", i, ", readgroups:", df$readgroups[i]))
}

fwrite(df, sep='\t', quote=FALSE, row.names=FALSE, col.names=names(df), file=BAM_METRICS)
