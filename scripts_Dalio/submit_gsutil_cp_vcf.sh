#!/bin/bash
#$ -l h_vmem=8g
#$ -j y
#$ -o /humgen/atgu1/fs03/dpalmer/logs/upload.out
#$ -w e
#$ -b n
#$ -l h_rt=48:00:00

source /broad/software/scripts/useuse
use .google-cloud-sdk

gsutil cp /seq/dax/BiPolar_CasesControls1_Exomes/Exome/v2/BiPolar_CasesControls1_Exomes.vcf.gz gs://raw_data_bipolar_dalio_w1_w2/bipolar_wes_dalio_W1_W2/BiPolar_CasesControls_Exomes.vcf.gz
