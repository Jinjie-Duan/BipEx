#!/bin/bash
#$ -l h_vmem=2g
#$ -j y
#$ -t 1-22

source /broad/software/scripts/useuse

use PLINK2

plink --bfile /psych/genetics_data/dpalmer/Dalio/plink/final_qc.chr${SGE_TASK_ID} \
      --indep 50 5 2 \
      --out /psych/genetics_data/dpalmer/Dalio/plink/pruned/17_final_qc.chr${SGE_TASK_ID}
      