#!/bin/bash
#$ -l h_vmem=2g
#$ -j y
#$ -t 1-22
#$ -l os=RedHat6

source /broad/software/scripts/useuse
mkdir /stanley/genetics/analysis/bipolar_dalio/plink/pruned

# This works with RedHat6 but this is now being phased out.

reuse PLINK2

plink --bfile /stanley/genetics/analysis/bipolar_dalio/plink/filterGT.chr${SGE_TASK_ID} \
      --indep 50 5 2\
      --out /stanley/genetics/analysis/bipolar_dalio/plink/pruned/04_chr${SGE_TASK_ID}
