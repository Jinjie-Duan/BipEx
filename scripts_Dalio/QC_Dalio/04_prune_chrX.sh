#!/bin/bash
#$ -l h_vmem=2g
#$ -j y
#$ -l os=RedHat6

source /broad/software/scripts/useuse

reuse PLINK2

plink --bfile /stanley/genetics/analysis/bipolar_dalio/plink/filterGT.chrX \
      --indep 50 5 2\
      --out /stanley/genetics/analysis/bipolar_dalio/plink/pruned/04_chrX
