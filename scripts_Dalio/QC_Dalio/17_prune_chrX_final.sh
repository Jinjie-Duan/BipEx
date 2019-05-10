#!/bin/bash
#$ -l h_vmem=2g
#$ -j y
#$ -o /humgen/atgu1/fs03/dpalmer/logs/17_prune.chrX.out

source /broad/software/scripts/useuse

reuse PLINK2

plink --bfile /psych/genetics_data/dpalmer/Dalio/plink/final_qc.chrX \
      --indep 50 5 2\
      --out /psych/genetics_data/dpalmer/Dalio/plink/pruned/17_final_qc.chrX
