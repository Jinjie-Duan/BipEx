#!/bin/bash
#$ -l h_vmem=8g
#$ -j y
#$ -cwd
#$ -o /humgen/atgu1/fs03/dpalmer/logs/get_bam_info.out
#$ -w e
#$ -b n
#$ -l h_rt=48:00:00

source /broad/software/scripts/useuse
use R-3.3

Rscript 00_get_bam_info.r
