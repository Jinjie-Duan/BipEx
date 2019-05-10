#!/bin/bash

cd /stanley/genetics/analysis/bipolar_dalio/plink/pruned
cat 04_chr1.prune.in 04_chr2.prune.in > 04_prune.keep.variant_list_tmp

for i in `seq 3 22`;
do
	cat 04_prune.keep.variant_list_tmp 04_chr${i}.prune.in > 04_prune.keep.variant_list
	mv 04_prune.keep.variant_list 04_prune.keep.variant_list_tmp
done

mv 04_prune.keep.variant_list_tmp 04_prune.keep.variant_list
