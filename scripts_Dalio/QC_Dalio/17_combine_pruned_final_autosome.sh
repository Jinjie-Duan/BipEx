#!/bin/bash

cd /psych/genetics_data/dpalmer/Dalio/plink/pruned
cat 17_final_qc.chr1.prune.in 17_final_qc.chr2.prune.in > 17_prune.final_qc.keep.variant_list_tmp
# cat 17_final_qc_no_URVs.chr1.prune.in 17_final_qc_no_URVs.chr2.prune.in > 17_prune.final_qc_no_URVs.keep.variant_list_tmp

for i in `seq 3 22`;
do
	cat 17_prune.final_qc.keep.variant_list_tmp 17_final_qc.chr${i}.prune.in > 17_prune.final_qc.keep.variant_list
	mv 17_prune.final_qc.keep.variant_list 17_prune.final_qc.keep.variant_list_tmp

	# cat 17_prune.final_qc_no_URVs.keep.variant_list_tmp 17_final_qc_no_URVs.chr${i}.prune.in > 17_prune.final_qc_no_URVs.keep.variant_list
	# mv 17_prune.final_qc_no_URVs.keep.variant_list 17_prune.final_qc_no_URVs.keep.variant_list_tmp
done

mv 17_prune.final_qc.keep.variant_list_tmp 17_prune.final_qc.keep.variant_list
# mv 17_prune.final_qc_no_URVs.keep.variant_list_tmp 17_prune.final_qc_no_URVs.keep.variant_list
