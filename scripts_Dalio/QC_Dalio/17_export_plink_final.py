import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

MT_HARDCALLS = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/bipolar_wes_dalio_W1_W2/filterGT.hardcalls.mt'
FINAL_PLINK_FILES = 'gs://dalio_bipolar_w1_w2_hail_02/data/plink/final_qc'

HIGH_LD_INTERVALS = 'gs://raw_data_bipolar_dalio_w1_w2/inputs/high_LD_regions_b37.interval_list'
FINAL_SAMPLE_LIST = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/16_final_qc.keep.sample_list'
# FINAL_SAMPLE_LIST_REMOVE_SINGLETON_OUTLIERS = 'gs://dalio_bipolar_w1_w2/data/samples/16_final_qc_remove_singleton_outliers.keep.sample_list'
FINAL_VARIANT_LIST = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/15_final_qc.keep.variant_list'

ht_final_samples = hl.import_table(FINAL_SAMPLE_LIST, no_header=True, key='f0')
ht_final_variants = hl.import_table(FINAL_VARIANT_LIST, types={'locus':hl.tlocus(), 'alleles':hl.tarray(hl.tstr)})
ht_final_variants = ht_final_variants.key_by(ht_final_variants.locus, ht_final_variants.alleles)

high_LD_intervals = hl.import_locus_intervals(HIGH_LD_INTERVALS)

mt = hl.read_matrix_table(MT_HARDCALLS)
mt = mt.filter_cols(hl.is_defined(ht_final_samples[mt.col_key]))
mt = mt.annotate_rows(in_high_LD = hl.is_defined(high_LD_intervals[mt.locus]))

mt = mt.filter_rows(hl.is_defined(ht_final_variants[mt.row_key]) & (~mt.in_high_LD))
mt = mt.filter_rows(mt.locus.in_x_nonpar() | mt.locus.in_autosome_or_par())
mt = hl.variant_qc(mt, name='qc')
mt = mt.filter_rows((mt.qc.AF[0] > 0.01) & (mt.qc.AF [0]< 0.99) & ((mt.qc.call_rate > 0.98) | mt.locus.in_x_nonpar() | mt.locus.in_x_par())).persist()

def rename_samples(mt, mapping):
    return mt.key_cols_by(s = hl.literal(mapping).get(mt.s, default=mt.s))

mt = rename_samples(mt, {'431-BG00852 D':'431-BG00852_D'})

for x in range(1,22):

	mt_chr = hl.filter_intervals(mt, [hl.parse_locus_interval(hl.str(x))])
	n_chr = mt_chr.count_rows()

	print('nVariants in chr')
	print(x)
	print(n_chr)

	hl.export_plink(mt_chr, FINAL_PLINK_FILES + '.chr' + str(x))

mt_chr = hl.filter_intervals(mt, [hl.parse_locus_interval('X')])
n_chr = mt_chr.count_rows()

print('nVariants in chr')
print('X')
print(n_chr)

hl.export_plink(mt_chr, FINAL_PLINK_FILES + '.chr' + 'X')
