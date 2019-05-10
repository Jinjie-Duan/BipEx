import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

# Names of .mt files.
MT_HARDCALLS = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/bipolar_wes_dalio_W1_W2/filterGT.hardcalls.mt'

# Read in the hard calls matrix table.
mt = hl.read_matrix_table(MT_HARDCALLS)

# Low complexity regions in the data.
LCRs = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/inputs/inputs_low_complexity_regions.interval_list'

INITIAL_VARIANT_QC_FILE  = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/02_prefilter_metrics.tsv'
INITIAL_VARIANT_LIST = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/02_prefilter.keep.variant_list'

# Import the interval lists for the LCRs.
LCR_intervals = hl.import_locus_intervals(LCRs)

# Annotate variants with flag indicating if they are in LCR or failed VQSR.
mt = mt.annotate_rows(fail_VQSR = hl.len(mt.filters) != 0)
mt = mt.annotate_rows(in_LCR = hl.is_defined(LCR_intervals[mt.locus]))

# Get information about the number of variants that were excluded.
fail_VQSR = mt.filter_rows(mt.fail_VQSR).count_rows()
in_LCR = mt.filter_rows(mt.in_LCR).count_rows()

print('nVariants failing VQSR:')
pprint(fail_VQSR)
print('nVariants in low complexity regions:')
pprint(in_LCR)

# Export variant annotations, and variant keytable.
mt_rows = mt.rows()
mt_rows.select(mt_rows.fail_VQSR, mt_rows.in_LCR).export(INITIAL_VARIANT_QC_FILE)

mt = mt.filter_rows(mt.fail_VQSR | mt.in_LCR, keep=False)
mt_rows_filter = mt.rows().select().export(INITIAL_VARIANT_LIST)

n_variants = hl.import_table(INITIAL_VARIANT_LIST).count()

print('nVariants after initial filter:')
print(n_variants)
