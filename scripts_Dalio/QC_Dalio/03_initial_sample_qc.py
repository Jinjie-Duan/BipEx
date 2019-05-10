import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

# Inputs
MT = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/bipolar_wes_dalio_W1_W2/filterGT.mt'
INITIAL_VARIANT_LIST = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/02_prefilter.keep.variant_list'
PHENOTYPES_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/phenotypes.ht'

# Outputs
INITIAL_SAMPLE_QC_FILE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/03_initial_sample_qc.tsv'

variants_to_filter = hl.import_table(INITIAL_VARIANT_LIST,
	types={'locus':hl.tlocus(), 'alleles':hl.tarray(hl.tstr)})
variants_to_filter = variants_to_filter.key_by(locus=variants_to_filter.locus, alleles=variants_to_filter.alleles)

sample_annotations = hl.read_table(PHENOTYPES_TABLE)

mt = hl.read_matrix_table(MT)
mt = mt.filter_rows(hl.is_defined(variants_to_filter[mt.row_key]))
mt = mt.annotate_cols(phenotype = sample_annotations[mt.s])
mt = hl.sample_qc(mt, name='qc')

mt.cols().select('phenotype', 'qc').flatten().export(output=INITIAL_SAMPLE_QC_FILE)
