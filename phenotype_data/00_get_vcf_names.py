import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

RAW_MT = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/dalio_bipolar_w1_w2/Dalio_W1_W2_GRCh38_exomes.mt'

mt = hl.read_matrix_table(RAW_MT)

mt.cols().select().export('gs://raw_data_bipolar_dalio_w1_w2_hail_02/dalio_bipolar_w1_w2/Dalio_W1_W2_GRCh38_exomes_sample_IDs.tsv')