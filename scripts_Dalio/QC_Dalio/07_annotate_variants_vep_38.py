import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

MT = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/bipolar_wes_dalio_W1_W2/filterGT_GRCh38_6_multi.mt'

ht = hl.read_matrix_table(MT).rows()

ht_vep = hl.vep(ht, "gs://hail-common/vep/vep/vep95-GRCh38-loftee-gcloud.json")
ht_vep.write("gs://dalio_bipolar_w1_w2_hail_02/data/annotations/vep_annotate.ht", overwrite=True)
