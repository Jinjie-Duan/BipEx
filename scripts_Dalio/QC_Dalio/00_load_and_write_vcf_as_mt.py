import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

ds = hl.import_vcf('gs://raw_data_bipolar_dalio_w1_w2/bipolar_wes_dalio_W1_W2/BiPolar_CasesControls_Exomes.vcf.bgz', reference_genome='GRCh37')
ds.write(output='gs://raw_data_bipolar_dalio_w1_w2_hail_02/dalio_bipolar_w1_w2/Dalio_W1_W2_exomes.mt', overwrite=True)
