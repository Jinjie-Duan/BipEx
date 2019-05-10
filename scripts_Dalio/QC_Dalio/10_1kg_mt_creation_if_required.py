import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

ds = hl.import_vcf('gs://raw_data_bipolar/data/ALL.1KG.qc.hardcalls.vcf.bgz', reference_genome='GRCh37')
ds.write(output='gs://raw_data_bipolar/data/ALL.1KG.qc.hardcalls.mt', overwrite=True)
