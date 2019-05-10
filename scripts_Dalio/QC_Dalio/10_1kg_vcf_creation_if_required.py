from __future__ import print_function
from hail import *

hc = HailContext(log='/tmp/hail.log')

# This is a hail 0.1 script used to generate a vcf for import into hail 0.2

VDS_1KG = 'gs://raw_data_bipolar/data/ALL.1KG.qc.hardcalls.vds'

vds_1kg = hc.read(VDS_1KG)
vds_1kg.export_vcf('gs://raw_data_bipolar/data/ALL.1KG.qc.hardcalls.vcf.bgz')
