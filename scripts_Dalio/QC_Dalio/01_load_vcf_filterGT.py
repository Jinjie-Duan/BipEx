import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

RAW_MT = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/dalio_bipolar_w1_w2/Dalio_W1_W2_exomes.mt'
MT = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/bipolar_wes_dalio_W1_W2/filterGT.mt'
MT_HARDCALLS = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/bipolar_wes_dalio_W1_W2/filterGT.hardcalls.mt'

mt = hl.read_matrix_table(RAW_MT)
mt = hl.split_multi_hts(mt)

mt = mt.filter_entries(
    hl.is_defined(mt.GT) & (
        (mt.GT.is_hom_ref() & ( ((mt.AD[0] / mt.DP) < 0.8) | (mt.GQ < 20) | (mt.DP < 10)) ) |
        (mt.GT.is_het() & ( ((mt.AD[0] + mt.AD[1]) / mt.DP < 0.8) | ((mt.AD[1] / mt.DP) < 0.2) | (mt.PL[0] < 20) | (mt.DP < 10) )) |
        (mt.GT.is_hom_var() & ( ((mt.AD[1] / mt.DP) < 0.8) | (mt.PL[0] < 20) | (mt.DP < 10) ))
    ), keep = False)

mt.write(MT, overwrite=True)
mt.select_entries(mt.GT).repartition(512).write(MT_HARDCALLS, overwrite=True)

mt = hl.read_matrix_table(MT_HARDCALLS)
n = mt.count()

pprint('nSamples:')
print(n[0])
pprint('nVariants:')
print(n[1])
