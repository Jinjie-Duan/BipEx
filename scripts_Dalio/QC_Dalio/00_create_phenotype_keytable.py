import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

BAM_METRICS = 'gs://dalio_bipolar_w1_w2/data/samples/bam_metrics.tsv'
PHENOFILE = 'gs://dalio_bipolar_w1_w2/data/samples/Dalio_phenotypes.tsv'
PHENOTYPES_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/phenotypes.ht'

ht_bam = hl.import_table(BAM_METRICS,
	key='Samples',
	types={'pct_chimeras': hl.tfloat64, 'pct_contamination': hl.tfloat64})

pheno_table = hl.import_table(PHENOFILE, key='Samples', quote="\"")

pheno_table = pheno_table.annotate(Phenotype = pheno_table.Primary_Disease_Coarse)
pheno_table.select(pheno_table.Project, pheno_table.Sex,
		pheno_table.PI, pheno_table.Phenotype, pheno_table.Primary_Disease,
		pheno_table.Location).join(ht_bam, how='left').write(PHENOTYPES_TABLE, overwrite=True)

ht = hl.read_table(PHENOTYPES_TABLE)
n_samples = ht.count()

print('')
print('nSamples: ', '{:,}'.format(n_samples))
pprint(ht.describe())
