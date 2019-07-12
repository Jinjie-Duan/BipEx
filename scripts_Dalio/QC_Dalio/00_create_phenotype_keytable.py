import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

PHENOFILE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/BIP_phenotype_information_cleaned_new_subtype_information_added.tsv'
PHENOTYPES_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/phenotypes.ht'

pheno_table = hl.import_table(PHENOFILE, quote="\"", key='SAMPLE_ALIAS',
	types={'PCT_CHIMERAS': hl.tfloat64, 'PCT_CONTAMINATION': hl.tfloat64})

pheno_table = pheno_table.drop(pheno_table.PHENOTYPE_FINE, pheno_table.PHENOTYPE_COARSE)

pheno_table = pheno_table.rename({'PHENOTYPE_COARSE_NEW': 'PHENOTYPE_COARSE', 'PHENOTYPE_FINE_NEW': 'PHENOTYPE_FINE'})

pheno_table.select(pheno_table.PROJECT_OR_COHORT, pheno_table.GENDER,
		pheno_table.PI, pheno_table.PHENOTYPE_COARSE, pheno_table.PHENOTYPE_FINE,
		pheno_table.LOCATION, pheno_table.INSTITUTION,
		pheno_table.PCT_CONTAMINATION, pheno_table.PCT_CHIMERAS).write(PHENOTYPES_TABLE, overwrite=True)

ht = hl.read_table(PHENOTYPES_TABLE)
n_samples = ht.count()

print('')
print('nSamples: ', '{:,}'.format(n_samples))
pprint(ht.describe())
