import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

ANNOTATION_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/data/annotations/gene.ht'

ht = hl.read_table(ANNOTATION_TABLE)
ht = ht.annotate(type = hl.case().when((ht.alleles[0].length() == 1) | (ht.alleles[1].length() == 1), "SNP")
	.when(ht.alleles[0].length() < ht.alleles[1].length(), "Insertion")
	.when(ht.alleles[0].length() > ht.alleles[1].length(), "Deletion")
	.default("No type"),
	inDiscovEHR = hl.is_defined(ht.discovEHR.info.AF))

n = ht.count()
n_inDiscovEHR = ht.aggregate(hl.agg.counter(ht.inDiscovEHR))
n_inGnomAD_nonpsych = ht.aggregate(hl.agg.counter(ht.inGnomAD_nonpsych))
n_have_gene = ht.aggregate(hl.agg.counter(hl.is_defined(ht.gene_transcript.gene_symbol)))
n_consequence_category = ht.aggregate(hl.agg.counter(ht.consequence_category))
n_most_severe_consequence = ht.aggregate(hl.agg.counter(ht.vep.most_severe_consequence))

print('nVariants:')
print(n[0])

print('discovEHR')
print(n_inDiscovEHR)

print('gnomAD_nonpsych')
print(n_inGnomAD_nonpsych)

print('have_gene')
print(n_have_gene)

print('consequence_category')
print(n_consequence_category)

print('most_severe_consequence')
print(n_most_severe_consequence)
