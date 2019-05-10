import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

MT_HARDCALLS = 'gs://raw_data_bipolar_dalio_w1_w2_hail_02/bipolar_wes_dalio_W1_W2/filterGT.hardcalls.mt'
GNOMAD_NONPSYCH_VARIANTS = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/gnomad.r2.0.1.nonpsych.variant_list'
ANNOTATION_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/data/annotations/gene.ht'

# MPC score.
MPC_SCORE = 'gs://raw_data_bipolar_dalio_w1_w2/inputs/fordist_constraint_official_mpc_values_v2.txt.bgz'

mpc_ht = hl.import_table(MPC_SCORE, impute=True)
mpc_ht = mpc_ht.annotate(locus = hl.locus(contig = mpc_ht.chrom, pos = mpc_ht.pos),
    alleles = [mpc_ht.ref, mpc_ht.alt])
mpc_ht = mpc_ht.key_by(mpc_ht.locus, mpc_ht.alleles).select('MPC')

# CADD
cadd_ht = hl.read_table("gs://hail-datasets-hail-data/CADD.1.4.GRCh37.ht")
# constraint information
constraint_ht = hl.import_table("gs://annotationdb/gene/constraint.tsv.bgz", impute=True).key_by('transcript')
# DiscovEHR
discovEHR_ht = hl.import_vcf("gs://annotationdb/discovEHR/GHS_Freeze_50.L3DP10.pVCF.frq.vcf.bgz")

mt = hl.read_matrix_table(MT_HARDCALLS).repartition(2048)

# mt = hl.vep(mt, "gs://raw_data_bipolar/inputs/vep_configuration.json")
# The following vep configuration file is public one hosted by hail.
mt = hl.vep(mt, "gs://hail-common/vep/vep/vep85-loftee-gcloud.json")
mt = mt.annotate_rows(mpc = mpc_ht[mt.row_key])
mt = mt.annotate_rows(cadd = cadd_ht[mt.row_key])
mt = mt.annotate_rows(discovEHR = discovEHR_ht.rows()[mt.row_key])

ptv = hl.set(["transcript_ablation", "splice_acceptor_variant",
              "splice_donor_variant", "stop_gained", "frameshift_variant"])

missense = hl.set(["stop_lost", "start_lost", "transcript_amplification",
                   "inframe_insertion", "inframe_deletion", "missense_variant",
                   "protein_altering_variant", "splice_region_variant"])

synonymous = hl.set(["incomplete_terminal_codon_variant", "stop_retained_variant", "synonymous_variant"])

non_coding = hl.set(["coding_sequence_variant", "mature_miRNA_variant", "5_prime_UTR_variant",
              "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant",
              "NMD_transcript_variant", "non_coding_transcript_variant", "upstream_gene_variant",
              "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant",
              "regulatory_region_ablation", "regulatory_region_amplification", "feature_elongation",
              "regulatory_region_variant", "feature_truncation", "intergenic_variant"])

mt = mt.annotate_rows(gene_transcript = hl.rbind( 
    hl.filter(mt.vep.transcript_consequences.filter(lambda tc: hl.set(tc.consequence_terms).contains(mt.vep.most_severe_consequence))),
    lambda transcripts: hl.or_else(hl.find(lambda x: x.canonical == 1, transcripts), transcripts[0])))

mt = mt.annotate_rows(consequence_category = 
    hl.case().when(ptv.contains(mt.vep.most_severe_consequence), "ptv")
             .when(missense.contains(mt.vep.most_severe_consequence) & 
                   (~hl.is_defined(mt.gene_transcript.polyphen_prediction) | 
                    ~hl.is_defined(mt.gene_transcript.sift_prediction) ), "other_missense")
             .when(missense.contains(mt.vep.most_severe_consequence) & 
                   (mt.gene_transcript.polyphen_prediction == "probably_damaging") & 
                   (mt.gene_transcript.sift_prediction == "deleterious"), "damaging_missense")
             .when(missense.contains(mt.vep.most_severe_consequence), "other_missense")
             .when(synonymous.contains(mt.vep.most_severe_consequence), "synonymous")
             .when(non_coding.contains(mt.vep.most_severe_consequence), "non_coding")
             .default("NA")
    )

# Now that we have the gene information for each variant, we can evaluate the constraint.

gnomad_nonpsych_variants_ht = hl.import_table(GNOMAD_NONPSYCH_VARIANTS, no_header=True)

gnomad_nonpsych_variants_ht = gnomad_nonpsych_variants_ht.annotate(**hl.parse_variant(gnomad_nonpsych_variants_ht.f0))
gnomad_nonpsych_variants_ht = gnomad_nonpsych_variants_ht.key_by(gnomad_nonpsych_variants_ht.locus, gnomad_nonpsych_variants_ht.alleles)

mt = mt.annotate_rows(inGnomAD_nonpsych = hl.is_defined(gnomad_nonpsych_variants_ht[mt.row_key]))
mt = mt.annotate_rows(gene_constraint = constraint_ht[mt.gene_transcript.transcript_id])

mt_rows = mt.rows().repartition(64)
mt_rows.write(ANNOTATION_TABLE, overwrite=True)



