#!/bin/bash

perl INSTALL.pl --AUTO ac -s homo_sapiens_refseq --ASSEMBLY GRCh37
perl INSTALL.pl --AUTO f -s homo_sapiens_refseq --ASSEMBLY GRCh37



VEP_CACHE="/home/gonzalo/.vep/"
VEP_REF="/home/gonzalo/.vep/homo_sapiens_refseq/105_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz"

perl ~/software/ensembl-vep \
--cache --offline --hgvs --refseq --dir $VEP_CACHE --dir_plugins $PLUGIN_DIR 





docker run -t -i -v ${VEP_CACHE}:/opt/vep/.vep ensemblorg/ensembl-vep ./vep --help


VEP_CACHE="/home/gonzalo/.vep/"
INPUT_DIR="/mnt/genetica4/gonzalo/pruebas_vep/"
threads="8"

docker run -t -i -u $(id -u):$(id -g) \
-v ${VEP_CACHE}:/opt/vep/.vep \
-v ${INPUT_DIR}:/opt/vep/working_dir \
ensemblorg/ensembl-vep ./vep \
--cache --offline --hgvs --refseq --dir /opt/vep/.vep \
-v --fork ${threads} --assembly GRCh37 --force_overwrite  --no_stats --tab --format vcf \
--fasta /opt/vep/.vep/homo_sapiens_refseq/105_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
--input_file /opt/vep/working_dir/head_21-4834-CES178.final.vcf \
--output_file /opt/vep/working_dir/my_output.tsv


docker run -t -i -u $(id -u):$(id -g) \
-v ${VEP_CACHE}:/opt/vep/.vep \
-v ${INPUT_DIR}:/opt/vep/working_dir \
ensemblorg/ensembl-vep ./vep \
--cache --offline --hgvs --refseq --dir /opt/vep/.vep \
-v --fork ${threads} --assembly GRCh37 --force_overwrite  --no_stats --tab --format vcf --check_existing \
--fasta /opt/vep/.vep/homo_sapiens_refseq/105_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
--input_file /opt/vep/working_dir/21-4834-CES178.final.vcf \
--output_file /opt/vep/working_dir/21-4834-CES178.final.annotated.tsv


docker run -t -i -u $(id -u):$(id -g) \
-v ${VEP_CACHE}:/opt/vep/.vep \
-v ${INPUT_DIR}:/opt/vep/working_dir \
ensemblorg/ensembl-vep ./vep \
--cache --offline --hgvs --refseq --dir /opt/vep/.vep \
-v --fork ${threads} --assembly GRCh37 --force_overwrite  --no_stats --vcf --format vcf \
--fasta /opt/vep/.vep/homo_sapiens_refseq/105_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
--input_file /opt/vep/working_dir/21-4834-CES178.final.vcf \
--output_file /opt/vep/working_dir/21-4834-CES178.final.annotated.vcf




docker run -t -i ensemblorg/ensembl-vep perl -e 'use Bio::DB::HTS::Faidx'




VEP_CACHE="/home/gonzalo/.vep/"
INPUT_DIR="/mnt/genetica4/gonzalo/pruebas_vep/"
threads="8"


docker run -t -i -u $(id -u):$(id -g) \
-v ${VEP_CACHE}:/opt/vep/.vep \
-v ${INPUT_DIR}:/opt/vep/working_dir \
ensemblorg/ensembl-vep ./vep \
--cache --offline --hgvs --refseq --dir /opt/vep/.vep \
-v --fork ${threads} --assembly GRCh37 --force_overwrite  --no_stats --tab --format vcf --check_existing \
--individual all --use_transcript_ref \
--fasta /opt/vep/.vep/homo_sapiens_refseq/105_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
--input_file /opt/vep/working_dir/head_21-4834-CES178.final.vcf \
--output_file /opt/vep/working_dir/head_21-4834-CES178.final.annotated.tsv \
--custom /opt/vep/working_dir/vcf_to_annotate.vcf.gz,,vcf,exact,0


docker run -t -i -u $(id -u):$(id -g) \
-v ${VEP_CACHE}:/opt/vep/.vep \
-v ${INPUT_DIR}:/opt/vep/working_dir \
ensemblorg/ensembl-vep ./vep \
--cache --offline --hgvs --refseq --dir /opt/vep/.vep \
-v --fork ${threads} --assembly GRCh37 --force_overwrite  --no_stats --vcf --format vcf --check_existing \
--individual all --use_transcript_ref \
--fasta /opt/vep/.vep/homo_sapiens_refseq/105_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
--input_file /opt/vep/working_dir/head_21-4834-CES178.final.vcf \
--output_file /opt/vep/working_dir/head_21-4834-CES178.final.annotated.vcf







input_vcf="/mnt/genetica4/gonzalo/ejemplo_trio2/CES_130_30-07-2020/snvs/CES_130_30-07-2020_2021_11_26_16_50_09.final.vcf"
input_vcf="/mnt/genetica4/gonzalo/pruebas_vep/head_21-4834-CES178.final.vcf"


# Header
bcftools view -h ${input_vcf} | grep "##" > vcf_to_annotate.vcf
for sample in $(bcftools query -l ${input_vcf})
do
echo "##INFO=<ID=${sample}_GT,Number=.,Type=String,Description=\"${sample} Genotype\">" >> vcf_to_annotate.vcf
echo "##INFO=<ID=${sample}_AD,Number=.,Type=String,Description=\"${sample} Allelic depths for the ref and alt alleles in the order listed\">" >> vcf_to_annotate.vcf
echo "##INFO=<ID=${sample}_DP,Number=.,Type=String,Description=\"${sample} Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">" >> vcf_to_annotate.vcf
done
bcftools view -h ${input_vcf} | grep "#CHROM" | cut -f1-8 >> vcf_to_annotate.vcf

# NEW INFO
bcftools query -f '[;%SAMPLE\_GT=%GT][;%SAMPLE\_AD=%AD][;%SAMPLE\_DP=%DP]\n' ${input_vcf} | sed 's/;//1' | sed 's/,/_/g' > new_info.txt
bcftools view -H ${input_vcf} | cut -f1-7 > old_info.txt
paste -d '\t' old_info.txt new_info.txt >> vcf_to_annotate.vcf
rm new_info.txt old_info.txt

bgzip -c vcf_to_annotate.vcf > vcf_to_annotate.vcf.gz
tabix -p vcf vcf_to_annotate.vcf.gz


fields=""
for sample in $(bcftools query -l ${input_vcf}); do fields="$(echo "${fields},${sample}_GT,${sample}_AD,${sample}_DP")"; done




VEP_CACHE="/home/gonzalo/.vep/"
INPUT_DIR="/mnt/genetica4/gonzalo/pruebas_vep/"
threads="8"
PLUGIN_DIR="/mnt/genetica7/vep_plugins/"



# hg38
dbscSNV="dbscSNV1.1_GRCh38.txt.gz"
LoFtool="LoFtool_scores.txt"
ExACpLI="ExACpLI_values.txt"
dbNSFP="dbNSFP4.3a_grch38.gz"
MaxEntScan="maxEntScan"
CADD_INDELS="CADD_v.6_hg38_InDels.tsv.gz"
CADD_SNVS="CADD_v.6_hg38_gnomad.genomes.r3.0.indel.tsv.gz"
Kaviar="Kaviar-160204-Public-hg38.vcf.gz"
CCRS_DB="ccrs.all.v2.20180420.hg38.bed.gz"
DENOVO_DB="denovo-db.non-ssc-samples.variants.hg38.chr.vcf.gz"
CLINVAR="clinvar.GRCh38.vcf.gz"
GNOMADg="gnomad.genomes.v3.1.2.hg38.vcf.bgz"
GNOMADe="gnomad.exomes.r2.1.1.hg38.vcf.bgz"
GNOMADg_cov="gnomad.genomes.v3.0.1.coverage.hg38.positions.vcf.gz"
GNOMADe_cov="gnomad.exomes.r2.1.1.coverage.hg38.positions.vcf.gz"
CSVS="CSVS_vcf/all.hg38.vcf.gz"
MutScore="mutscore-v1.0-hg38.vcf.gz"
MAF_FJD_COHORT="MAFdb_AN20_latest.hg38.vcf.gz"
SpliceAI_SNV="spliceai_scores.raw.snv.hg38.vcf.gz"
# SpliceAI_INDEL="spliceai_scores.raw.indel.hg38.vcf.gz"


# hg19
dbscSNV="dbscSNV1.1_GRCh37.txt.gz"
LoFtool="LoFtool_scores.txt"
ExACpLI="ExACpLI_values.txt"
dbNSFP="dbNSFP4.3a_grch37.gz"
MaxEntScan="maxEntScan"
CADD_INDELS="CADD_v.6_hg19_InDels.tsv.gz"
CADD_SNVS="CADD_v.6_hg19_whole_genome_SNVs.tsv.gz"
Kaviar="Kaviar-160204-Public-hg19.vcf.gz"
CCRS_DB="ccrs.all.v2.20180420.hg19.bed.gz"
DENOVO_DB="denovo-db.non-ssc-samples.variants.hg19.chr.vcf.gz"
CLINVAR="clinvar.GRCh37.vcf.gz"
GNOMADg="gnomad.genomes.v3.1.2.hg19.vcf.bgz"
GNOMADe="gnomad.exomes.r2.1.1.hg19.vcf.bgz"
GNOMADg_cov="gnomad.genomes.v3.0.1.coverage.hg37.positions.vcf.gz"
GNOMADe_cov="gnomad.exomes.r2.1.1.coverage.hg37.positions.vcf.gz"
CSVS="CSVS_vcf/all.hg19.vcf.gz"
MutScore="mutscore-v1.0-hg19.vcf.gz"
MAF_FJD_COHORT="MAFdb_AN20_latest.hg37.vcf.gz"
SpliceAI_SNV="spliceai_scores.raw.snv.hg19.vcf.gz"
#SpliceAI_INDEL="spliceai_scores.raw.indel.hg19.vcf.gz"


SpliceAI_INDEL="head_21-4834-CES178.final.INDEL.spliceai.vcf.gz"
SAMPLE_INFO="vcf_to_annotate.vcf.gz"


docker run -t -i -u $(id -u):$(id -g) \
-v ${VEP_CACHE}:/opt/vep/.vep \
-v ${INPUT_DIR}:/opt/vep/working_dir \
-v ${PLUGIN_DIR}:/opt/vep/Plugins \
ensemblorg/ensembl-vep ./vep \
--cache --offline --dir_cache /opt/vep/.vep --dir_plugins /opt/vep/.vep/Plugins \
--refseq --species homo_sapiens --assembly GRCh37 --force_overwrite --use_transcript_ref \
--verbose --fork ${threads} --no_stats --tab --format vcf \
--fasta /opt/vep/.vep/homo_sapiens_refseq/105_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
--input_file /opt/vep/working_dir/head_21-4834-CES178.final.vcf \
--output_file /opt/vep/working_dir/head_21-4834-CES178.final.annotated.15.tsv \
--check_existing --canonical --numbers --hgvs --biotype --regulatory --symbol --protein \
--sift p --polyphen p --allele_number --variant_class --pubmed \
--plugin dbscSNV,/opt/vep/Plugins/${dbscSNV} \
--plugin LoFtool,/opt/vep/Plugins/${LoFtool} \
--plugin ExACpLI,/opt/vep/Plugins/${ExACpLI} \
--plugin dbNSFP,/opt/vep/Plugins/${dbNSFP},\
LRT_pred,M-CAP_pred,MetaLR_pred,MetaSVM_pred,MutationAssessor_pred,MutationTaster_pred,PROVEAN_pred,\
FATHMM_pred,MetaRNN_pred,PrimateAI_pred,DEOGEN2_pred,BayesDel_addAF_pred,BayesDel_noAF_pred,ClinPred_pred,\
LIST-S2_pred,Aloft_pred,fathmm-MKL_coding_pred,fathmm-XF_coding_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,\
phyloP30way_mammalian,phastCons30way_mammalian,\
GERP++_RS,Interpro_domain,GTEx_V8_gene,GTEx_V8_tissue \
--plugin MaxEntScan,/opt/vep/Plugins/${MaxEntScan} \
--plugin CADD,/opt/vep/Plugins/${CADD_INDELS},/opt/vep/Plugins/${CADD_SNVS} \
--custom /opt/vep/Plugins/${Kaviar},kaviar,vcf,exact,0,AF,AC,AN \
--custom /opt/vep/Plugins/${CCRS_DB},gnomAD_exomes_CCR,bed,overlap,0 \
--custom /opt/vep/Plugins/${DENOVO_DB},denovoVariants,vcf,exact,0,SAMPLE_CT \
--custom /opt/vep/Plugins/${CLINVAR},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
--custom /opt/vep/Plugins/${GNOMADg},gnomADg,vcf,exact,0,AF,AC,AN,nhomalt,popmax,AF_popmax,AC_popmax,AF_nfe,AC_nfe \
--custom /opt/vep/Plugins/${GNOMADe},gnomADe,vcf,exact,0,AF,AC,AN,nhomalt,popmax,AF_popmax,AC_popmax,AF_nfe,AC_nfe \
--custom /opt/vep/Plugins/${GNOMADg_cov},gnomADg_cov,vcf,overlap,0,median,perc_20x \
--custom /opt/vep/Plugins/${GNOMADe_cov},gnomADe_cov,vcf,overlap,0,median,perc_20x \
--custom /opt/vep/Plugins/${CSVS},CSVS,vcf,exact,0,AF,AC, \
--custom /opt/vep/Plugins/${MAF_FJD_COHORT},FJD_MAF,vcf,exact,0,AF,AC \
--custom /opt/vep/Plugins/${MutScore},Mut,vcf,exact,0,Score \
--custom /opt/vep/Plugins/${SpliceAI_SNV},SpliceAI_SNV,vcf,exact,0,SpliceAI \
--custom /opt/vep/working_dir/${SpliceAI_INDEL},SpliceAI_INDEL,vcf,exact,0,SpliceAI \
--custom /opt/vep/working_dir/${SAMPLE_INFO},SAMPLE,vcf,exact,0${fields}


head -n 1000 ${INPUT_DIR}/head_21-4834-CES178.final.annotated.13.tsv | grep "#Uploaded_variation" -n | sed 's/:.*//'


#==========#
# SpliceAI #
#==========#
REFERENCE_DIR="/mnt/genetica7/references/"
WORKING_DIR="/mnt/genetica4/gonzalo/pruebas_vep/"

docker run -it \
-u $(id -u):$(id -g) \
-v ${WORKING_DIR}:/mnt/io_dir \
-v ${REFERENCE_DIR}:/mnt/references \
broadinstitute/gatk \
./gatk SelectVariants --tmp-dir /mnt/io_dir \
-R /mnt/references/hg19.fa \
-V /mnt/io_dir/head_21-4834-CES178.final.vcf \
--select-type-to-include INDEL \
-O /mnt/io_dir/head_21-4834-CES178.final.INDEL.vcf

docker run -it -u $(id -u):$(id -g) \
-v ${REFERENCE_DIR}:/home/ref/ \
-v ${WORKING_DIR}:/home/working_dir/ \
skysbiodocker/spliceai \
spliceai -R /home/ref/hg19.fa -A grch37 \
-I /home/working_dir/head_21-4834-CES178.final.INDEL.vcf \
-O /home/working_dir/head_21-4834-CES178.final.INDEL.spliceai.vcf

bgzip head_21-4834-CES178.final.INDEL.spliceai.vcf
tabix -p vcf head_21-4834-CES178.final.INDEL.spliceai.vcf.gz



#=========#
# AutoMap #
#=========#
docker run -it -u $(id -u):$(id -g) \
-v ${WORKING_DIR}:/home/working_dir/ \
automap bash AutoMap_v1.0.sh \
--vcf /home/working_dir/21-4834-CES178.final.vcf \
--out /home/working_dir/ --genome hg19




##########
# QUITAR #
##########

# Solo se annota con merge:
# --uniprot --domains --gene_phenotype

# Ya están anotadas las frecuancias con gnomad exomas y genomas
# dbNSFP gnomAD_exomes_AF,gnomAD_exomes_NFE_AF,1000Gp3_AF,1000Gp3_EUR_AF,ExAC_AF,ExAC_EAS_AF,ExAC_NFE_AF,ExAC_Adj_AF,

# Information is already annotated in Existing_variation 
# dbNSFP rs_dbSNP150

# phyloP30way_mammalian_rankscore,phastCons30way_mammalian_rankscore,GERP++_RS_rankscore,

##########
# AÑADIR #
##########

# dbNSFP
# El aa de referencia y alternativo y posición está en HGVSp
5	aaref: reference amino acid
		"." if the variant is a splicing site SNP (2bp on each end of an intron)
6	aaalt: alternative amino acid
		"." if the variant is a splicing site SNP (2bp on each end of an intron)
12	aapos: amino acid position as to the protein.
		"-1" if the variant is a splicing site SNP (2bp on each end of an intron). 
		Multiple entries separated by ";", corresponding to Ensembl_proteinid

14	Ensembl_geneid: Ensembl gene id
15	Ensembl_transcriptid: Ensembl transcript ids (Multiple entries separated by ";")
16	Ensembl_proteinid: Ensembl protein ids
		Multiple entries separated by ";",  corresponding to Ensembl_transcriptids
17	Uniprot_acc: Uniprot accession number matching the Ensembl_proteinid
		Multiple entries separated by ";".
30	refcodon: reference codon


68	VEST4_score: VEST 4.0 score. Score ranges from 0 to 1. The larger the score the more likely
		the mutation may cause functional change. 
		Multiple scores separated by ";", corresponding to Ensembl_transcriptid.
		Please note this score is free for non-commercial use. For more details please refer to 
		http://wiki.chasmsoftware.org/index.php/SoftwareLicense. Commercial users should contact 
		the Johns Hopkins Technology Transfer office.
79	MetaRNN_pred: Prediction of our MetaRNN based ensemble prediction score,"T(olerated)" or
		"D(amaging)". The score cutoff between "D" and "T" is 0.5. The rankscore cutoff between 
		"D" and "T" is 0.6149.
83	REVEL_score: REVEL is an ensemble score based on 13 individual scores for predicting the
		pathogenicity of missense variants. Scores range from 0 to 1. The larger the score the more 
		likely the SNP has damaging effect. "REVEL scores are freely available for non-commercial use.  
		For other uses, please contact Weiva Sieh" (weiva.sieh@mssm.edu)
		Multiple entries are separated by ";", corresponding to Ensembl_transcriptid.
90	MVP_score: A pathogenicity prediction score for missense variants using deep learning approach.
		The range of MVP score is from 0 to 1. The larger the score, the more likely the variant is 
		pathogenic. The authors suggest thresholds of 0.7 and 0.75 for separating damaging vs tolerant 
		variants in constrained genes (ExAC pLI >=0.5) and non-constrained genes (ExAC pLI<0.5), respectively. 
		Details see doi: http://dx.doi.org/10.1101/259390
		Multiple entries are separated by ";", corresponding to Ensembl_transcriptid.
92	MPC_score: A deleteriousness prediction score for missense variants based on regional missense
		constraint. The range of MPC score is 0 to 5. The larger the score, the more likely the variant is 
		pathogenic. Details see doi: http://dx.doi.org/10.1101/148353.
		Multiple entries are separated by ";", corresponding to Ensembl_transcriptid.
96	PrimateAI_pred: Prediction of PrimateAI score based on the authors recommendation, "T(olerated)" or
		"D(amaging)". The score cutoff between "D" and "T" is 0.80.
99	DEOGEN2_pred: Prediction of DEOGEN2 score based on the authors recommendation, "T(olerated)" or
		"D(amaging)". The score cutoff between "D" and "T" is 0.5.
102	BayesDel_addAF_pred: Prediction of BayesDel_addAF score based on the authors recommendation, "T(olerated)" or
		"D(amaging)". The score cutoff between "D" and "T" is 0.0692655.
105	BayesDel_noAF_pred: Prediction of BayesDel_noAF score based on the authors recommendation, "T(olerated)" or
		"D(amaging)". The score cutoff between "D" and "T" is -0.0570105.		
108	ClinPred_pred: Prediction of ClinPred score based on the authors recommendation, "T(olerated)" or "D(amaging)".
		The score cutoff between "D" and "T" is 0.5.
111	LIST-S2_pred: Prediction of LIST-S2 score based on the authors recommendation, "T(olerated)" or "D(amaging)".
		The score cutoff between "D" and "T" is 0.85.
116	Aloft_pred: final classification predicted by ALoFT;
		values can be Tolerant, Recessive or Dominant
		multiple values separated by ";", corresponding to Ensembl_proteinid.
128	fathmm-MKL_coding_pred: If a fathmm-MKL_coding_score is >0.5 (or rankscore >0.28317) 
		the corresponding nsSNV is predicted as "D(AMAGING)"; otherwise it is predicted as "N(EUTRAL)".
132	fathmm-XF_coding_pred: If a fathmm-XF_coding_score is >0.5, the corresponding nsSNV is predicted
		as "D(AMAGING)"; otherwise it is predicted as "N(EUTRAL)".
632	clinvar_clnsig: clinical significance by clinvar
		Possible values: Benign, Likely_benign, Likely_pathogenic, Pathogenic, drug_response, 
		histocompatibility. A negative score means the score is for the ref allele
633	clinvar_trait: the trait/disease the clinvar_clnsig referring to
634	clinvar_review: ClinVar Review Status summary
		Possible values:  no assertion criteria provided, criteria provided, single submitter,
		criteria provided, multiple submitters, no conflicts, reviewed by expert panel, practice guideline


#################
# Documentación #
#################
phyloP20way_mammalian
160	phyloP30way_mammalian: phyloP (phylogenetic p-values) conservation score based on the
		multiple alignments of 30 mammalian genomes (including human). The larger the score, 
		the more conserved the site. Scores range from -20 to 1.312 in dbNSFP.

phyloP20way_mammalian_rankscore
161	phyloP30way_mammalian_rankscore: phyloP30way_mammalian scores were ranked among all
		phyloP30way_mammalian scores in dbNSFP. The rankscore is the ratio of the rank of the 
		score over the total number of phyloP30way_mammalian scores in dbNSFP.

phastCons20way_mammalian
166	phastCons30way_mammalian: phastCons conservation score based on the multiple alignments
		of 30 mammalian genomes (including human). The larger the score, the more conserved 
		the site. Scores range from 0 to 1. 

phastCons20way_mammalian_rankscore
167	phastCons30way_mammalian_rankscore: phastCons30way_mammalian scores were ranked among
		all phastCons30way_mammalian scores in dbNSFP. The rankscore is the ratio of the rank 
		of the score over the total number of phastCons30way_mammalian scores in dbNSFP.

GERP++_RS
156	GERP++_RS: GERP++ RS score, the larger the score, the more conserved the site. Scores range from
		-12.3 to 6.17.

GERP++_RS_rankscore
157	GERP++_RS_rankscore: GERP++ RS scores were ranked among all GERP++ RS scores in dbNSFP.
		The rankscore is the ratio of the rank of the score over the total number of GERP++ RS 
		scores in dbNSFP.

LRT_pred 
52	LRT_pred: LRT prediction, D(eleterious), N(eutral) or U(nknown), which is not solely
		determined by the score.

MutationTaster_pred
56	MutationTaster_pred: MutationTaster prediction, "A" ("disease_causing_automatic"),
		"D" ("disease_causing"), "N" ("polymorphism") or "P" ("polymorphism_automatic"). The 
		score cutoff between "D" and "N" is 0.5 for MTnew and 0.31733 for the rankscore.

MutationAssessor_pred
61	MutationAssessor_pred: MutationAssessor s functional impact of a variant -
		predicted functional, i.e. high ("H") or medium ("M"), or predicted non-functional,
		i.e. low ("L") or neutral ("N"). The MAori score cutoffs between "H" and "M", 
		"M" and "L", and "L" and "N", are 3.5, 1.935 and 0.8, respectively. The rankscore cutoffs 
		between "H" and "M", "M" and "L", and "L" and "N", are 0.9307, 0.52043 and 0.19675, 
		respectively.

FATHMM_pred
64	FATHMM_pred: If a FATHMMori score is <=-1.5 (or rankscore >=0.81332) the corresponding nsSNV
		is predicted as "D(AMAGING)"; otherwise it is predicted as "T(OLERATED)".
		Multiple predictions separated by ";", corresponding to Ensembl_proteinid.

PROVEAN_pred
67	PROVEAN_pred: If PROVEANori <= -2.5 (rankscore>=0.54382) the corresponding nsSNV is
		predicted as "D(amaging)"; otherwise it is predicted as "N(eutral)". 
		Multiple predictions separated by ";", corresponding to Ensembl_proteinid.

MetaLR_pred
75	MetaLR_pred: Prediction of our MetaLR based ensemble prediction score,"T(olerated)" or
		"D(amaging)". The score cutoff between "D" and "T" is 0.5. The rankscore cutoff between 
		"D" and "T" is 0.81101.

MetaSVM_pred
72	MetaSVM_pred: Prediction of our SVM based ensemble prediction score,"T(olerated)" or
		"D(amaging)". The score cutoff between "D" and "T" is 0. The rankscore cutoff between
		"D" and "T" is 0.82257.

M-CAP_pred
82	M-CAP_pred: Prediction of M-CAP score based on the authors recommendation, "T(olerated)" or
		"D(amaging)". The score cutoff between "D" and "T" is 0.025.

Interpro_domain
640	Interpro_domain: domain or conserved site on which the variant locates. Domain annotations come from Interpro database. The number in the brackets following a specific domain is the count of times Interpro assigns the variant position to that domain, typically coming from different predicting databases. Multiple entries separated by ";".

GTEx_V8_gene
641	GTEx_V8_gene: target gene of the (significant) eQTL SNP

GTEx_V8_tissue
642	GTEx_V8_tissue: tissue type of the expression data with which the eQTL/gene pair is detected















#==============#
# Installation #
#==============#

VEP_CACHE="/home/gonzalo/.vep/"


docker run -t -i -u $(id -u):$(id -g) -v ${VEP_CACHE}:/opt/vep/.vep ensemblorg/ensembl-vep perl INSTALL.pl -a cf -s homo_sapiens_refseq -y GRCh37 --CACHEDIR /opt/vep/.vep
#docker run -t -i -u $(id -u):$(id -g) -v ${VEP_CACHE}:/opt/vep/.vep ensemblorg/ensembl-vep perl INSTALL.pl -a p -s homo_sapiens_refseq -y GRCh37 --CACHEDIR /opt/vep/.vep -g dbscSNV,LoFtool,ExACpLI,dbNSFP,MaxEntScan,CADD,SpliceAI
docker run -t -i -u $(id -u):$(id -g) -v ${VEP_CACHE}:/opt/vep/.vep ensemblorg/ensembl-vep perl INSTALL.pl -a p -s homo_sapiens_refseq -y GRCh37 --CACHEDIR /opt/vep/.vep -g all


# CADD (April 11, 2020)
# CADD is a tool for scoring the deleteriousness of single nucleotide variants as well as insertion/deletions variants in the human genome.
# CADD scores are freely available for all non-commercial applications. If you are planning on using them in a commercial application, please obtain a license.
wget https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh37/whole_genome_SNVs.tsv.gz
wget https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh37/whole_genome_SNVs.tsv.gz.tbi
wget https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh37/InDels.tsv.gz
wget https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh37/InDels.tsv.gz.tbi
mv whole_genome_SNVs.tsv.gz CADD_v.6_hg19_whole_genome_SNVs.tsv.gz
mv whole_genome_SNVs.tsv.gz.tbi CADD_v.6_hg19_whole_genome_SNVs.tsv.gz.tbi
mv InDels.tsv.gz CADD_v.6_hg19_InDels.tsv.gz
mv InDels.tsv.gz.tbi CADD_v.6_hg19_InDels.tsv.gz.tbi

wget https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz
wget https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz.tbi
wget https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.indel.tsv.gz
wget https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.indel.tsv.gz.tbi
mv whole_genome_SNVs.tsv.gz CADD_v.6_hg38_whole_genome_SNVs.tsv.gz
mv whole_genome_SNVs.tsv.gz.tbi CADD_v.6_hg38_whole_genome_SNVs.tsv.gz.tbi
mv gnomad.genomes.r3.0.indel.tsv.gz CADD_v.6_hg38_gnomad.genomes.r3.0.indel.tsv.gz
mv gnomad.genomes.r3.0.indel.tsv.gz.tbi CADD_v.6_hg38_gnomad.genomes.r3.0.indel.tsv.gz.tbi



# dbNSFP (February 18, 2022)
# dbNSFP is a database developed for functional prediction and annotation of all potential non-synonymous single-nucleotide variants (nsSNVs) in the human genome.
# It compiles prediction scores from 38 prediction algorithms (SIFT, SIFT4G, Polyphen2-HDIV, Polyphen2-HVAR, LRT, MutationTaster2, MutationAssessor, FATHMM, MetaSVM, MetaLR, MetaRNN, CADD, CADD_hg19, VEST4, PROVEAN, FATHMM-MKL coding, FATHMM-XF coding, fitCons x 4, LINSIGHT, DANN, GenoCanyon, Eigen, Eigen-PC, M-CAP, REVEL, MutPred, MVP, MPC, PrimateAI, GEOGEN2, BayesDel_addAF, BayesDel_noAF, ClinPred, LIST-S2, ALoFT), 9 conservation scores (PhyloP x 3, phastCons x 3, GERP++, SiPhy and bStatistic) and other related information including allele frequencies observed in the 1000 Genomes Project phase 3 data, UK10K cohorts data, ExAC consortium data, gnomAD data and the NHLBI Exome Sequencing Project ESP6500 data, various gene IDs from different databases, functional descriptions of genes, gene expression and gene interaction information, etc.
wget ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFP4.3a.zip
unzip dbNSFP4.3a.zip
zcat dbNSFP4.3a_variant.chr1.gz | head -n1 > h
zgrep -h -v ^#chr dbNSFP4.3a_variant.chr* | awk '$8 != "." ' | sort -T /mnt/genetica6/tmp_sort -k8,8 -k9,9n - | cat h - | bgzip -c > dbNSFP4.3a_grch37.gz
tabix -s 8 -b 9 -e 9 dbNSFP4.3a_grch37.gz

zgrep -h -v ^#chr dbNSFP4.3a_variant.chr* | sort -T /mnt/genetica6/tmp_sort -k1,1 -k2,2n - | cat h - | bgzip -c > dbNSFP4.3a_grch38.gz
tabix -s 1 -b 2 -e 2 dbNSFP4.3a_grch38.gz

rm dbNSFP4.3a_variant.chr* tryhg* try.vcf search_dbNSFP43a* LICENSE.txt
cut -f1,14,20,26,27,30,31,32,33,36,37,95 dbNSFP4.3_gene.complete.txt > dbNSFP4.3_gene.complete.pvm.txt

"""R
dbNSFPgenepath = "/home/gonzalo/tblab/mnt/genetica7/vep_plugins_gene/dbNSFP4.3_gene.complete.pvm.txt"
dbNSFP_gene = read.delim(dbNSFPgenepath, header = TRUE, stringsAsFactors = F, quote = "")

dbNSFP_gene$MGI_mouse_phenotype_filt = gsub("\\(.*?\\)","", dbNSFP_gene$MGI_mouse_phenotype, perl = T)

dbNSFP_gene = write.table(dbNSFP_gene, dbNSFPgenepath, col.names = TRUE, row.names = F, sep = "\t", quote = F)
"""



# dbscSNV (April 12, 2015)
# dbscSNV includes all potential human SNVs within splicing consensus regions (−3 to +8 at the 5’ splice site and −12 to +2 at the 3’ splice site), i.e. scSNVs, related functional annotations and two ensemble prediction scores for predicting their potential of altering splicing.
wget ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbscSNV1.1.zip
unzip dbscSNV1.1.zip
head -n1 dbscSNV1.1.chr1 > h
cat dbscSNV1.1.chr* | grep -v ^chr | cat h - | bgzip -c > dbscSNV1.1_GRCh37.txt.gz
tabix -s 1 -b 2 -e 2 -c c dbscSNV1.1_GRCh37.txt.gz

cat dbscSNV1.1.chr* | grep -v ^chr | sort -k5,5 -k6,6n | cat h - | awk '$5 != "."' | bgzip -c > dbscSNV1.1_GRCh38.txt.gz
tabix -s 5 -b 6 -e 6 -c c dbscSNV1.1_GRCh38.txt.gz

rm dbscSNV1.1.chr* h


# ExACpLI
# A VEP plugin that adds the probabililty of a gene being loss-of-function intolerant (pLI) to the VEP output.
# The closer pLI is to 1, the more likely the gene is loss-of-function (LoF) intolerant.
# wget ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/functional_gene_constraint/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt
wget https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/105/ExACpLI_values.txt
# mv ExACpLI_values.txt ${VEP_CACHE}/Plugins/


# LoFtool
# LoFtool provides a rank of genic intolerance and consequent susceptibility to disease based on the ratio of Loss-of-function (LoF) to synonymous mutations for each gene
# The lower the LoFtool gene score percentile the most intolerant is the gene to functional variation.
wget https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/105/LoFtool_scores.txt
# mv LoFtool_scores.txt ${VEP_CACHE}/Plugins/


# MaxEntScan
# 
# http://genes.mit.edu/burgelab/maxent/download/
wget http://hollywood.mit.edu/burgelab/maxent/download/fordownload.tar.gz
gunzip fordownload.tar.gz
mv fordownload.tar maxEntScan.tar
tar -xvf maxEntScan.tar


# SpliceAI
# A VEP plugin that retrieves pre-calculated annotations from SpliceAI. 
# SpliceAI is a deep neural network, developed by Illumina, Inc that predicts splice junctions from an arbitrary pre-mRNA transcript sequence.
# Delta score of a variant, defined as the maximum of (DS_AG, DS_AL, DS_DG, DS_DL), ranges from 0 to 1 and can be interpreted as the probability of the variant being splice-altering. The author-suggested cutoffs are:
#   0.2 (high recall)
#   0.5 (recommended)
#   0.8 (high precision)
# The output includes the gene symbol, delta scores (DS) and delta positions (DP) for acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
# Son 422GB !!!!


# DENOVO_DB (V1.6.1)
# De novo germline variants were obtained from denovo-db at
wget https://denovo-db.gs.washington.edu/denovo-db.non-ssc-samples.variants.vcf.gz
gunzip denovo-db.non-ssc-samples.variants.vcf.gz
bgzip denovo-db.non-ssc-samples.variants.vcf 
tabix -p vcf denovo-db.non-ssc-samples.variants.vcf.gz
mv denovo-db.non-ssc-samples.variants.vcf.gz denovo-db.non-ssc-samples.variants.hg19.vcf.gz
mv denovo-db.non-ssc-samples.variants.vcf.gz.tbi denovo-db.non-ssc-samples.variants.hg19.vcf.gz.tbi


# CCRs
# Constrained Coding Regions
wget https://s3.us-east-2.amazonaws.com/ccrs/ccrs/ccrs.autosomes.v2.20180420.bed.gz
wget https://s3.us-east-2.amazonaws.com/ccrs/ccrs/ccrs.xchrom.v2.20180420.bed.gz
cat ccrs.autosomes.v2.20180420.bed.gz ccrs.xchrom.v2.20180420.bed.gz > ccrs.all.v2.20180420.bed.gz 
tabix -p bed ccrs.all.v2.20180420.bed.gz
rm ccrs.autosomes.v2.20180420.bed.gz ccrs.xchrom.v2.20180420.bed.gz


# gnomAD genome noVEP sites
# wget http://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh37/variation_genotype/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz
# wget http://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh37/variation_genotype/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz.tbi


# Kaviar (February 29, 2016 (version 160204-Public))
wget http://s3-us-west-2.amazonaws.com/kaviar-160204-public/Kaviar-160204-Public-hg19.vcf.tar
tar -xvf Kaviar-160204-Public-hg19.vcf.tar
bgzip Kaviar-160204-Public-hg19.vcf
tabix -p vcf Kaviar-160204-Public-hg19.vcf.gz
mv Kaviar-160204-Public/vcfs/Kaviar-160204-Public-hg19.vcf.gz* .
rm -r Kaviar-160204-Public Kaviar-160204-Public-hg19.vcf.tar

wget http://s3-us-west-2.amazonaws.com/kaviar-160204-public/Kaviar-160204-Public-hg38.vcf.tar
tar -xvf Kaviar-160204-Public-hg38.vcf.tar
mv Kaviar-160204-Public/vcfs/Kaviar-160204-Public-hg38.vcf.gz* .
rm -r Kaviar-160204-Public Kaviar-160204-Public-hg38.vcf.tar




# ClinVar
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi
mv clinvar.vcf.gz clinvar.GRCh37.vcf.gz
mv clinvar.vcf.gz.tbi clinvar.GRCh37.vcf.gz.tbi

wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
mv clinvar.vcf.gz clinvar.GRCh38.vcf.gz
mv clinvar.vcf.gz.tbi clinvar.GRCh38.vcf.gz.tbi




# MutScore
wget https://storage.googleapis.com/rivolta_mutscore/mutscore-v1.0-hg19.tsv.gz
vcfFromBed="mutscore-v1.0-hg19.vcf.gz"
bcftools view -h gnomad.exomes.r2.1.1.hg19.vcf.bgz | head -n 1 | bgzip -c > ${vcfFromBed}
bcftools view -h gnomad.exomes.r2.1.1.hg19.vcf.bgz | grep "##contig=" | bgzip -c >> ${vcfFromBed}
echo '##INFO=<ID=Score,Number=A,Type=String,Description="MutScore">' | bgzip -c  >> ${vcfFromBed}
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" | bgzip -c >> ${vcfFromBed}
zcat mutscore-v1.0-hg19.tsv.gz | awk -F"\t" '{print "chr"$1"\t"$2"\t.\t"$3"\t"$4"\t.\t.\tScore="$5}' | bgzip -c >> ${vcfFromBed}
tabix -p vcf ${vcfFromBed}


wget https://storage.googleapis.com/rivolta_mutscore/mutscore-v1.0-hg38.tsv.gz
vcfFromBed="mutscore-v1.0-hg38.vcf.gz"
bcftools view -h gnomad.genomes.v3.1.2.hg38.vcf.bgz | head -n 1 | bgzip -c > ${vcfFromBed}
bcftools view -h gnomad.genomes.v3.1.2.hg38.vcf.bgz | grep "##contig=" | bgzip -c >> ${vcfFromBed}
echo '##INFO=<ID=Score,Number=A,Type=String,Description="MutScore">' | bgzip -c  >> ${vcfFromBed}
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" | bgzip -c >> ${vcfFromBed}
zcat mutscore-v1.0-hg19.tsv.gz | awk -F"\t" '{print "chr"$1"\t"$2"\t.\t"$3"\t"$4"\t.\t.\tScore="$5}' | bgzip -c >> ${vcfFromBed}
tabix -p vcf ${vcfFromBed}





# gnomad
# Keep confidence sites

bcftools annotate -x ^INFO/AC,INFO/AN,INFO/AF,INFO/nhomalt,\
INFO/popmax,INFO/AC_popmax,INFO/AN_popmax,INFO/AF_popmax,INFO/nhomalt_popmax,\
INFO/AC_nfe,INFO/AN_nfe,INFO/AF_nfe,INFO/nhomalt_nfe \
--threads 8 -O z -o gnomad.genomes.r2.1.1.sites.1.annotfilt.vcf.gz gnomad.genomes.r2.1.1.sites.1.vcf.bgz

bcftools view -f 'PASS' --threads 8 -O z -o gnomad.genomes.r2.1.1.sites.1.annotfilt.pass.vcf.gz gnomad.genomes.r2.1.1.sites.1.annotfilt.vcf.gz


##INFO=<ID=AC,Number=A,Type=Integer,Description="Alternate allele count for samples">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in samples">
##INFO=<ID=AF,Number=A,Type=Float,Description="Alternate allele frequency in samples">
##INFO=<ID=nhomalt,Number=A,Type=Integer,Description="Count of homozygous individuals in samples">

##INFO=<ID=popmax,Number=A,Type=String,Description="Population with maximum AF">
##INFO=<ID=AC_popmax,Number=A,Type=Integer,Description="Allele count in the population with the maximum AF">
##INFO=<ID=AN_popmax,Number=A,Type=Integer,Description="Total number of alleles in the population with the maximum AF">
##INFO=<ID=AF_popmax,Number=A,Type=Float,Description="Maximum allele frequency across populations (excluding samples of Ashkenazi, Finnish, and indeterminate ancestry)">
##INFO=<ID=nhomalt_popmax,Number=A,Type=Integer,Description="Count of homozygous individuals in the population with the maximum allele frequency">

##INFO=<ID=AC_nfe,Number=A,Type=Integer,Description="Alternate allele count for samples of Non-Finnish European ancestry">
##INFO=<ID=AN_nfe,Number=1,Type=Integer,Description="Total number of alleles in samples of Non-Finnish European ancestry">
##INFO=<ID=AF_nfe,Number=A,Type=Float,Description="Alternate allele frequency in samples of Non-Finnish European ancestry">
##INFO=<ID=nhomalt_nfe,Number=A,Type=Integer,Description="Count of homozygous individuals in samples of Non-Finnish European ancestry">

# todo en la general, poner la frecuencia y el AC en las otras dos. En la pop max anotar población




#================================#
# Download gnomAD v3.1.2 genomes #
#================================#

# Download
chrs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX)

# Processing (reduce columns)
chrs=(chr1 chr2 chr3 chr4 chr5)
# chrs=(chr6 chr7)
chrs=(chr7 chr8 chr9 chr10 chr11 chr12) 
chrs=(chr13 chr14 chr15 chr16 chr17) 
# chrs=(chr16 chr17) 
chrs=(chr18 chr19 chr20) 
chrs=(chr21 chr22 chrX chrY) 
# chrs=(chr22 chrX chrY)


for chr in "${chrs[@]}"
do
echo "https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.${chr}.vcf.bgz"
wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.${chr}.vcf.bgz
wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.${chr}.vcf.bgz.tbi

md5sum gnomad.genomes.v3.1.2.sites.${chr}.vcf.bgz

bcftools annotate -x ^INFO/AC,INFO/AN,INFO/AF,INFO/nhomalt,\
INFO/popmax,INFO/AC_popmax,INFO/AN_popmax,INFO/AF_popmax,INFO/nhomalt_popmax,\
INFO/AC_nfe,INFO/AN_nfe,INFO/AF_nfe,INFO/nhomalt_nfe \
--threads 8 -O z -o gnomad.genomes.v3.1.2.sites.${chr}.annotfilt.vcf.bgz gnomad.genomes.v3.1.2.sites.${chr}.vcf.bgz
done

65b21b95252786012721de95458a90e3  gnomad.genomes.v3.1.2.sites.chr1.vcf.bgz   OK DELETED
6ad213299e21959d7b91787928d6f8b5  gnomad.genomes.v3.1.2.sites.chr2.vcf.bgz   OK DELETED
13e380299b1cb61d2d18a45a522e8810  gnomad.genomes.v3.1.2.sites.chr3.vcf.bgz-  OK DELETED
90021c65fa7dabfed1a404b713db503f  gnomad.genomes.v3.1.2.sites.chr4.vcf.bgz-  OK DELETED
fc17883bc83787524b7872ba3002a05f  gnomad.genomes.v3.1.2.sites.chr5.vcf.bgz-  OK DELETED
3f1a3910c97e3625fe12962d2ed69a9a  gnomad.genomes.v3.1.2.sites.chr6.vcf.bgz   OK DELETED
98fa44053cf0f907b9ae52ac461bb717  gnomad.genomes.v3.1.2.sites.chr7.vcf.bgz   OK DELETED
5ac0337761d358832aa4adb041d1c3b7  gnomad.genomes.v3.1.2.sites.chr8.vcf.bgz   OK DELETED
1a6fcfc564cd1d8e58db3e330b743303  gnomad.genomes.v3.1.2.sites.chr9.vcf.bgz-  OK DELETED
f84d74e58e9eb568bdfe349e6daa2645  gnomad.genomes.v3.1.2.sites.chr10.vcf.bgz- OK DELETED
843ec265623593200dd64d56a5a045b3  gnomad.genomes.v3.1.2.sites.chr11.vcf.bgz- OK DELETED
31f76395874ec62a7998e9fb3b2850a1  gnomad.genomes.v3.1.2.sites.chr12.vcf.bgz- OK DELETED
78a27d411d3512dbbc6eb4b120391cb8  gnomad.genomes.v3.1.2.sites.chr13.vcf.bgz- OK DELETED
e8273eef9760aac6ad9f799a52082d5d  gnomad.genomes.v3.1.2.sites.chr14.vcf.bgz- OK DELETED
89d6cd19e5c7621f627fab6577cd0e87  gnomad.genomes.v3.1.2.sites.chr15.vcf.bgz- OK DELETED
33ef1541b49d0347d6141609a287bdc7  gnomad.genomes.v3.1.2.sites.chr16.vcf.bgz- OK DELETED
d13b9287a89e11d8119a93bc0f094c77  gnomad.genomes.v3.1.2.sites.chr17.vcf.bgz- OK DELETED
e95595ca7a90edc528a09b209ea9bc4c  gnomad.genomes.v3.1.2.sites.chr18.vcf.bgz- OK DELETED
32308281570563a5bc71a815ac11154f  gnomad.genomes.v3.1.2.sites.chr19.vcf.bgz- OK DELETED
5cb14ed936340875c8ef664dbe6e772d  gnomad.genomes.v3.1.2.sites.chr20.vcf.bgz- OK DELETED
92c4b4f1bd63a7bd64ac8378d88d86e2  gnomad.genomes.v3.1.2.sites.chr21.vcf.bgz- OK DELETED
35f4d25b924ab2e9520c725dd9699ee5  gnomad.genomes.v3.1.2.sites.chr22.vcf.bgz  OK DELETED
040080a18046533728fa60800eedcf4b  gnomad.genomes.v3.1.2.sites.chrX.vcf.bgz
112ed724f8d1cf13b031006def03f55e  gnomad.genomes.v3.1.2.sites.chrY.vcf.bgz





# for chr in "${chrs[@]}"
# do
# md5sum gnomad.genomes.v3.1.2.sites.${chr}.vcf.bgz
# done


# for chr in "${chrs[@]}"
# do
# md5sum gnomad.genomes.v3.1.2.sites.${chr}.vcf.bgz

# bcftools annotate -x ^INFO/AC,INFO/AN,INFO/AF,INFO/nhomalt,\
# INFO/popmax,INFO/AC_popmax,INFO/AN_popmax,INFO/AF_popmax,INFO/nhomalt_popmax,\
# INFO/AC_nfe,INFO/AN_nfe,INFO/AF_nfe,INFO/nhomalt_nfe \
# --threads 8 -O z -o gnomad.genomes.v3.1.2.sites.${chr}.annotfilt.vcf.bgz gnomad.genomes.v3.1.2.sites.${chr}.vcf.bgz

# # bcftools view -f 'PASS' --threads 8 -O z -o gnomad.genomes.v3.1.2.sites.${chr}.annotfilt.pass.vcf.bgz gnomad.genomes.v3.1.2.sites.${chr}.annotfilt.vcf.bgz
# done




# Annotate FILTER column into INFO column

chrs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9) 
chrs=(chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22)
chrs=(chrX chrY)

for chr in "${chrs[@]}"
do 
echo ${chr}
# Header
input_vcf="gnomad.genomes.v3.1.2.sites.${chr}.annotfilt.vcf.bgz"
output_vcf="gnomad.genomes.v3.1.2.sites.${chr}.output.vcf.bgz"
bcftools view -h ${input_vcf} | grep "##" | bgzip -c > ${output_vcf}
echo "##INFO=<ID=filt,Number=.,Type=String,Description=\"Filter\">" | bgzip -c >> ${output_vcf}
bcftools view -h ${input_vcf} | grep "#CHROM" | cut -f1-8 | bgzip -c >> ${output_vcf}

# NEW INFO
bcftools query -f 'filt=%FILTER\n' ${input_vcf} | sed 's/;/_/g' | sed 's/,/_/g' > new_info.${chr}.txt
bcftools view -H ${input_vcf} | cut -f1-8 > old_info.${chr}.txt
paste -d ';' old_info.${chr}.txt new_info.${chr}.txt | bgzip -c >> ${output_vcf}
rm new_info.${chr}.txt old_info.${chr}.txt

tabix -p vcf ${output_vcf}
done


# Concatenation
bcftools concat -O z -o gnomad.genomes.v3.1.2.hg38.vcf.bgz *output.vcf.bgz
tabix -p vcf gnomad.genomes.v3.1.2.hg38.vcf.bgz



#====================================================#
# Liftover gnomAD v3.1.2 genomes from hg 38 to hg 37 #
#====================================================#

# Download chain and genome to permorm the liftover
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
gunzip hg38ToHg19.over.chain.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz
bgzip -c hg19.fa > hg19.fa.gz

# create index file for fasta file
samtools faidx hg19.fa.gz


# Download GATK docker image
docker pull broadinstitute/gatk

# Generate fasta dictionary
docker run -it \
-u $(id -u):$(id -g) \
-v /mnt/genetica7/references:/mnt/references \
broadinstitute/gatk \
./gatk CreateSequenceDictionary -R /mnt/references/hg19.fa.gz


# Perform the liftover
docker run -it \
-u $(id -u):$(id -g) \
-v /mnt/genetica7/vep_plugins_hg38/:/mnt/io_dir \
-v /mnt/genetica7/references:/mnt/references \
broadinstitute/gatk \
./gatk LiftoverVcf \
I=/mnt/io_dir/gnomad.genomes.v3.1.2.hg38.vcf.bgz \
O=/mnt/io_dir/gnomad.genomes.v3.1.2.hg19.vcf.bgz \
CHAIN=/mnt/references/hg38ToHg19.over.chain \
REJECT=/mnt/io_dir/rejected_variants.vcf \
R=/mnt/references/hg19.fa.gz
















#===============================#
# Download gnomAD v2.1.1 exomes #
#===============================#

# Download
wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz
md5sum gnomad.exomes.r2.1.1.sites.vcf.bgz




# Processing (reduce columns)
bcftools annotate -x ^INFO/AC,INFO/AN,INFO/AF,INFO/nhomalt,\
INFO/popmax,INFO/AC_popmax,INFO/AN_popmax,INFO/AF_popmax,INFO/nhomalt_popmax,\
INFO/AC_nfe,INFO/AN_nfe,INFO/AF_nfe,INFO/nhomalt_nfe \
--threads 8 -O z -o gnomad.exomes.r2.1.1.sites.annotfilt.vcf.bgz gnomad.exomes.r2.1.1.sites.vcf.bgz

#bcftools view -f 'PASS' --threads 8 -O z -o gnomad.exomes.r2.1.1.sites.annotfilt.pass.vcf.bgz gnomad.exomes.r2.1.1.sites.annotfilt.vcf.bgz

# Annotate FILTER column into INFO column
	# Header
input_vcf="gnomad.exomes.r2.1.1.sites.annotfilt.vcf.bgz"
output_vcf="gnomad.exomes.r2.1.1.sites.annotfilt.infofilt.vcf"
bcftools view -h ${input_vcf} | grep "##" > ${output_vcf}
echo "##INFO=<ID=filt,Number=.,Type=String,Description=\"Filter\">" >> ${output_vcf}
bcftools view -h ${input_vcf} | grep "#CHROM" | cut -f1-8 >> ${output_vcf}

	# NEW INFO
bcftools query -f 'filt=%FILTER\n' ${input_vcf} | sed 's/;/_/g' | sed 's/,/_/g' > new_info.txt
bcftools view -H ${input_vcf} | cut -f1-8 > old_info.txt
paste -d ';' old_info.txt new_info.txt >> ${output_vcf}
rm new_info.txt old_info.txt


# rm gnomad.exomes.r2.1.1.sites.vcf.bgz

awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' ${output_vcf} | bgzip -c > gnomad.exomes.r2.1.1.hg19.vcf.bgz
tabix -p vcf gnomad.exomes.r2.1.1.hg19.vcf.bgz
rm  ${output_vcf} 




#===================================================#
# Liftover gnomAD v2.1.1 exomes from hg 37 to hg 38 #
#===================================================#

# Download chain and genome to permorm the liftover
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
gunzip hg19ToHg38.over.chain.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
bgzip hg38.fa -c > hg38.fa.gz

# create index file for fasta file
samtools faidx hg38.fa.gz


# Generate fasta dictionary
docker run -it \
-u $(id -u):$(id -g) \
-v /mnt/genetica7/references:/mnt/references \
broadinstitute/gatk \
./gatk CreateSequenceDictionary -R /mnt/references/hg38.fa.gz


# Perform the liftover
docker run -it \
-u $(id -u):$(id -g) \
-v /mnt/genetica7/vep_plugins:/mnt/io_dir \
-v /mnt/genetica7/references:/mnt/references \
broadinstitute/gatk \
./gatk LiftoverVcf \
I=/mnt/io_dir/gnomad.exomes.r2.1.1.hg19.vcf.bgz \
O=/mnt/io_dir/gnomad.exomes.r2.1.1.hg38.vcf.bgz \
CHAIN=/mnt/references/hg19ToHg38.over.chain \
REJECT=/mnt/io_dir/rejected_variants_exomes.vcf \
R=/mnt/references/hg38.fa.gz
















#=========================================#
# Download gnomAD v3.0.1 genomes coverage #
#=========================================#

# Download coverage file (looks like it is 1-based and bed files are 0 based) https://www.biostars.org/p/84686/
wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.0.1/coverage/genomes/gnomad.genomes.r3.0.1.coverage.summary.tsv.bgz

# Filter desired columns
zcat gnomad.genomes.r3.0.1.coverage.summary.tsv.bgz | tail -n +2 | awk -F"[:\t]" '{print $1"\t"$2-1"\t"$2"\t"$4"\t"$10}' | grep -v "chrM" | sort -T /mnt/genetica6/tmp_sort -k1,1 -k2,2n -k3,3n -t$'\t' | bgzip -c > gnomad.genomes.v3.0.1.coverage.hg38.bed.gz
tabix -p bed gnomad.genomes.v3.0.1.coverage.hg38.bed.gz
# mv gnomad.exomes.v3.0.1.coverage.hg38.tsv.gz gnomad.exomes.r2.1.1.coverage.hg38.1based.bed.gz
# zcat gnomad.exomes.r2.1.1.coverage.hg38.1based.bed.gz | awk -F"\t" '{print $1"\t"$2-1"\t"$3"\t"$4"\t"$5}' | grep -v "chrM" | bgzip -c > gnomad.genomes.v3.0.1.coverage.hg38.bed.gz
# tabix -p bed gnomad.genomes.v3.0.1.coverage.hg38.bed.gz


# Get nucleotides
bedtools getfasta -fi /mnt/genetica7/references/hg38.fa -bed gnomad.genomes.v3.0.1.coverage.hg38.bed.gz -bedOut | bgzip -c > gnomad.genomes.v3.0.1.coverage.hg38.positions.bed.gz

# Transform bed file to vcf
vcfFromBed="gnomad.genomes.v3.0.1.coverage.hg38.positions.vcf.gz"
bcftools view -h gnomad.genomes.v3.1.2.hg38.vcf.bgz | head -n 1 | bgzip -c > ${vcfFromBed}
bcftools view -h gnomad.genomes.v3.1.2.hg38.vcf.bgz | grep "##contig=" | bgzip -c >> ${vcfFromBed}
echo '##INFO=<ID=median,Number=A,Type=String,Description="Median depth">' | bgzip -c  >> ${vcfFromBed}
echo '##INFO=<ID=perc_20x,Number=A,Type=String,Description="% of samples over 20x coverage">' | bgzip -c >> ${vcfFromBed}
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"| bgzip -c  >> ${vcfFromBed}
zcat gnomad.genomes.v3.0.1.coverage.hg38.positions.bed.gz | awk -F"\t" '{print $1 "\t" $3 "\t.\t" $6 "\t.\t.\t.\tmedian=" $4";perc_20x="$5}' | bgzip -c >> ${vcfFromBed}
tabix -p vcf ${vcfFromBed}

# Perform the liftover
docker run -it \
-u $(id -u):$(id -g) \
-v /mnt/genetica6/gonzalo:/mnt/io_dir \
-v /mnt/genetica7/references:/mnt/references \
broadinstitute/gatk \
./gatk LiftoverVcf \
I=/mnt/io_dir/gnomad.genomes.v3.0.1.coverage.hg38.positions.vcf.gz \
O=/mnt/io_dir/gnomad.genomes.v3.0.1.coverage.hg37.positions.vcf.gz \
CHAIN=/mnt/references/hg38ToHg19.over.chain \
REJECT=/mnt/io_dir/rejected_variants_genomes.coverage.vcf \
R=/mnt/references/hg19.fa

# rm gnomad.genomes.v3.0.1.coverage.hg38.bed.gz gnomad.genomes.v3.0.1.coverage.hg38.positions.bed.gz









#========================================#
# Download gnomAD v2.1.1 exomes coverage #
#========================================#

# Download coverage file
wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/2.1/coverage/exomes/gnomad.exomes.coverage.summary.tsv.bgz

# Filter desired columns
zcat gnomad.exomes.coverage.summary.tsv.bgz | tail -n +2 | awk -F"\t" '{print chr$1"\t"$2-1"\t"$2"\t"$4"\t"$9}' | grep -v "chrM" | sort -k1,1 -k2,2n -k3,3n -t$'\t' | bgzip -c > gnomad.exomes.r2.1.1.coverage.hg37.bed.gz
tabix -p bed gnomad.exomes.r2.1.1.coverage.hg37.bed.gz
# mv gnomad.exomes.r2.1.1.coverage.hg37.tsv.gz gnomad.exomes.r2.1.1.coverage.hg37.1based.bed.gz
# zcat gnomad.exomes.r2.1.1.coverage.hg37.1based.bed.gz | awk -F"\t" '{print "chr"$1"\t"$2-1"\t"$3"\t"$4"\t"$5}' | grep -v "chrM" | bgzip -c > gnomad.exomes.r2.1.1.coverage.hg37.bed.gz
# tabix -p bed gnomad.exomes.r2.1.1.coverage.hg37.bed.gz


# Get nucleotides
bedtools getfasta -fi /mnt/genetica7/references/hg19.fa -bed gnomad.exomes.r2.1.1.coverage.hg37.bed.gz -bedOut | bgzip -c > gnomad.exomes.r2.1.1.coverage.hg37.positions.bed.gz

# Transform bed file to vcf
vcfFromBed="gnomad.exomes.r2.1.1.coverage.hg37.positions.vcf.gz"
bcftools view -h gnomad.exomes.r2.1.1.hg19.vcf.bgz | head -n 1 | bgzip -c > ${vcfFromBed}
bcftools view -h gnomad.exomes.r2.1.1.hg19.vcf.bgz | grep "##contig=" | bgzip -c >> ${vcfFromBed}
echo '##INFO=<ID=median,Number=A,Type=String,Description="Median depth">' | bgzip -c  >> ${vcfFromBed}
echo '##INFO=<ID=perc_20x,Number=A,Type=String,Description="% of samples over 20x coverage">' | bgzip -c >> ${vcfFromBed}
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" | bgzip -c >> ${vcfFromBed}
zcat gnomad.exomes.r2.1.1.coverage.hg37.positions.bed.gz | awk -F"\t" '{print $1 "\t" $3 "\t.\t" $6 "\t.\t.\t.\tmedian=" $4";perc_20x="$5}' | bgzip -c >> ${vcfFromBed}
tabix -p vcf ${vcfFromBed}


# Perform the liftover
docker run -it \
-u $(id -u):$(id -g) \
-v /mnt/genetica7/vep_plugins:/mnt/io_dir \
-v /mnt/genetica7/references:/mnt/references \
broadinstitute/gatk \
./gatk LiftoverVcf \
I=/mnt/io_dir/gnomad.exomes.r2.1.1.coverage.hg37.positions.vcf.gz \
O=/mnt/io_dir/gnomad.exomes.r2.1.1.coverage.hg38.positions.vcf.gz \
CHAIN=/mnt/references/hg19ToHg38.over.chain \
REJECT=/mnt/io_dir/rejected_variants_exomes.coverage.vcf \
R=/mnt/references/hg38.fa

rm gnomad.exomes.r2.1.1.coverage.hg37.positions.bed.gz gnomad.exomes.r2.1.1.coverage.hg37.bed.gz











#====================================#
# LiftOver denovo variants and  ccrs #
#====================================#

# Perform the liftover to denovo variants
bcftools view -h denovo-db.non-ssc-samples.variants.hg19.vcf.gz | bgzip -c > denovo-db.non-ssc-samples.variants.hg19.chr.vcf.gz
bcftools view -H denovo-db.non-ssc-samples.variants.hg19.vcf.gz | awk '{print "chr"$0}' | bgzip -c >> denovo-db.non-ssc-samples.variants.hg19.chr.vcf.gz
tabix -p vcf denovo-db.non-ssc-samples.variants.hg19.chr.vcf.gz
docker run -it \
-u $(id -u):$(id -g) \
-v /mnt/genetica7/vep_plugins:/mnt/io_dir \
-v /mnt/genetica7/references:/mnt/references \
broadinstitute/gatk \
./gatk LiftoverVcf \
I=/mnt/io_dir/denovo-db.non-ssc-samples.variants.hg19.chr.vcf.gz \
O=/mnt/io_dir/denovo-db.non-ssc-samples.variants.hg38.chr.vcf.gz \
CHAIN=/mnt/references/hg19ToHg38.over.chain \
REJECT=/mnt/io_dir/rejected_variants_denovo-db.non-ssc-samples.vcf \
R=/mnt/references/hg38.fa


# Perform the liftover
awk -F"\t" '{print "chr"$1 "\t" $2 "\t" $3 "\t" $4}' ccrs.all.v2.20180420.bed | grep -v "chr#chrom" > ccrs.all.v2.20180420.hg19.bed
# Subir ccrs.all.v2.20180420.hg19.bed a https://genome.ucsc.edu/cgi-bin/hgLiftOver
bgzip ccrs.all.v2.20180420.hg19.bed
tabix -p bed ccrs.all.v2.20180420.hg19.bed.gz
bgzip ccrs.all.v2.20180420.hg38.bed
zcat ccrs.all.v2.20180420.hg38.bed.gz | sort -k1,1 -k2,2n -k3,3n -t$'\t' | bgzip -c > ccrs.all.v2.20180420.hg38.bed.gz_tmp
mv ccrs.all.v2.20180420.hg38.bed.gz_tmp ccrs.all.v2.20180420.hg38.bed.gz
tabix -p bed ccrs.all.v2.20180420.hg38.bed.gz




##INFO=<ID=AC,Number=A,Type=Integer,Description="Alternate allele count for samples">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in samples">
##INFO=<ID=AF,Number=A,Type=Float,Description="Alternate allele frequency in samples">
##INFO=<ID=nhomalt,Number=A,Type=Integer,Description="Count of homozygous individuals in samples">


#==========================#
# CSVS (spanish frequency) #
#==========================#

# Transform bed file to vcf

for file in /mnt/genetica7/vep_plugins/CSVS_frequencies/all* /mnt/genetica7/vep_plugins/CSVS_frequencies/Healthy.csv
do
echo ${file}
group_name="$(basename ${file} .csv)"
echo ${group_name}

csv_file="${file}"
vcf_hg19="/mnt/genetica7/vep_plugins/CSVS_vcf/${group_name}.hg19.vcf.gz"
vcf_hg38="/mnt/genetica7/vep_plugins/CSVS_vcf/${group_name}.hg38.vcf.gz"

vcfFromBed="${vcf_hg19}_tmp"
bcftools view -h gnomad.exomes.r2.1.1.hg19.vcf.bgz | head -n 1 | bgzip -c > ${vcfFromBed}
bcftools view -h gnomad.exomes.r2.1.1.hg19.vcf.bgz | grep "##contig=" | bgzip -c >> ${vcfFromBed}
echo '##INFO=<ID=AC,Number=A,Type=Integer,Description="Alternate allele count for samples">' | bgzip -c  >> ${vcfFromBed}
echo '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in samples">' | bgzip -c >> ${vcfFromBed}
echo '##INFO=<ID=AF,Number=A,Type=Float,Description="Alternate allele frequency in samples">' | bgzip -c >> ${vcfFromBed}
echo '##INFO=<ID=nhomalt,Number=A,Type=Integer,Description="Count of homozygous individuals in samples">' | bgzip -c >> ${vcfFromBed}
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" | bgzip -c >> ${vcfFromBed}
tail ${csv_file} -n +2 | grep -v "MT" | awk -F"\t" '$4 !~ /,/ {print "chr"$1 "\t" $2 "\t" $5 "\t" $3 "\t" $4 "\t.\t.\t" "AC="$7+($8*2) ";AN="($6+$7+$8)*2 ";AF="$11 ";nhomalt="$8}' | bgzip -c >> ${vcfFromBed}
bcftools sort -O z -o ${vcf_hg19}_sorted ${vcfFromBed}

bcftools norm -m- -c ws -O z -o ${vcf_hg19} -f /mnt/genetica7/references/hg19.fa ${vcf_hg19}_sorted
tabix -p vcf ${vcf_hg19}

rm ${vcf_hg19}_sorted ${vcfFromBed}

# Perform the liftover
echo "Perform the liftover"
docker run -it \
-u $(id -u):$(id -g) \
-v /mnt/genetica7/vep_plugins:/mnt/io_dir \
-v /mnt/genetica7/references:/mnt/references \
broadinstitute/gatk \
./gatk LiftoverVcf \
I=/mnt/io_dir/CSVS_vcf/${group_name}.hg19.vcf.gz \
O=/mnt/io_dir/CSVS_vcf/${group_name}.hg38.vcf.gz \
CHAIN=/mnt/references/hg19ToHg38.over.chain \
REJECT=/mnt/io_dir/rejected_variants_${group_name}.vcf \
R=/mnt/references/hg38.fa

done




#==================#
# LiftOver MAF_FJD #
#==================#
bcftools view -h MAFdb_AN20_latest.vcf.gz | sed 's/\tFORMAT//' | bgzip -c > MAFdb_AN20_latest.hg37.vcf.gz
bcftools view -H MAFdb_AN20_latest.vcf.gz | bgzip -c >> MAFdb_AN20_latest.hg37.vcf.gz
tabix -p vcf MAFdb_AN20_latest.hg37.vcf.gz

docker run -it \
-u $(id -u):$(id -g) \
-v /mnt/genetica7/vep_plugins:/mnt/io_dir \
-v /mnt/genetica7/references:/mnt/references \
broadinstitute/gatk \
./gatk LiftoverVcf \
I=/mnt/io_dir/MAFdb_AN20_latest.hg37.vcf.gz \
O=/mnt/io_dir/MAFdb_AN20_latest.hg38.vcf.gz \
CHAIN=/mnt/references/hg19ToHg38.over.chain \
REJECT=/mnt/io_dir/rejected_variants_MAFdb_AN20.vcf \
R=/mnt/references/hg38.fa



65b21b95252786012721de95458a90e3  gnomad.genomes.v3.1.2.sites.chr1.vcf.bgz
6ad213299e21959d7b91787928d6f8b5  gnomad.genomes.v3.1.2.sites.chr2.vcf.bgz
13e380299b1cb61d2d18a45a522e8810  gnomad.genomes.v3.1.2.sites.chr3.vcf.bgz
90021c65fa7dabfed1a404b713db503f  gnomad.genomes.v3.1.2.sites.chr4.vcf.bgz
fc17883bc83787524b7872ba3002a05f  gnomad.genomes.v3.1.2.sites.chr5.vcf.bgz
3f1a3910c97e3625fe12962d2ed69a9a  gnomad.genomes.v3.1.2.sites.chr6.vcf.bgz
98fa44053cf0f907b9ae52ac461bb717  gnomad.genomes.v3.1.2.sites.chr7.vcf.bgz
5ac0337761d358832aa4adb041d1c3b7  gnomad.genomes.v3.1.2.sites.chr8.vcf.bgz
1a6fcfc564cd1d8e58db3e330b743303  gnomad.genomes.v3.1.2.sites.chr9.vcf.bgz
f84d74e58e9eb568bdfe349e6daa2645  gnomad.genomes.v3.1.2.sites.chr10.vcf.bgz
843ec265623593200dd64d56a5a045b3  gnomad.genomes.v3.1.2.sites.chr11.vcf.bgz
31f76395874ec62a7998e9fb3b2850a1  gnomad.genomes.v3.1.2.sites.chr12.vcf.bgz
78a27d411d3512dbbc6eb4b120391cb8  gnomad.genomes.v3.1.2.sites.chr13.vcf.bgz
e8273eef9760aac6ad9f799a52082d5d  gnomad.genomes.v3.1.2.sites.chr14.vcf.bgz
89d6cd19e5c7621f627fab6577cd0e87  gnomad.genomes.v3.1.2.sites.chr15.vcf.bgz
33ef1541b49d0347d6141609a287bdc7  gnomad.genomes.v3.1.2.sites.chr16.vcf.bgz
d13b9287a89e11d8119a93bc0f094c77  gnomad.genomes.v3.1.2.sites.chr17.vcf.bgz
e95595ca7a90edc528a09b209ea9bc4c  gnomad.genomes.v3.1.2.sites.chr18.vcf.bgz
32308281570563a5bc71a815ac11154f  gnomad.genomes.v3.1.2.sites.chr19.vcf.bgz
5cb14ed936340875c8ef664dbe6e772d  gnomad.genomes.v3.1.2.sites.chr20.vcf.bgz
92c4b4f1bd63a7bd64ac8378d88d86e2  gnomad.genomes.v3.1.2.sites.chr21.vcf.bgz
35f4d25b924ab2e9520c725dd9699ee5  gnomad.genomes.v3.1.2.sites.chr22.vcf.bgz
040080a18046533728fa60800eedcf4b  gnomad.genomes.v3.1.2.sites.chrX.vcf.bgz




-bash-4.2$ ll /home/proyectos/bioinfo/references/VEPdbs -h
total 257G
-rw-rwxr-- 1 lodela  bioinfo 3,1G may 31  2019 af-only-gnomad.raw.sites.hg19.vcf.gz
-rw-rwxr-- 1 lodela  bioinfo 2,4M may 31  2019 af-only-gnomad.raw.sites.hg19.vcf.gz.tbi
-rw-rwxr-- 1 lodela  bioinfo  69M may 31  2019 af-only-gnomad.raw.sites.hg19.vcf.idx
-rwxrwxrwx 1 lodela  bioinfo  603 abr 25  2019 dbConc.sh
-rw-rwxr-- 1 lodela  bioinfo  20G may  8  2019 dbNSFP3.5a_hg19.gz
-rw-rwxr-- 1 lodela  bioinfo 838K may  8  2019 dbNSFP3.5a_hg19.gz.tbi
-rw-rwxr-- 1 lodela  bioinfo  76K ago  6  2017 dbNSFP3.5a.readme.txt
-rw-rwxr-- 1 lodela  bioinfo  12G ago  4  2017 dbNSFP3.5a_variant.chr1
-rw-rwxr-- 1 lodela  bioinfo 4,5G ago  4  2017 dbNSFP3.5a_variant.chr10
-rw-rwxr-- 1 lodela  bioinfo 6,8G ago  4  2017 dbNSFP3.5a_variant.chr11
-rw-rwxr-- 1 lodela  bioinfo 6,2G ago  4  2017 dbNSFP3.5a_variant.chr12
-rw-rwxr-- 1 lodela  bioinfo 2,1G ago  4  2017 dbNSFP3.5a_variant.chr13
-rw-rwxr-- 1 lodela  bioinfo 3,7G ago  4  2017 dbNSFP3.5a_variant.chr14
-rw-rwxr-- 1 lodela  bioinfo 4,0G ago  4  2017 dbNSFP3.5a_variant.chr15
-rw-rwxr-- 1 lodela  bioinfo 4,9G ago  4  2017 dbNSFP3.5a_variant.chr16
-rw-rwxr-- 1 lodela  bioinfo 6,6G ago  4  2017 dbNSFP3.5a_variant.chr17
-rw-rwxr-- 1 lodela  bioinfo 1,9G ago  4  2017 dbNSFP3.5a_variant.chr18
-rw-rwxr-- 1 lodela  bioinfo 7,3G ago  4  2017 dbNSFP3.5a_variant.chr19
-rw-rwxr-- 1 lodela  bioinfo 8,7G ago  4  2017 dbNSFP3.5a_variant.chr2
-rw-rwxr-- 1 lodela  bioinfo 2,8G ago  4  2017 dbNSFP3.5a_variant.chr20
-rw-rwxr-- 1 lodela  bioinfo 1,2G ago  4  2017 dbNSFP3.5a_variant.chr21
-rw-rwxr-- 1 lodela  bioinfo 2,4G ago  4  2017 dbNSFP3.5a_variant.chr22
-rw-rwxr-- 1 lodela  bioinfo 6,8G ago  4  2017 dbNSFP3.5a_variant.chr3
-rw-rwxr-- 1 lodela  bioinfo 4,6G ago  4  2017 dbNSFP3.5a_variant.chr4
-rw-rwxr-- 1 lodela  bioinfo 5,3G ago  4  2017 dbNSFP3.5a_variant.chr5
-rw-rwxr-- 1 lodela  bioinfo 5,9G ago  4  2017 dbNSFP3.5a_variant.chr6
-rw-rwxr-- 1 lodela  bioinfo 5,4G ago  4  2017 dbNSFP3.5a_variant.chr7
-rw-rwxr-- 1 lodela  bioinfo 3,9G ago  4  2017 dbNSFP3.5a_variant.chr8
-rw-rwxr-- 1 lodela  bioinfo 4,6G ago  4  2017 dbNSFP3.5a_variant.chr9
-rw-rwxr-- 1 lodela  bioinfo  17M ago  4  2017 dbNSFP3.5a_variant.chrM
-rw-rwxr-- 1 lodela  bioinfo 4,0G ago  4  2017 dbNSFP3.5a_variant.chrX
-rw-rwxr-- 1 lodela  bioinfo 184M ago  4  2017 dbNSFP3.5a_variant.chrY
-rw-rwxr-- 1 lodela  bioinfo  86M ago  6  2017 dbNSFP3.5_gene
-rw-rwxr-- 1 lodela  bioinfo 115M ago  5  2017 dbNSFP3.5_gene.complete
-rw-rwxr-- 1 lodela  bioinfo  17G may  7  2019 dbNSFPv3.5a.zip
-rwxrwxr-x 1 lodela  bioinfo 350M abr 23  2019 dbscSNV1.1_GRCh37.txt.gz
-rwxrwxr-x 1 lodela  bioinfo 677K abr 23  2019 dbscSNV1.1_GRCh37.txt.gz.tbi
-rwxrwxr-x 1 lodela  bioinfo 863K abr 23  2019 ExAC.r0.3.1.sites.vep.vcf.gz.tbi
-rwxrwxr-x 1 rromero bioinfo  26G sep 25  2019 gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz
-rwxr-xr-x 1 rromero bioinfo 2,7M jun 12  2019 gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz.tbi
-rw-rwxr-- 1 lodela  bioinfo 4,4K may  8  2019 h
-rwxr-xr-x 1 rromero bioinfo 589M jun 13  2019 InDels.tsv.gz
-rwxr-xr-x 1 rromero bioinfo 1,7M jun 14  2019 InDels.tsv.gz.tbi
drwxr-xr-x 3 rromero bioinfo 4,0K jun  7  2019 Kaviar-160204-Public
-rw-rwxr-- 1 lodela  bioinfo 3,2K may 31  2013 LICENSE.txt
drwxr-xr-x 3 rromero bioinfo 4,0K jun 12  2019 maxEntScan
-rw-rwxr-- 1 lodela  bioinfo  27K ago  6  2017 search_dbNSFP35a.class
-rw-rwxr-- 1 lodela  bioinfo  40K ago  6  2017 search_dbNSFP35a.java
-rw-rwxr-- 1 lodela  bioinfo 152K ago  6  2017 search_dbNSFP35a.readme.pdf
-rw-rwxr-- 1 lodela  bioinfo  254 mar 20  2016 tryhg18.in
-rw-rwxr-- 1 lodela  bioinfo  256 abr 12  2015 tryhg19.in
-rw-rwxr-- 1 lodela  bioinfo  254 abr 12  2015 tryhg38.in
-rw-rwxr-- 1 lodela  bioinfo 5,9K nov 30  2016 try.vcf
-rwxr-xr-x 1 rromero bioinfo  78G jun 13  2019 whole_genome_SNVs.tsv.gz
-rwxr-xr-x 1 rromero bioinfo 2,6M jun 13  2019 whole_genome_SNVs.tsv.gz.tbi




-bash-4.2$ ll /usr/local/bioinfo/vep/ensembl-vep-release-103/Plugins -h
total 1,3M
-rw-r--r-- 1 root noruser 7,8K feb  1  2021 AncestralAllele.pm
-rw-r--r-- 1 root noruser 3,5K feb  1  2021 Blosum62.pm
-rw-r--r-- 1 root noruser 3,6K feb  1  2021 CADD.pm
-rw-r--r-- 1 root noruser 4,8K feb  1  2021 Carol.pm
-rw-r--r-- 1 root noruser 1,4K feb  1  2021 CCDSFilter.pm
-rw-r--r-- 1 root noruser 8,3K feb  1  2021 Condel.pm
drwxr-xr-x 3 root noruser 4,0K mar  9  2021 config
-rw-r--r-- 1 root noruser 9,5K feb  1  2021 Conservation.pm
-rw-r--r-- 1 root noruser 2,1K feb  1  2021 CONTRIBUTING.md
-rw-r--r-- 1 root noruser 3,4K feb  1  2021 CSN.pm
-rw-r--r-- 1 root noruser 3,9K feb  1  2021 DAS.pm
-rwxr-xr-x 1 root noruser  15K feb  1  2021 dbNSFP.pm
-rw-r--r-- 1 root noruser   80 feb  1  2021 dbNSFP_replacement_logic
-rw-r--r-- 1 root noruser 6,2K feb  1  2021 dbscSNV.pm
-rw-r--r-- 1 root noruser 8,9K feb  1  2021 DisGeNET.pm
-rwxr-xr-x 1 root noruser 4,8K feb  1  2021 Downstream.pm
-rw-r--r-- 1 root noruser  13K feb  1  2021 Draw.pm
-rw-r--r-- 1 root noruser 3,4K feb  1  2021 ExACpLI.pm
-rw-r--r-- 1 root noruser 459K feb  1  2021 ExACpLI_values.txt
-rw-r--r-- 1 root noruser 9,3K feb  1  2021 ExAC.pm
-rw-r--r-- 1 root noruser 2,7K feb  1  2021 FATHMM_MKL.pm
-rw-r--r-- 1 root noruser 4,3K feb  1  2021 FATHMM.pm
-rw-r--r-- 1 root noruser 2,6K feb  1  2021 FlagLRG.pm
-rw-r--r-- 1 root noruser 3,9K feb  1  2021 FunMotifs.pm
-rw-r--r-- 1 root noruser  57K feb  1  2021 G2P.pm
-rw-r--r-- 1 root noruser 9,5K feb  1  2021 GeneSplicer.pm
-rw-r--r-- 1 root noruser 6,5K feb  1  2021 gnomADc.pm
-rw-r--r-- 1 root noruser 2,4K feb  1  2021 GO.pm
-rw-r--r-- 1 root noruser 3,7K feb  1  2021 GXA.pm
-rw-r--r-- 1 root noruser 1,7K feb  1  2021 HGVSIntronOffset.pm
-rw-r--r-- 1 root noruser 2,1K feb  1  2021 HGVSReferenceBase.pm
-rw-r--r-- 1 root noruser  11K feb  1  2021 LD.pm
-rw-r--r-- 1 root noruser  12K feb  1  2021 LICENSE
-rwxr-xr-x 1 root noruser  11K feb  1  2021 LocalID.pm
-rw-r--r-- 1 root noruser 3,1K feb  1  2021 LoFtool.pm
-rw-r--r-- 1 root noruser 233K feb  1  2021 LoFtool_scores.txt
-rw-r--r-- 1 root noruser 3,3K feb  1  2021 LOVD.pm
-rw-r--r-- 1 root noruser  14K feb  1  2021 Mastermind.pm
-rw-r--r-- 1 root noruser  27K feb  1  2021 MaxEntScan.pm
-rw-r--r-- 1 root noruser 2,9K feb  1  2021 miRNA.pm
-rw-r--r-- 1 root noruser 3,2K feb  1  2021 MPC.pm
-rw-r--r-- 1 root noruser 4,0K feb  1  2021 MTR.pm
-rw-r--r-- 1 root noruser 4,0K feb  1  2021 NearestExonJB.pm
-rwxr-xr-x 1 root noruser 3,4K feb  1  2021 NearestGene.pm
-rw-r--r-- 1 root noruser 3,1K feb  1  2021 neXtProt_headers.txt
-rw-r--r-- 1 root noruser  13K feb  1  2021 neXtProt.pm
-rw-r--r-- 1 root noruser 1,6K feb  1  2021 NonSynonymousFilter.pm
-rw-r--r-- 1 root noruser  13K feb  1  2021 Phenotypes.pm
-rw-r--r-- 1 root noruser  44K feb  1  2021 plugin_config.txt
-rwxr-xr-x 1 root noruser 7,3K feb  1  2021 PolyPhen_SIFT.pm
-rw-r--r-- 1 root noruser 3,1K feb  1  2021 PON_P2.pm
-rw-r--r-- 1 root noruser 7,5K feb  1  2021 PostGAP.pm
-rw-r--r-- 1 root noruser 4,7K feb  1  2021 ProteinSeqs.pm
-rw-r--r-- 1 root noruser 3,0K feb  1  2021 RankFilter.pm
-rw-r--r-- 1 root noruser 1,3K feb  1  2021 README
-rw-r--r-- 1 root noruser 1,1K feb  1  2021 README.md
-rw-r--r-- 1 root noruser 6,9K feb  1  2021 ReferenceQuality.pm
-rw-r--r-- 1 root noruser 5,7K feb  1  2021 RefSeqHGVS.pm
-rw-r--r-- 1 root noruser 3,2K feb  1  2021 REVEL.pm
-rw-r--r-- 1 root noruser 3,2K feb  1  2021 SameCodon.pm
-rw-r--r-- 1 root noruser 8,9K feb  1  2021 satMutMPRA.pm
-rw-r--r-- 1 root noruser 2,5K feb  1  2021 SingleLetterAA.pm
-rw-r--r-- 1 root noruser 9,1K feb  1  2021 SpliceAI.pm
-rw-r--r-- 1 root noruser 5,0K feb  1  2021 SpliceRegion.pm
-rw-r--r-- 1 root noruser  10K feb  1  2021 StructuralVariantOverlap.pm
-rw-r--r-- 1 root noruser 6,5K feb  1  2021 SubsetVCF.pm
-rw-r--r-- 1 root noruser 1,8K feb  1  2021 TSSDistance.pm




(base) gonzalo@tblabserver:~$ docker run -t -i -u $(id -u):$(id -g) -v ${VEP_CACHE}:/opt/vep/.vep ensemblorg/ensembl-vep perl INSTALL.pl -a p -s homo_sapiens_refseq -y GRCh37 --CACHEDIR /opt/vep/.vep -g dbscSNV,LoFtool,ExACpLI,dbNSFP,MaxEntScan,CADD,SpliceAI
WARNING: Unable to carry out version check for 'ensembl-vep'

 - installing "dbscSNV"
 - This plugin requires installation
 - This plugin requires data
 - See /opt/vep/.vep/Plugins/dbscSNV.pm for details
 - OK

 - installing "LoFtool"
 - This plugin requires data
 - See /opt/vep/.vep/Plugins/LoFtool.pm for details
 - OK

 - installing "ExACpLI"
 - This plugin requires data
 - See /opt/vep/.vep/Plugins/ExACpLI.pm for details
 - OK

 - installing "dbNSFP"
 - This plugin requires installation
 - This plugin requires data
 - See /opt/vep/.vep/Plugins/dbNSFP.pm for details
 - OK

 - installing "MaxEntScan"
 - This plugin requires installation
 - See /opt/vep/.vep/Plugins/MaxEntScan.pm for details
 - OK

 - installing "CADD"
 - This plugin requires data
 - See /opt/vep/.vep/Plugins/CADD.pm for details
 - OK

 - installing "SpliceAI"
 - add "--plugin SpliceAI" to your VEP command to use this plugin
 - OK

NB: One or more plugins that you have installed will not work without installation or downloading data; see logs above

All done




















docker pull biocontainers/bwa:v0.7.17_cv1
docker pull broadinstitute/gatk
docker pull ensemblorg/ensembl-vep

docker pull gonumo/automap
docker pull gonumo/basespace






docker run -u $(id -u):$(id -g) -it \
-v /mnt/genetica7/references/:/mnt/genetica7/references/ \
bwa_samtools bwa index /mnt/genetica7/references/hg19.fa.gz



docker run -u $(id -u):$(id -g) -it \
-v /mnt/genetica7/references/:/mnt/genetica7/references/ \
broadinstitute/gatk ./gatk \
IndexFeatureFile -I /mnt/genetica7/references/1000G_phase1.indels.hg19.sites.vcf





docker run -u $(id -u):$(id -g) -it \
-v /mnt/genetica7/references/:/mnt/genetica7/references/ \
broadinstitute/gatk ./gatk \
IndexFeatureFile -I /mnt/genetica7/references/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf






docker run -u $(id -u):$(id -g) -it \
-v /mnt/genetica7/references/:/mnt/genetica7/references/ \
broadinstitute/gatk ./gatk \
IndexFeatureFile -I /mnt/genetica7/references/dbsnp_138.hg19.vcf




# Prepare idex for dragen
# https://gatk.broadinstitute.org/hc/en-us/articles/4407897446939--How-to-Run-germline-single-sample-short-variant-discovery-in-DRAGEN-mode
docker run -u $(id -u):$(id -g) -it \
-v /mnt/genetica7/references/:/mnt/genetica7/references/ \
broadinstitute/gatk ./gatk \
ComposeSTRTableFile \
-R /mnt/genetica7/references/hg19.fa.gz \
-O /mnt/genetica7/references/hg19.str.zip


docker run -u $(id -u):$(id -g) -it \
-v /mnt/genetica7/references/:/mnt/genetica7/references/ \
broadinstitute/gatk ./gatk \
ComposeSTRTableFile \
-R /mnt/genetica7/references/hg38.fa.gz \
-O /mnt/genetica7/references/hg38.str.zip

















################################
# Reference genome preparation #
################################


# Download GATK docker image
docker pull broadinstitute/gatk
docker pull biocontainers/bwa:v0.7.17_cv1

printf "chr1\nchr2\nchr3\nchr4\nchr5\nchr6\nchr7\nchr8\nchr9\nchr10\nchr11\nchr12\nchr13\nchr14\nchr15\nchr16\nchr17\nchr18\nchr19\nchr20\nchr21\nchr22\nchrX\nchrY\nchrM\n" > canonical_chr.txt

workingdir="/mnt/genetica7/references_masked"

########
# hg19 #
########

# Download reference genome
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz

# Select primary assembly chromosomes
samtools faidx -r canonical_chr.txt hg19.fa > new_hg19.fa
mv new_hg19.fa hg19.fa
bgzip hg19.fa
samtools faidx hg19.fa.gz


# Generate fasta dictionary and gzi
docker run -it \
-u $(id -u):$(id -g) \
-v ${workingdir}:${workingdir} \
-w ${workingdir} \
broadinstitute/gatk \
gatk CreateSequenceDictionary -R hg19.fa.gz


# Prepare idex for dragen
# https://gatk.broadinstitute.org/hc/en-us/articles/4407897446939--How-to-Run-germline-single-sample-short-variant-discovery-in-DRAGEN-mode
docker run -it \
-u $(id -u):$(id -g) \
-v ${workingdir}:${workingdir} \
-w ${workingdir} \
broadinstitute/gatk \
gatk ComposeSTRTableFile \
-R ${workingdir}/hg19.fa.gz \
-O ${workingdir}/hg19.str_table.tsv


# Generate BWA index
docker run -it \
-u $(id -u):$(id -g) \
-v ${workingdir}:${workingdir} \
-w ${workingdir} \
bwa_samtools \
bwa index ${workingdir}/hg19.fa.gz






########
# hg38 #
########

# Download reference genome
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
bgzip hg38.fa -c > hg38.fa.gz


# Select primary assembly chromosomes
samtools faidx -r canonical_chr.txt hg38.fa > new_hg38.fa
mv new_hg38.fa hg38.fa
bgzip hg38.fa
samtools faidx hg38.fa.gz


# Generate fasta dictionary and gzi
docker run -it \
-u $(id -u):$(id -g) \
-v ${workingdir}:${workingdir} \
-w ${workingdir} \
broadinstitute/gatk \
gatk CreateSequenceDictionary -R hg38.fa.gz


# Download gatk legacy bundles for Base Quality Score Recalibration (BQSR)
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.indels.hg19.sites.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.indels.hg19.sites.vcf.idx.gz

wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx.gz

wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.vcf.idx.gz





# Prepare idex for dragen
# https://gatk.broadinstitute.org/hc/en-us/articles/4407897446939--How-to-Run-germline-single-sample-short-variant-discovery-in-DRAGEN-mode
docker run -it \
-u $(id -u):$(id -g) \
-v ${workingdir}:${workingdir} \
-w ${workingdir} \
broadinstitute/gatk \
gatk ComposeSTRTableFile \
-R ${workingdir}/hg38.fa.gz \
-O ${workingdir}/hg38.str_table.tsv


# Generate BWA index
docker run -it \
-u $(id -u):$(id -g) \
-v ${workingdir}:${workingdir} \
-w ${workingdir} \
bwa_samtools \
bwa index ${workingdir}/hg38.fa.gz


# Download gatk legacy bundles for Base Quality Score Recalibration (BQSR)
# Data obtained from: https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi

wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi






# Install singularity
echo 'export PATH=/usr/local/go/bin:$PATH' >> ~/.bashrc && \
source ~/.bashrc



export VERSION=3.11.0 && # adjust this as necessary \
wget https://github.com/sylabs/singularity/releases/download/v${VERSION}/singularity-ce-${VERSION}.tar.gz && \
tar -xzf singularity-ce-${VERSION}.tar.gz && \
cd singularity-ce-${VERSION}


./mconfig && \
make -C builddir && \
sudo make -C builddir install





docker tag manta tblabfjd/manta:1.6.0
docker tag annotsv tblabfjd/annotsv:3.1.1
docker tag automap tblabfjd/bioinfotools:1.0.0
docker tag convading tblabfjd/convading:1.2.1
docker tag bwa_samtools tblabfjd/bwa:0.7.17-r1198-dirty
docker tag basespace tblabfjd/basespace:1.5.1

docker push tblabfjd/manta:1.6.0
docker push tblabfjd/annotsv:3.1.1
docker push tblabfjd/bioinfotools:2.0.0
# docker push tblabfjd/convading:1.2.1
docker push tblabfjd/bwa:0.7.17-r1198-dirty
docker push tblabfjd/basespace:1.5.1



singularity pull gatk_4.4.0.0.sif docker://broadinstitute/gatk:4.4.0.0
# singularity pull gatk_4.3.0.0.sif docker://broadinstitute/gatk:4.3.0.0
singularity pull ensembl-vep_release_105.0.sif docker://ensemblorg/ensembl-vep:release_105.0
singularity pull deepvariant_1.4.0.sif docker://google/deepvariant:1.4.0

singularity pull manta_1.6.0.sif docker://tblabfjd/manta:1.6.0
singularity pull annotsv_3.1.1.sif docker://tblabfjd/annotsv:3.1.1
singularity pull bioinfotools_2.0.0.sif docker://tblabfjd/bioinfotools:2.0.0
# singularity pull convading_1.2.1.sif docker://tblabfjd/convading:1.2.1
singularity pull bwa_0.7.17-r1198-dirty.sif docker://tblabfjd/bwa:0.7.17-r1198-dirty
singularity pull basespace_1.5.1.sif docker://tblabfjd/basespace:1.5.1

