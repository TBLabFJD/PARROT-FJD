#!/bin/sh

###########################
### TASK: SNV annotaton ###
###########################


###############
## ARGUMENTS ##
###############


threads=$1
VCF_IN=$2
softwareFile=$3
interval=$4 # must be a chromosome 



#####################
## MODULES AND DBs ##
#####################




module load vep/release103
VEP="/usr/local/bioinfo/vep/ensembl-vep-release-103/vep"
FILTER_VEP='/usr/local/bioinfo/vep/ensembl-vep-release-103/filter_vep'
VEP_CACHE='/usr/local/bioinfo/vep/ensembl-vep/t/testdata/cache/homo_sapiens/'
VEP_FASTA="/home/proyectos/bioinfo/references/VEPfasta/Homo_sapiens.GRCh37.dna.primary_assembly.fa"
PLUGIN_DIR=/usr/local/bioinfo/vep/ensembl-vep-release-103/Plugins
PLUGIN_DBS="/home/proyectos/bioinfo/references/VEPdbs"
#dbNSFP_DB="${PLUGIN_DBS}/dbNSFP3.5a_hg19.gz"
CCS_DB="/home/proyectos/bioinfo/references/CCS/ccrs.autosomes.v2.20180420.bed.gz"
DENOVO_DB="/home/proyectos/bioinfo/references/denovodb/denovo-db.non-ssc-samples.variants.vcf.gz"
MAF_FJD_COHORT="/home/proyectos/bioinfo/fjd/MAF_FJD_v3.0/db/latest/MAFdb_AN20_latest.vcf.gz"

title="MAPPING"

if [ ! -f $softwareFile ] || ! grep -q $title $softwareFile  ; then
	
	printf "VEP ANNOTATION:\n" >> ${softwareFile}
	module list 2>> ${softwareFile}

fi










###############
## VARIABLES ##
###############


# folder
SNV="$(dirname "${VCF_IN}")" 

# sample name
file="$(basename "${VCF_IN}")" 
name="$(echo $file | cut -d "." -f 1)"


# files 
VCF_FILTERED="${SNV}/${name}${interval}.final.vcf"
VCF_ANNOTATED="${SNV}/${name}${interval}.annotated.vcf"
VCF_ANNOTATED_POPFILT="${SNV}/${name}${interval}.annotated.MAFfiltered.vcf" 





##############
## COMMANDS ##
##############




printf '\n\nFiltering by quality and chromosome...\n'



start=`date +%s`


if [ "$interval" != "" ]; then

	perl $FILTER_VEP \
	-i ${VCF_IN} -o ${VCF_FILTERED} \
	--filter "(FILTER = PASS) and (DP > 10) and (CHROM = $interval)" --force_overwrite 


else

	perl $FILTER_VEP \
	-i ${VCF_IN} -o ${VCF_FILTERED} \
	--filter "(FILTER = PASS) and (DP > 10) and (CHROM in chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY)" --force_overwrite 

fi


if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf '\nVEP FILTERING PASS, DP and chromosome for '${name}' DONE\n' 


else
	printf "\nERROR: PROBLEMS WITH VEP FILTERING FOR QUALITY AND CHR"
	exit 1
fi







printf '\n\nVEP annotation...\n'



perl $VEP \
--cache --offline --hgvs --refseq --dir $VEP_CACHE --dir_plugins $PLUGIN_DIR --v --fork $threads --assembly GRCh37 --fasta $VEP_FASTA --force_overwrite  --no_stats \
--biotype --regulatory --protein --symbol --allele_number --numbers --domains --uniprot --variant_class \
--canonical --vcf \
--sift p --polyphen p --af --max_af \
--format vcf \
--pubmed \
--plugin dbscSNV,$PLUGIN_DBS/dbscSNV1.1_GRCh37.txt.gz \
--plugin LoFtool,$PLUGIN_DIR/LoFtool_scores.txt \
--plugin ExACpLI,$PLUGIN_DIR/ExACpLI_values.txt \
--plugin dbNSFP,${PLUGIN_DBS}/dbNSFP3.5a_hg19.gz,gnomAD_exomes_AF,gnomAD_exomes_NFE_AF,1000Gp3_AF,1000Gp3_EUR_AF,ExAC_AF,ExAC_EAS_AF,ExAC_NFE_AF,ExAC_Adj_AF,rs_dbSNP150,phyloP20way_mammalian,phyloP20way_mammalian_rankscore,phastCons20way_mammalian,phastCons20way_mammalian_rankscore,GERP++_RS,GERP++_RS_rankscore,LRT_pred,MutationTaster_pred,MutationAssessor_pred,FATHMM_pred,PROVEAN_pred,MetaLR_pred,MetaSVM_pred,M-CAP_pred,Interpro_domain,GTEx_V6p_gene,GTEx_V6p_tissue \
--plugin MaxEntScan,$PLUGIN_DBS/maxEntScan \
--custom ${PLUGIN_DBS}/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz,gnomADg,vcf,exact,0,AF_NFE,AF_Male,AF_Female,Hom,POPMAX,AF_POPMAX \
--custom ${PLUGIN_DBS}/Kaviar-160204-Public/vcfs/Kaviar-160204-Public-hg19.vcf.gz,kaviar,vcf,exact,0,AF,AC,AN \
--custom ${CCS_DB},gnomAD_exomes_CCR,bed,overlap,0,ccr_pct \
--plugin CADD,${PLUGIN_DBS}/InDels.tsv.gz,${PLUGIN_DBS}/whole_genome_SNVs.tsv.gz \
--custom ${MAF_FJD_COHORT},FJD_MAF,vcf,exact,0,AF,AN \
--custom ${DENOVO_DB},denovoVariants,vcf,exact,0,SAMPLE_CT \
-i ${VCF_FILTERED} -o ${VCF_ANNOTATED}


if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf '\nVEP ANNOTATION for '${name}' DONE\n' 

else
	printf "\nERROR: PROBLEMS WITH VEP ANNOTATION"
	exit 1
fi









printf '\n\nFiltering by population frequencies...\n'


perl $FILTER_VEP \
-i $VCF_ANNOTATED -o $VCF_ANNOTATED_POPFILT \
--filter "MAX_AF < 0.05 or not MAX_AF" \
--force_overwrite



if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf '\nVEP FREQUENCY FILTERING for '${name}' DONE\n' 
	#rm $VCF_ANNOTATED

else
	printf "\nERROR: PROBLEMS WITH VEP FREQUENCY FILTERING"
	exit 1
fi







end=`date +%s`
runtime=$((end-start))
printf '\nExecuting time: '$runtime 












