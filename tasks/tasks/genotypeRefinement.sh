#!/bin/sh

############################
### TASK: SNV filtering ###
###########################

# Two options: CNN or Hard Filtering


###############
## ARGUMENTS ##
###############

local=$1
run=$2
MDAP=$3
ped=$4
fasta=$5
softwarePath=$6




#####################
## MODULES AND DBs ##
#####################



module load gatk/4.2.0
source ${softwarePath}/pipeline.config

eval "$(${conda_bin} shell.bash hook)"

source activate gatk

alias gatk="java -jar ${gatkPath_path}"
#alias gatk='java -jar /usr/local/bioinfo/gatk/4.2.0/gatk-package-4.2.0.0-local.jar'	


softwareFile="${MDAP}/software_${run}.txt"
title="GENOTYPE REFINEMENT:"
if [ ! -f $softwareFile ] || ! grep -q $title $softwareFile  ; then
	printf "GENOTYPE REFINEMENT:\n" >> ${softwareFile}
	module list 2>> ${softwareFile}

fi







###############
## VARIABLES ##
###############

# snv results folder
SNV="${MDAP}/snvs"
mkdir $SNV

# Temporal Folder
TMP=$MDAP/${run}_tmp
mkdir $TMP

# Input Vcf
VCF=$SNV/${run}.gatkLabeled.vcf


##############
## COMMANDS ##
##############



echo -e "\n\n\n- GENOTYPE POSTERIOR RECALCULATION BASED ON FAMILY PRIOR INFORMATION (GATK) "
echo -e "-------------------------------------------------------------------------------\n"


if [ "$ped" != "null" ]; then

	start=`date +%s`

	# Recalculate genotype posterior probabilities based on family prior information.
	# Compute the most likely genotype combination of trios and parent/child pairs given their genotype likelihoods and a mutation prior;

	gatk CalculateGenotypePosteriors \
	   -R $fasta \
	   -V $VCF \
	   -ped $ped \
	   -O $SNV/filtered_INDEL_SNP_data_GP_$run.vcf \
	   --skip-population-priors  
	s1="$?"


	# Filter variants based on GQ: genotypes with GQ < 20 based on the posteriors are flagged for posterior filtering. 
	gatk VariantFiltration \
	-R $fasta \
	-V $SNV/filtered_INDEL_SNP_data_GP_$run.vcf \
	--filter-expression "GQ < 20" \
	--filter-name "lowGQ" \
	-O $SNV/filtered_INDEL_SNP_data_GP_GQfiltered_$run.vcf 
	s2="$?"


	# Annotation of High and Low Confidence De Novo mutations: 

	gatk VariantAnnotator \
		-R $fasta \
		-V $SNV/filtered_INDEL_SNP_data_GP_GQfiltered_$run.vcf  \
		-O ${VCF}_tmp \
		-A PossibleDeNovo  \
		-A StrandBiasBySample  \
		-A AS_FisherStrand \
		-ped $ped 
	s3="$?"



	if [ "$s1" = "0"  ] &&  [ "$s2" = "0" ] &&  [ "$s3" = "0" ] ; then
		echo -e  '\nEXIT STATUS: 0'
		echo -e '\nPosterior genotyping and de novo mutations annotation DONE\n'
		rm  $SNV/filtered_INDEL_SNP_data_GP_$run.vcf*
		rm  $SNV/filtered_INDEL_SNP_data_GP_GQfiltered_$run.vcf* 
		mv ${VCF}_tmp ${VCF}
		rm  ${VCF}_tmp*idx



	else
		echo -e  "ERROR: PROBLEMS WITH posterior genotyping and de novo mutations annotation"
		exit 1
	fi


	end=`date +%s`
	runtime=$((end-start))
	echo -e  '\nExecuting time: '$runtime 



else

	echo -e "\n NOT INPUT PED FILE: genotype posterior calculation skipped"

fi





# Removing temporal forder

rm -r $TMP