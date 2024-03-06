#!/bin/sh

#########################
### TASK: Genotyping ###
#########################

# Perform joint genotyping on:

# 1. a singular sample by providing a single-sample GVCF
# 2. a combined multi-sample GVCF
# 3. on GenomicsDB workspace created with GenomicsDBImport


###############
## ARGUMENTS ##
###############

local=$1
run=$2
MDAP=$3
sample=$4  # sample or run name if multi-sample vcf
input=$5
output=$6
fasta=$7  # ref Fasta
ped=$8
softwarePath=$9

#####################
## MODULES AND DBs ##
#####################



module load samtools/1.9
module load picard/2.18.9
module load gatk/4.2.0
module load bedtools/2.27.0
module load R/R
source ${softwarePath}/pipeline.config

alias picard="java -jar ${picard_path}"
alias gatk="java -jar ${gatkPath_path}"

#alias picard='java -jar /usr/local/bioinfo/picard-tools/2.18.9/picard.jar'
#alias gatk='java -jar /usr/local/bioinfo/gatk/4.2.0/gatk-package-4.2.0.0-local.jar'	


softwareFile="${MDAP}/software_${run}.txt"
title="SNV CALLING"
if [ ! -f $softwareFile ] || ! grep -q $title $softwareFile  ; then
	printf "SNV CALLING:\n" >> ${softwareFile}
	module list 2>> ${softwareFile}

fi







###########
## INPUT ##
###########


# Temporal Folder

TMP=$MDAP/${sample}_tmp
mkdir $TMP





#############
## COMMAND ##
#############



printf "\n\n\n- GENOTYPECALLER (GATK) "
printf "\n--------------------------\n"


#GenotypeGVCFs into final VCF

printf  '\nUsing GATK GenotypeGVCFs'
start=`date +%s`
echo $ped

if [ "$ped" != "null" ]; then

	gatk GenotypeGVCFs --tmp-dir $TMP \
		-R $fasta \
		-V $input \
		-O $output \
		-ped $ped

else

	gatk GenotypeGVCFs --tmp-dir $TMP \
		-R $fasta \
		-V $input \
		-O $output

fi






if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf  '\nGATK GenotypeGVCFs for '${sample}' DONE'

else
	printf "\nERROR: PROBLEMS WITH GENOTYPECALLER"
	exit 1
fi


end=`date +%s`
runtime=$((end-start))
printf '\nExecuting time: '$runtime 



# Removing temporal forder

rm -r $TMP
