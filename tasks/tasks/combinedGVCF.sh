#!/bin/sh

############################
### TASK: Combined GVCFs ###
############################

# Combine per-sample gVCF files produced by HaplotypeCaller into a multi-sample gVCF file


###############
## ARGUMENTS ##
###############

local=$1
MDAP=$2
run=$3
fasta=$4 # Ref Fasta
softwarePath=$5



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






###############
## VARIABLES ##
###############



# genotyping folder

GEN=$MDAP/genotyping 


# Temporal Folder

TMP=$MDAP/${run}_tmp
mkdir $TMP




##############
## COMMANDS ##
##############




echo -e "\n\n\n- COMBINED GVCFS (GATK) "
echo -e "-----------------------------\n"

echo -e 'mkdir combined_gvcf_data'
echo -e '\nUsing GATK COMBINEGVCFs for combining GVCFs'
start=`date +%s`

gatk CombineGVCFs --tmp-dir $TMP \
-R $fasta \
--variant $GEN/my_list_of_gvcfs_files_to_combine_$run.list \
-O $GEN/${run}.g.vcf


if [ "$?" = "0" ]; then
	echo -e  '\nEXIT STATUS: 0'
	echo -e '\nGATK COMBINEGVCFs DONE'

else
	echo -e  "ERROR: PROBLEMS WITH COMBINEGVCFs"
	exit 1
fi

end=`date +%s`
runtime=$((end-start))
echo -e  '\nExecuting time: '$runtime 


