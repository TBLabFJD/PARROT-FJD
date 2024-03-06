#!/bin/sh

#########################
### TASK: Merge GVCFs ###
#########################

# Combines multiple variant files into a single variant file (used for combination of interval-based vcf from WGS into one vcf)


###########
## ARGUMENTS ##
###########

filestocombine=$1
output=$2
run=$3
MDAP=$4




###########
## MODULES ##
###########



module load samtools/1.9
module load picard/2.18.9
module load gatk/4.2.0
module load bedtools/2.27.0
module load R/R
alias picard='java -jar /usr/local/bioinfo/picard-tools/2.18.9/picard.jar'
alias gatk='java -jar /usr/local/bioinfo/gatk/4.2.0/gatk-package-4.2.0.0-local.jar'	


softwareFile="${MDAP}/software_${run}.txt"
title="SNV CALLING"
if [ ! -f $softwareFile ] || ! grep -q $title $softwareFile  ; then

	printf "SNV CALLING:\n" >> ${softwareFile}
	module list 2>> ${softwareFile}

fi






###########
## INPUT ##
###########


### RUN:


# Merge GVCFs generated per-interval for the same sample


echo -e "\n\n\n- merge VCFS/GVCF (GATK) "
echo -e "-----------------------------\n"

echo -e '\nUsing GATK mergeGVCF for merging VCFs/GVCFs from different intervals/chromosomes'
start=`date +%s`

gatk MergeVcfs  \
-I $filestocombine \
-O $output


if [ "$?" = "0" ]; then
	echo -e  '\nEXIT STATUS: 0'
	echo -e '\nGATK mergeGVCF DONE'

else
	echo -e  "ERROR: PROBLEMS WITH mergeGVCF"
	exit 1
fi

end=`date +%s`
runtime=$((end-start))
echo -e  '\nExecuting time: '$runtime 


