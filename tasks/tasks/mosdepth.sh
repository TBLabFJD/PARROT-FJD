#!/bin/sh

########################################################################################
### TASK: Bam coverage calculation for QC (in SNP and CNV) and internal MAF database ###
########################################################################################


###############
## ARGUMENTS ##
###############

local=$1
run=$2
MDAP=$3
sample=$4
bamfile=$5
panel=$6
padding=$7
analysis=$8
mycountthreshold=$9


#####################
## MODULES AND DBs ##
#####################


		
module load miniconda/3.6
module load bedtools/2.27.0

softwareFile="${MDAP}/software_${run}.txt"
title="COVERAGE QC"

if [ ! -f $softwareFile ] || ! grep -q $title $softwareFile  ; then

	printf "COVERAGE QC:\n" >> ${softwareFile}
	module list 2>> ${softwareFile}	
	mosdepth --version >> ${softwareFile}

fi










###############
## VARIABLES ##
###############


# snv results folder
QC="${MDAP}/qc"
mkdir $QC


# Percentage of covered regions 
covPerc=0.9




##############
## COMMANDS ##
##############


# 0. Definition of read threshold to check quality. QC parameters will be different dependind on SNV or CNVs

if [ "$analysis" != "cnv" ]; then readthreshold=10; else readthreshold=$mycountthreshold; fi





# 1. QC to check coverage for input bams. Writing results to File

cov=0

if [ "$panel" != "genome" ]; then


	# run mosdepth to check QC at panel regions
	mosdepth -x --no-per-base --by $panel ${QC}/${sample}_qc_$analysis ${bamfile}
	s1="$?"

	cov=$(awk -v reads="$readthreshold" '{if($1=="total" && $2==reads){print $3}}' ${QC}/${sample}_qc_${analysis}.mosdepth.region.dist.txt)
	s2="$?"



else

	# run mosdepth to check QC (whole genome)
	mosdepth -x --no-per-base ${QC}/${sample} ${bamfile}
	s1="$?"

	cov=$(awk -v reads="$readthreshold" '{if($1=="total" && $2==reads){print $3}}' ${QC}/${sample}_qc_${analysis}.mosdepth.global.dist.txt)
	s2="$?"



fi




if [ "$s1" = "0" ]  &&  [ "$s2" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf '\nCOVERAGE QC for '${sample}' DONE\n' 
	
	#rm ${QC}/${sample}_qc*

	pass=$(awk 'BEGIN{if ('$cov'>'$covPerc') print 0}')

	if [ "$pass" = "0" ]; then
		echo -e $sample"\t"$readthreshold"\t"$cov"\tPASSED\t"$run"\t"$bamfile >> ${QC}/minCovFilterResults_${analysis}_${run}.txt
		printf 'COVERAGE LEVELS for '${sample}' ARE OK.\n'
	else
		printf 'COVERAGE LEVELS for '${sample}' ARE LOW.\n'
		echo -e $sample"\t"$readthreshold"\t"$cov"\tFAILED\t"$run"\t"$bamfile >> ${QC}/minCovFilterResults_${analysis}_${run}.txt
		
		printf '\nERROR ***************SKIPPING '${analysis}' ANALYSIS FOR SAMPLE '${sample}'.\n'
		exit 2
	fi
else
	printf "\nERROR: PROBLEMS WITH COVERAGE QC"
	exit 1
fi









# 2. If coverage is OK for sample and in case of SNVs analysis we need to compute coverage for MAF ddbb


if [ "$analysis" != "cnv" ]; then


	mosdepth --quantize 10: -n -x  ${QC}/${sample}_cov ${bamfile}
	s1="$?"

	# in case of WES or panels we generate bed file with padding to consider only those regions

	if [ "$panel" != "genome" ]; then
		
		# create panel with padding regions: regions where SNVs are called
		bedpadding="$QC/$(basename -- $panel)_padding_${sample}.bed"
		awk  -v pad="$padding" '{print $1"\t"$2-pad"\t"$3+pad"\t"$4}' $panel | bedtools merge -i - > $bedpadding
		s2="$?"

		# we just keep coverage at padding regions
		bedtools intersect -a ${bedpadding} -b ${QC}/${sample}_cov.quantized.bed.gz -wb | awk '{print $1"\t"$2"\t"$3"\t"$NF}'  > ${QC}/${sample}_padding.quantized.bed
		s3="$?"

	fi
	




	if [ "$s1" = "0" ]  &&  [ "$s2" = "0" ]  &&  [ "$s3" = "0" ]; then
		printf '\nEXIT STATUS: 0'
		printf '\nTARGET COVERAGE FOR '${sample}' DONE\n' 
		rm ${QC}/${sample}_cov*
		rm $bedpadding



	else
		printf "\nERROR: PROBLEMS WITH TARGET COVERAGE CALCULATION"
		exit 1
	fi





fi








