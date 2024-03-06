#!/bin/sh

#########################
### TASK: SNV calling ###
#########################

# Haplotipe Caller
# Genotype Caller


###############
## ARGUMENTS ##
###############

local=$1
run=$2
MDAP=$3
sample=$4
fasta=$5 
bamfile=$6
intervals=$7
panel=$8
padding=$9
cvcf=${10}
removebam=${11}
mem=${12}
softwarePath=${13}

#####################
## MODULES AND DBs ##
#####################


module load samtools/1.9
module load picard/2.18.9
module load gatk/4.2.0
module load bedtools/2.27.0
module load R/R
source ${softwarePath}/pipeline.config

alias picard="java -jar -Xmx${mem}g ${picard_path}"
alias gatk="java -jar -Xmx${mem}g ${gatkPath_path}"
#alias picard="java -jar  -Xmx${mem}g /usr/local/bioinfo/picard-tools/2.18.9/picard.jar"
#alias gatk="java -jar -Xmx${mem}g /usr/local/bioinfo/gatk/4.2.0/gatk-package-4.2.0.0-local.jar"

softwareFile="${MDAP}/software_${run}.txt"
title="SNV CALLING"
if [ ! -f $softwareFile ] || ! grep -q $title $softwareFile  ; then
	printf "SNV CALLING:\n" >> ${softwareFile}
	module list 2>> ${softwareFile}

fi






###########
## INPUT ##
###########



# genotyping folder

GEN=$MDAP/genotyping 
mkdir $GEN



# Temporal Folder

TMP=$MDAP/${sample}_tmp
mkdir $TMP





#############
## COMMAND ##
#############




printf "\n\n\n- HAPLOTYPECALLER (GATK)"





printf '\nGATK HaplotypeCallerGVCF for '${sample}' STARTS'
start=`date +%s`

if [ "$cvcf" != "True" ]; then 

	if [ "$intervals" != "True" ]; then
		gatk HaplotypeCaller --tmp-dir $TMP \
		-R $fasta \
		-I $bamfile \
		-bamout $GEN/${sample}.bamout.bam \
		-O $GEN/${sample}.vcf \
		--annotate-with-num-discovered-alleles true 
	else
		gatk HaplotypeCaller --tmp-dir $TMP \
		-R $fasta \
		-I $bamfile \
		-bamout $GEN/${sample}.bamout.bam \
		-O $GEN/${sample}.vcf \
		--annotate-with-num-discovered-alleles true \
		-L $panel -ip $padding
	fi

	if [ "$?" = "0" ]; then
		printf '\nEXIT STATUS: 0'
		printf '\nGATK HaplotypeCaller in VCF mode for '${sample}' DONE\n' 
		print $removebam
		if [ "$removebam" = "single" ]; then
		 	basename=${bamfile%.*}
		 	rm ${basename}.bam
		 	rm ${basename}.bai
		fi

	else
		printf "\nERROR: PROBLEMS WITH HAPLOTYPECALLER"
		exit 1
	fi


else

	if [ "$intervals" != "True" ]; then
		gatk HaplotypeCaller --tmp-dir $TMP \
		-R $fasta \
		-I $bamfile \
		-ERC GVCF \
		-bamout $GEN/${sample}_bamout.bam \
		-O $GEN/${sample}.g.vcf \
		-G StandardAnnotation \
		-G AS_StandardAnnotation \
		-G StandardHCAnnotation \
		-A FisherStrand -A StrandOddsRatio -A RMSMappingQuality -A MappingQualityRankSumTest -A ReadPosRankSumTest -A DepthPerSampleHC -A BaseQualityRankSumTest -A ExcessHet \
		--annotate-with-num-discovered-alleles true 
	else
		gatk HaplotypeCaller --tmp-dir $TMP \
		-R $fasta \
		-I $bamfile \
		-ERC GVCF \
		-bamout $GEN/${sample}_bamout.bam \
		-O $GEN/${sample}.g.vcf \
		-G StandardAnnotation \
		-G AS_StandardAnnotation \
		-G StandardHCAnnotation \
		-A FisherStrand -A StrandOddsRatio -A RMSMappingQuality -A MappingQualityRankSumTest -A ReadPosRankSumTest -A DepthPerSampleHC -A BaseQualityRankSumTest -A ExcessHet \
		--annotate-with-num-discovered-alleles true \
		-L $panel -ip $padding
	fi

	if [ "$?" = "0" ]; then
		printf '\nEXIT STATUS: 0'
		printf '\nGATK HaplotypeCaller in GVCF mode for '${sample}' DONE\n' 

		if [ "$removebam" = "single" ]; then
		 	basename=${bamfile%.*}
		 	rm ${basename}.bam
		 	rm ${basename}.bai
		fi

	else
		printf "\nERROR: PROBLEMS WITH HAPLOTYPECALLER"
		exit 1
	fi


	printf '\n'$GEN/${sample}.g.vcf >> $GEN/my_list_of_gvcfs_files_to_combine_$run.list

fi





end=`date +%s`
runtime=$((end-start))
printf '\nExecuting time: '$runtime 



# remove temporal folder

rm -r $TMP















# printf "\n\n\n- GENOTYPECALLER (GATK) "
# printf "\n--------------------------\n"


# #GenotypeGVCFs into final VCF

# printf  '\nUsing GATK GenotypeGVCFs for final VCF'
# start=`date +%s`

# gatk GenotypeGVCFs --tmp-dir=$TMP \
# 	-R $fasta \
# 	-V $GEN/${sample}.g.vcf \
# 	-G StandardAnnotation \
# 	-O $GEN/${sample}.vcf 

# if [ "$?" = "0" ]; then
# 	printf '\nEXIT STATUS: 0'
# 	printf  '\nGATK GenotypeGVCFs for '${sample}' DONE'

# else
# 	printf "\nERROR: PROBLEMS WITH GENOTYPECALLER"
# 	exit 1
# fi


# end=`date +%s`
# runtime=$((end-start))
# printf '\nExecuting time: '$runtime 




