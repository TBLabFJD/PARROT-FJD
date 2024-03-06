#!/bin/sh

######################################
### Task: Mapping and unmapped BAM ###
######################################

# Mapping Fastq
# Generate uBam
# Merge bam and ubam

###############
## ARGUMENTS ##
###############

run=$1
MDAP=$2
sample=$3
forward=$4
reverse=$5
ref=$6
threads=$7
softwarePath=$8

#####################
## MODULES AND DBs ##
#####################



module load bwa/0.7.17
module load gatk/4.2.0 
module load samtools/1.9
source ${softwarePath}/pipeline.config

alias picard="java -jar ${picard_path}"
alias gatk="java -jar ${gatkPath_path}"


softwareFile="${MDAP}/software_${run}.txt"
title="MAPPING"
if [ ! -f $softwareFile ] || ! grep -q $title $softwareFile  ; then
	printf "MAPPING:\n" >> ${softwareFile}
	module list 2>> ${softwareFile}

fi






###############
## VARIABLES ##
###############



# Alignment fastq folder
AF=$MDAP/alignmentfastq 
mkdir $AF


# Temporal Folder
TMP=$MDAP/${sample}_tmp
mkdir $TMP








#############
## COMMAND ##
#############




printf "\n- BWA "

start=`date +%s`

bwa mem -v 3 -t $threads -Y \
$ref \
$forward \
$reverse |  samtools view -1 > $AF/${sample}.mapped.bam


if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf  '\nBWA MEM '${sample}'  DONE'
else
	printf "\nERROR: PROBLEMS WITH MAPPING"
	exit 1
fi

end=`date +%s`
runtime=$((end-start))
printf '\nExecuting time: '$runtime 






printf "\n- FASTQ to uBWA "


start=`date +%s`


gatk FastqToSam --TMP_DIR $TMP \
--FASTQ ${forward} \
--FASTQ2 ${reverse} \
--OUTPUT $AF/${sample}.unmapped.bam \
--READ_GROUP_NAME ${sample} \
--SAMPLE_NAME ${sample} \
--PLATFORM illumina 


if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf  '\nuBAM generation '${sample}'  DONE'
	if [ "$cat" = "True" ] || [ "$basespace" = "True" ]; then
		rm $forward
		rm $reverse
	fi
else
	printf "\nERROR: PROBLEMS WITH FASTQ TO SAM"
	exit 1
fi

end=`date +%s`
runtime=$((end-start))
printf '\nExecuting time: '$runtime 












printf "\n- Merge uBAM and generate BWA-aligned BAM "


# Merge original input uBAM file with BWA-aligned BAM file

start=`date +%s`

gatk MergeBamAlignment --TMP_DIR $TMP \
  --VALIDATION_STRINGENCY SILENT \
  --EXPECTED_ORIENTATIONS FR \
  --ATTRIBUTES_TO_RETAIN X0 \
  --ALIGNED_BAM $AF/${sample}.mapped.bam \
  --UNMAPPED_BAM $AF/${sample}.unmapped.bam  \
  --OUTPUT $AF/${sample}.mapped.merged.bam \
  --REFERENCE_SEQUENCE $ref \
  --PAIRED_RUN true \
  --SORT_ORDER "unsorted" \
  --IS_BISULFITE_SEQUENCE false \
  --ALIGNED_READS_ONLY false \
  --CLIP_ADAPTERS false \
  --ADD_MATE_CIGAR true \
  --MAX_INSERTIONS_OR_DELETIONS -1 \
  --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
  --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
  --ALIGNER_PROPER_PAIR_FLAGS true \
  --UNMAP_CONTAMINANT_READS true


if [ "$?" = "0" ]; then
	printf '\nEXIT STATUS: 0'
	printf  '\nuBAM and BWA-aligned BAM merging '${sample}'  DONE'
	rm $AF/${sample}.mapped.bam
	rm $AF/${sample}.unmapped.bam

else
	printf "\nERROR: PROBLEMS WITH MERGE BAM ALIGNMENT"
	exit 1
fi

end=`date +%s`
runtime=$((end-start))
printf '\nExecuting time: '$runtime 



# remove temporal folder

rm -r $TMP

