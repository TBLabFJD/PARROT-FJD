#!/bin/bash
#SBATCH --account=bioinfo_serv
#SBATCH --partition=bioinfo
#SBATCH --job-name=wesXHMM   #job name
#SBATCH --mail-type=END # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=ldelafuente.lorena@gmail.com # Where to send mail        
#SBATCH --mem-per-cpu=5gb # Per processor memory
#SBATCH --cpus-per-task=1
#SBATCH -t 15:00:00     # Walltime
#SBATCH -o %j_XHMM.out # Name output file 
#SBATCH --error=%j_XHMM.err
##SBATCH --file=
##SBATCH --initaldir=


# modules

module load R/3.5.0
module load python/3.6.5
module load gatk/3.8
module load xhmm



bamlist=$1 
output=$2
tag=$3
REF="/home/proyectos/bioinfo/references/hg19/ucsc.hg19.fasta"
BED='/home/proyectos/bioinfo/lodela/SVbenchmark/CODEX2/exome_pull_down_targets/20130108.exome.targets_woXY.bed'


#java -jar /usr/local/bioinfo/gatk/3.8/GenomeAnalysisTK.jar \
#-T DepthOfCoverage \
#-I ${bamlist} -L $BED \
#-R $REF \
#-dt BY_SAMPLE -dcov 5000 -l INFO --omitDepthOutputAtEachBase --omitLocusTable \
#--minBaseQuality 0 --minMappingQuality 20 --start 1 --stop 5000 --nBins 200 \
#--includeRefNSites \
#--countType COUNT_FRAGMENTS \
#-o ${output}

# tagI=""

# for file in /home/proyectos/bioinfo/lodela/SVbenchmark/CODEX2/individual_copynumber/copy_number_variation_data/XHMM/*interval_summary;do

# 	tagI="${tagI} --GATKdepths $file"

# done

# echo $tagI




# # Combines GATK Depth-of-Coverage outputs for multiple samples (at same loci):

# xhmm --mergeGATKdepths -o ${output}/DATA.1000G.txt \
# $tagI


#Optionally, run GATK to calculate the per-target GC content and create a list of the targets with extreme GC content:

java -jar /usr/local/bioinfo/gatk/3.8/GenomeAnalysisTK.jar \
-T GCContentByInterval -L $BED \
-R $REF \
-o $output/DATA.locus_GC.txt

cat $output/DATA.locus_GC.txt | awk '{if ($2 < 0.1 || $2 > 0.9) print $1}' \
>  ${output}/extreme_gc_targets.txt




#Filters samples and targets and then mean-centers the targets: 
xhmm --matrix -r  ${output}/DATA.1000G.txt --centerData --centerType target \
-o ${output}/DATA.filtered_centered.${tag}.txt \
--outputExcludedTargets ${output}/DATA.filtered_centered.${tag}.txt.filtered_targets.txt \
--outputExcludedSamples ${output}/DATA.filtered_centered.${tag}.txt.filtered_samples.txt \
--minTargetSize 10 --maxTargetSize 10000 \
--minMeanTargetRD 10 --maxMeanTargetRD 500 \
--minMeanSampleRD 25 --maxMeanSampleRD 200 \
--maxSdSampleRD 150 \
--excludeTargets ${output}/extreme_gc_targets.txt 


#Runs PCA on mean-centered data:
xhmm --PCA -r ${output}/DATA.filtered_centered.${tag}.txt --PCAfiles ${output}/DATA.${tag}_PCA
echo "PHASE: 1"
echo $0

#Normalizes mean-centered data using PCA information:
xhmm --normalize -r  ${output}/DATA.filtered_centered.${tag}.txt --PCAfiles ${output}/DATA.${tag}_PCA \
--normalizeOutput ${output}/DATA.PCA_normalized.txt \
--PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7
echo "PHASE: 2"


#Filters and z-score centers (by sample) the PCA-normalized data:
xhmm --matrix -r  ${output}/DATA.PCA_normalized.txt --centerData --centerType sample --zScoreData \
-o  ${output}/DATA.PCA_normalized.filtered.sample_zscores.${tag}.txt \
--outputExcludedTargets  ${output}/DATA.PCA_normalized.filtered.sample_zscores.${tag}.txt.filtered_targets.txt \
--outputExcludedSamples  ${output}/DATA.PCA_normalized.filtered.sample_zscores.${tag}.txt.filtered_samples.txt \
--maxSdTargetRD 30
echo "PHASE: 3"
echo $0


#Filters original read-depth data to be the same as filtered, normalized data:
xhmm --matrix -r ${output}/DATA.1000G.txt \
--excludeTargets ${output}/DATA.filtered_centered.${tag}.txt.filtered_targets.txt \
--excludeTargets ${output}/DATA.PCA_normalized.filtered.sample_zscores.${tag}.txt.filtered_targets.txt \
--excludeSamples ${output}/DATA.filtered_centered.${tag}.txt.filtered_samples.txt \
--excludeSamples ${output}/DATA.PCA_normalized.filtered.sample_zscores.${tag}.txt.filtered_samples.txt \
-o ${output}/DATA.same_filtered.${tag}.txt

echo "PHASE: 4"
echo $0


#Discovers CNVs in normalized data:
xhmm --discover -p ${output}/params.txt \
-r ${output}/DATA.PCA_normalized.filtered.sample_zscores.${tag}.txt -R ${output}/DATA.same_filtered.${tag}.txt \
-c ${output}/DATA.xcnv -a ${output}/DATA.aux_xcnv -s ${output}/DATA

echo $0
echo "PHASE: 5"

#NOTE: A description of the .xcnv file format can be found below.


#Genotypes discovered CNVs in all samples:
xhmm --genotype -p ${output}/params.txt \
-r ${output}/DATA.PCA_normalized.filtered.sample_zscores.${tag}.txt -R ${output}/DATA.same_filtered.${tag}.txt \
-g ${output}/DATA.xcnv -F $REF \
-v ${output}/DATA.vcf


echo "PHASE: 6"
echo $0
