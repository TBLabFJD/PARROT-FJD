#!/bin/bash
#SBATCH --account=bioinfo_serv
#SBATCH --partition=bioinfo
#SBATCH --job-name=wesGATK   #job name
#SBATCH --mail-type=FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ldelafuente.lorena@gmail.com # Where to send mail        
#SBATCH --mem-per-cpu=3gb # Per processor memory
#SBATCH --cpus-per-task=1
#SBATCH -t 15:00:00     # Walltime
#SBATCH -o %j_GATK.out # Name output file 
#SBATCH --error=%j_GATK.err
##SBATCH --file=
##SBATCH --initaldir=


#https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-common-and-rare-germline-copy-number-variants


# modules
module load gatk/4.1.5.0
eval "$(/usr/local/miniconda/python-3.6/bin/conda shell.bash hook)"
source activate gatk
#alias gatk='java -jar /usr/local/bioinfo/gatk/4.1.5.0/gatk-package-4.1.5.0-local.jar'   
#alias gatk='java -jar /usr/local/bioinfo/gatk/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar'   



# arguments

BED='/home/proyectos/bioinfo/lodela/SVbenchmark/CODEX2/exome_pull_down_targets/20130108.exome.targets_woXY.bed'
REF="/home/proyectos/bioinfo/references/hg19/ucsc.hg19.fasta"
REFdict="/home/proyectos/bioinfo/references/hg19/ucsc.hg19.dict"
input="/home/proyectos/bioinfo/lodela/SVbenchmark/CODEX2/bams_ok/"
output="/home/proyectos/bioinfo/lodela/SVbenchmark/CODEX2/individual_copynumber/copy_number_variation_data/gatk"
priors="/home/proyectos/bioinfo/references/gCNVgatk/contig_ploidy_priors.tsv"



 # 1. Collect raw counts data with PreprocessIntervals and CollectReadCounts


# For exome data, pad target regions, e.g. with 250 bases.

bedsuffix=$(basename $BED .bed)

pBED=${output}/${bedsuffix}.preprocessed.interval_list
GCannot=${output}/${bedsuffix}.annotated.mytsv
fBED=${output}/cohort.gc.filtered.interval_list

# echo $pBED
# echo $GCannot

# gatk PreprocessIntervals \
#         -R $REF \
#         -L $BED \
#         --bin-length 0 \
#         -imr OVERLAPPING_ONLY \
#         -O ${pBED}


#Count reads per bin using CollectReadCounts

tagI=""
for bam in ${input}/*bam; do

    countfile=$(basename $bam _alignment.bam)
    countoutput=${output}/${countfile}.tsv
    # gatk CollectReadCounts \
    #     -L ${pBED} \
    #     -R $REF \
    #     -imr OVERLAPPING_ONLY \
    #     -I $bam \
    #     --format TSV \
    #     -O $countoutput

    tagI="${tagI} -I ${countoutput}"

done


## 2. (Optional) Annotate intervals with features and subset regions of interest with FilterIntervals



## AnnotateIntervals with GC content

# gatk AnnotateIntervals \
#     -L ${pBED} \
#     -R $REF \
#     -imr OVERLAPPING_ONLY \
#     -O $GCannot
#This produces a four-column table where the fourth column gives the fraction of GC content.



## FilterIntervals based on GC-content and cohort extreme counts

# gatk FilterIntervals \
#     -L ${pBED} \
#     --annotated-intervals $GCannot \
#     $tagI  \
#     -imr OVERLAPPING_ONLY \
#     -O ${fBED}




# 3. Call autosomal and allosomal contig ploidy with DetermineGermlineContigPloidy

scatteroutput=${output}/scatter
mkdir $scatteroutput


# gatk IntervalListTools \
#         --INPUT ${fBED} \
#         --SUBDIVISION_MODE INTERVAL_COUNT \
#         --SCATTER_CONTENT 5000 \
#         --OUTPUT $scatteroutput 



# gatk DetermineGermlineContigPloidy \
#         -L ${fBED} \
#         --interval-merging-rule OVERLAPPING_ONLY \
#         $tagI \
#         --contig-ploidy-priors ${priors} \
#         --output $output \
#         --output-prefix ploidy \
#         --verbosity DEBUG



# 4. Call copy number variants with GermlineCNVCaller


gCNVoutput=${output}/cohort
mkdir $gCNVoutput
tagM=""
tagC=""

#for dir in ${scatteroutput}/*; do

#dir=$1
#prefix=$(basename $dir)
#echo $prefix

# gatk GermlineCNVCaller \
#         --run-mode COHORT \
#         -L $dir/scattered.interval_list \
#         $tagI \
#         --contig-ploidy-calls ploidy-calls \
#         --annotated-intervals $GCannot  \
#         --interval-merging-rule OVERLAPPING_ONLY \
#         --output $gCNVoutput \
#         --output-prefix $prefix \
#         --verbosity DEBUG

#tagM="${tagM} --model-shard-path ${gCNVoutput}/${prefix}-model"
#tagC="${tagC} --calls-shard-path ${gCNVoutput}/${prefix}-model"

#done



# # 5. Call copy number segments and consolidate sample results with PostprocessGermlineCNVCalls

n=$1
ls -ltrh

tagC=""
for file in cohort/temp_00*-calls; do tagC="${tagC} --calls-shard-path ${file}"; done
tagM=""
for file in cohort/temp_00*-model; do tagM="${tagM} --model-shard-path ${file}"; done


gatk PostprocessGermlineCNVCalls \
    $tagM \
    $tagC \
    --allosomal-contig chrX --allosomal-contig chrY \
    --contig-ploidy-calls ploidy-calls \
    --sample-index $n \
    --output-genotyped-intervals genotyped-intervals-cohort90-twelve-${n}.vcf.gz \
    --output-genotyped-segments genotyped-segments-cohort90-${n}.vcf.gz \
    --sequence-dictionary ${REFdict} \
    --output-denoised-copy-ratios sample_${n}_denoised_copy_ratios.tsv \
    --autosomal-ref-copy-number 2











