# PARROT-FJD
Pipeline of Analysis and Research of Rare diseases Optimized in Tblab - Fundación Jiménez Díaz. This is a germline variant calling pipeline implemented in Nextflow which performs mapping, SNV/INDEL calling and annotation, and CNV calling and annotation for targeted sequencing (gene panels and WES) and whole genome sequencing. 

<p align="center">
  <img width="300" height="800" src="https://github.com/TBLabFJD/PARROT-FJD/assets/48798983/309dfb16-ef3f-4885-bdb1-02a77b89f414")
">
</p>


## How to run this pipeline
The different tasks previously mention are divided into different workflows which are specified usig the `--analysis` flag followed by the corresponding letters:
 - D (Download): It downloads the FASTQ files from BaseSpace. If CNV calling or no samples are specified all samples from a project will be downloaded. 
 - M (Mapping): Specified FASTQ files from a directory (or the ones downloaded) are mapped into analysis-ready BAM files.
 - S (SNV/INDEL calling): Specified BAM files (or the ones just mapped) are used for SNV and INDEL calling using GATK by default. Dragen and DeepVariant are also available. The variant caller can be selected with the parameter `--vc_tools`. The options are: 
     - `gatk`: Haplotypecaller
     - `dragen`
     - `deepvariant`
     - `all` (equivalent to: gatk,dragen,deepvariant). Using this option, resulting vcfs will be merged into a single-sample VCF (single VCFs from each variant caller are also available)
   More than one tool can be chosen using "," (`--vc_tools gatk,dragen`)
 - A (Annotation of SNVs and INDELS). Specified VCF files from a directory (or the ones just generated in the SNV/INDEL calling step) are annotated and transformed into a TSV file.
 - C (CNV calling and annotation). Specified BAM files (or the ones just mapped) are used for CNV calling using Exomedepth, Convading, Panelcn.mops and GATK, and annotation using AnnotSV. In the case of analycing WGS (`--capture G`) the variant calling used is Manta.

Mapping and variant calling processes can be parallelized to speed up the analysis. These option can be activated using the parameters `--parallel_mapping true` and `--parallel_calling true`.

`--parallel_mapping true`: FASTP will be executed to split FASTQ files in three chunks that will be mapped in parallel.   
`--parallel_calling true`: BAM file will be split by chromosomes in smaller BAM files that will be processed in parallel. 

For using this options, using `--cpus-per-task=44` is recommended. 

You can generate and keep a cram file out of your bam, when running either MS or just S. The cram is generated inside the out folder: /out/cram/. By default: --keep_cram false
`--keep_cram true`: generate and keep cram file

You can generate the mosdepth bed file from your bam, when running either MS or just S. The bed is generated inside the out folder: /out/qc/mosdepth_cov/. By default: --mosdepth_bed false
`--mosdepth_bed true`: generate the mosdepth.bed file needed to update the db of allele frequencies. 


There are different profiles available depending on the reference release to use, where to run it, and type of contenerization:

Mandatory to choose one:
 - hg19: use the reference genome hg19
 - hg38: use the reference genome hg38

Mandatory to choose one:
 - tblabserver: run pipeline in the server just for the canonical chromosomes
 - tblabserver_allcontigs: run pipeline in the server just for all contigs
 - uam: run pipeline in the CCC (UAM) just for the canonical chromosomes
 - uam_allcontigs: run pipeline in the CCC (UAM) just for all contigs

Mandatory to choose one:
 - docker: run pipeline using docker containers
 - singularity: run pipeline using singularity containers

Mandatory when running in the CCC (UAM):
 - batch: run pipeline in for slurm executor in the CCC (UAM)


Profiles to use in the CCC (UAM): 
 - `-profile hg19,singularity,uam,batch` 
 - `-profile hg38,singularity,uam,batch`
 - `-profile hg19,singularity,uam_allcontigs,batch` 
 - `-profile hg38,singularity,uam_allcontigs,batch`

Profiles to use in the server:
 - `-profile hg19,singularity,tblabserver` 
 - `-profile hg38,singularity,tblabserver`
 - `-profile hg19,docker,tblabserver` 
 - `-profile hg38,docker,tblabserver`
 - `-profile hg19,singularity,tblabserver_allcontigs` 
 - `-profile hg38,singularity,tblabserver_allcontigs`
 - `-profile hg19,docker,tblabserver_allcontigs` 
 - `-profile hg38,docker,tblabserver_allcontigs`

Check the file nexflow.config to see the description of all arguments.

## Examples

Perform a complete analysis from downloading samples from Basespace to SNV/INDEL and CNV calling and annotation.
```
nextflow run /home/gonzalo/nextflowtest/NextVariantFJD/main.nf \
-profile hg38,singularity,tblabserver --analysis DMSAC \
--input project_name \
--output /output/path/ \
--bed /path/to/captured/regions.bed \
-with-report report.html
```


Perform an analysis from mapping to SNV/INDEL calling and annotation.
```
nextflow run /home/gonzalo/nextflowtest/NextVariantFJD/main.nf \
-profile hg38,singularity,tblabserver --analysis MSA \
--input /input/path/to/fastq/ \
--output /output/path/ \
-with-report report.html
```


Perform a complete analysis from downloading samples from Basespace to SNV/INDEL and CNV calling and annotation. Variant calling analysis is performed to a subset of samples and genes.
```
nextflow run /home/gonzalo/nextflowtest/NextVariantFJD/main.nf \
-profile hg38,singularity,tblabserver --analysis DMSAC \
--input project_name \
--output /output/path/ \
--bed /path/to/captured/regions.bed \
--samples /path/to/samplefile.txt \
--genelist /path/to/genelist.txt \
-with-report report.html
```


Perform a SNV/INDEL calling and annotation analysis. Variant calling analysis is performed to a subset of samples. Genes are prioritized using genelista and glowgenes.
```
nextflow run /home/gonzalo/nextflowtest/NextVariantFJD/main.nf \
-profile hg38,singularity,tblabserver --analysis SA \
--input /input/path/to/vcfs/ \
--output /output/path/ \
--samples /path/to/samplefile.txt \
--genelist /path/to/genelist.txt \
--glowgenes /path/to/glowgenes_result.txt \
-with-report report.html
```

## Versions

GATK 4.4.0.0
VEP release 105
deepvariant v1.4.0
annotsv 3.1.1
manta 1.6.0
BWA 0.7.17-r1198-dirty
bcftools 1.15
bedtools 2.27.1

R version 4.2.3
python 3.6


