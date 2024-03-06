# NextVariantFJD
This is a germline variant calling pipeline implemented in Nextflow which performs mapping, SNV/INDEL calling and annotation, and CNV calling and annotation for targeted sequencing (gene panels and WES) and whole genome sequencing. These differen


## How to run this pipeline
The different tasks previously mention are divided into different workflows which are specified usig the `--analysis` flag followed by the corresponding letters:
 - D (Download): It downloads the FASTQ files from BaseSpace. If CNV calling or no samples are specified all samples from a project will be downloaded. 
 - M (Mapping): Specified FASTQ files from a directory (or the ones downloaded) are mapped into analysis-ready BAM files.
 - S (SNV/INDEL calling): Specified BAM files (or the ones just mapped) are used for SNV and INDEL calling using GATK, Dragen and DeepVariant workflows and merged into a single-sample VCF (single VCFs from each variant caller are also available)
 - A (Annotation of SNVs and INDELS). Specified VCF files from a directory (or the ones just generated in the SNV/INDEL calling step) are annotated and transformed into a TSV file.
 - C (CNV calling and annotation). Specified BAM files (or the ones just mapped) are used for CNV calling using Exomedepth, Convading, Panelcn.mops and GATK, and annotation using AnnotSV. In the case of analycing WGS (`--capture G`) the variant calling used is Manta.


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
