#!/bin/bash




module load python/2.7.15
module load perl
source ~/.Renviron


bed="/home/proyectos/bioinfo/lodela/SVbenchmark/CODEX2/exome_pull_down_targets/20130108.exome.targets.bed"


sbatch  --account=bioinfo_serv --partition=fastbioinfo exomeDepthWES.sh /home/proyectos/bioinfo/lodela/SVbenchmark/CODEX2/bams_ok/ /home/proyectos/bioinfo/lodela/SVbenchmark/CODEX2/individual_copynumber $bed

sbatch  --account=bioinfo_serv --partition=fastbioinfo codex2WES.sh /home/proyectos/bioinfo/lodela/SVbenchmark/CODEX2/bams_ok/ /home/proyectos/bioinfo/lodela/SVbenchmark/CODEX2/individual_copynumber $bed














