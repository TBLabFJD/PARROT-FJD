#!/bin/bash
#SBATCH --account=bioinfo_serv
#SBATCH --partition=bioinfo
#SBATCH --job-name=wesED   #job name
#SBATCH --mail-type=END # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=ldelafuente.lorena@gmail.com # Where to send mail        
#SBATCH --mem-per-cpu=5gb # Per processor memory
#SBATCH --cpus-per-task=5
#SBATCH -t 15:00:00     # Walltime
#SBATCH -o %j_ED.out # Name output file 
#SBATCH --error=%j_ED.err
##SBATCH --file=
##SBATCH --initaldir=

# modules
module load R
source ~/.Renviron

# arguments

bam=$1
output=$2
panel=$3

# command line

Rscript --vanilla /home/proyectos/bioinfo/lodela/VariantCallingFJD/utilities/exomeDepth_WES.R -d $bam -o $output -b $panel -n ED -s all



