#!/bin/bash
#SBATCH --account=bioinfo_serv
#SBATCH --partition=bioinfo
#SBATCH --job-name=wesC2   #job name
#SBATCH --mail-type=END # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=ldelafuente.lorena@gmail.com # Where to send mail        
#SBATCH --mem-per-cpu=5gb # Per processor memory
#SBATCH --cpus-per-task=4
#SBATCH -t 15:00:00     # Walltime
#SBATCH -o %j_C2.out # Name output file 
#SBATCH --error=%j_C2.err
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

Rscript --vanilla /home/proyectos/bioinfo/lodela/VariantCallingFJD/utilities/CODEX2_WES_loadedInput_null.R -d $bam -o $output -n C2 -b $panel 



