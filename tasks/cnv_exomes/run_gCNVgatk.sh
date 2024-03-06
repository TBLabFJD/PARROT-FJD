#!/bin/bash
output="/home/proyectos/bioinfo/lodela/SVbenchmark/CODEX2/individual_copynumber/copy_number_variation_data/gatk"

scatteroutput=${output}/scatter


for n in {1..89}; do 


	echo $n


	sbatch  --account=bioinfo_serv --partition=fastbioinfo /home/proyectos/bioinfo/lodela/VariantCallingFJD/utilities/gCNV_gatk.sh $n


done

