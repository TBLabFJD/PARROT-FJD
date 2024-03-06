#!/bin/sh

############################
### TASK: SNV filtering ###
###########################

# Two options: CNN or Hard Filtering


###############
## ARGUMENTS ##
###############

local=$1
run=$2
MDAP=$3
name=$4 # sample name or run name
REF=$5
ftype=$6 # when multi VCF is used only "HF" can be used
cvcf=$7
softwarePath=$8
#skipmapping=$9    
#INPUT=${10}




#####################
## MODULES AND DBs ##
#####################



module load gatk/4.2.0
source $softwarePath/pipeline.config

eval "$(${conda_bin} shell.bash hook)"

source activate gatk

alias gatk="java -jar ${gatkPath_path}"


#alias gatk='java -jar /usr/local/bioinfo/gatk/4.2.0/gatk-package-4.2.0.0-local.jar'	
#CNN_REF=/home/proyectos/bioinfo/references/CNNgatk


softwareFile="${MDAP}/software_${run}.txt"
title="SNV FILTERING"
if [ ! -f $softwareFile ] || ! grep -q $title $softwareFile  ; then

	printf "SNV FILTERING:\n" >> ${softwareFile}
	module list 2>> ${softwareFile}

fi









###############
## VARIABLES ##
###############

# Ref Fasta
fasta=$REF


# genotyping folder
GEN=$MDAP/genotyping 

# snv results folder
SNV="${MDAP}/snvs"
mkdir $SNV

# Temporal Folder
TMP=$MDAP/${name}_tmp
mkdir $TMP







##############
## COMMANDS ##
##############


# IMPORTANT: CNN method can be only applied to single samples since the default models were trained on single-sample VCFs. 


if [ "$ftype" != "HF" ]; then



	printf "\n\n\n- Convolutional Neural Network (CNN)"
	printf "\n---------------------------------------\n"


	# 0. Set bam folder 

	if [ "$skipmapping" != "True" ]; then bamfolder=$MDAP/bams; else bamfolder=${INPUT}; fi


	# 1. Annotate a VCF with scores

	gatk CNNScoreVariants \
	--disable-avx-check \
	-I ${bamfolder}/${name}.bam \
	-V $GEN/${name}.vcf  \
	-O $SNV/${name}.cnnScored.vcf  \
	-R $fasta \
	--tensor-type read_tensor

	#echo $SNV/${name}.cnnScored.vcf
	#echo ${CNN_REF}/b37_hapmap_3.3.b37.chr.vcf.gz


	# 2. Filtering variants according to CNN score.

	gatk FilterVariantTranches \
	-V  $SNV/${name}_cnnScored.vcf \
	--resource ${cnn_hapmap_path} \
	--resource ${cnn_1000G_path} \
	--info-key CNN_2D \
	--snp-tranche 99.95 \
	--indel-tranche 99.4 \
	-O $SNV/${name}.gatkLabeled.vcf \
	--invalidate-previous-filters


	#bcftools view -i 'FILTER="."' $SNV/${name}_cnnLabeled.vcf  > $SNV/${name}_gatkFiltered.vcf
	#s6="$?"


	#--resource ${CNN_REF}/b37_1000G_omni2.5.b37.vcf.gz \


	# # compare 2 vcf files

	# gatk Concordance \
	# -truth gs://gatk-tutorials/$WORKSHOP/2-germline/CNNScoreVariants/vcfs/hg001_na12878_b37_truth.vcf.gz \
	# -eval /home/jupyter-user/CNN/Output/my_2d_filtered.vcf \
	# -L 20:1000000-1432828 \
	# -S /home/jupyter-user/CNN/Output/2d_filtered_concordance.txt




else




	printf "\n\n\n- Hard Filtering (classical way GATK)"
	printf "\n---------------------------------------\n"


	#HARD FILTERING
	#First step extacting the SNP's
	#Second step extracting the INDEL's


	### SNPs

	#1.Extract the SNP's from the call set.
	#printf "\nExtract the SNP's from the call set."
	start=`date +%s`

	gatk SelectVariants --tmp-dir $TMP \
	-R $fasta \
	-V $GEN/${name}.vcf \
	--select-type-to-include SNP \
	-O $SNV/${name}.snp.vcf
	s1="$?"


	#2.Apply the filters to the SNP's callset.

	gatk VariantFiltration --tmp-dir $TMP \
	-V $SNV/${name}.snp.vcf \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
	-O $SNV/${name}.filtered.snp.vcf
	s2="$?"




	### INDELs

	#3. Extract the INDELS from the ORIGINAL call set.
	gatk SelectVariants --tmp-dir $TMP \
	-R $fasta \
	-V $GEN/${name}.vcf \
	--select-type-to-include INDEL \
	-O $SNV/${name}.indel.vcf
	s3="$?"



	#4.Apply the filters to the INDEL's callset.
	gatk VariantFiltration --tmp-dir $TMP \
	-V $SNV/${name}.indel.vcf \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
	-O $SNV/${name}.filtered.indel.vcf
	s4="$?"


	# Combine Variants after using SNPS and INDELS filtering into a single file and get it ready for annotation.

	gatk MergeVcfs --TMP_DIR $TMP \
	-R $REF \
	-I $SNV/${name}.filtered.snp.vcf \
	-I $SNV/${name}.filtered.indel.vcf \
	-O $SNV/${name}.gatkLabeled.vcf
	s5="$?"

	#bcftools view -i 'FILTER="PASS"' $SNV/${name}_hLabeled.vcf  > $SNV/${name}_gatkFiltered.vcf
	#s6="$?"


	if [ "$s1" = "0"  ] &&  [ "$s2" = "0" ] &&  [ "$s3" = "0" ] &&  [ "$s4" = "0" ] &&  [ "$s5" = "0" ]; then
		printf '\nEXIT STATUS: 0'
		printf  '\nHARD FILTERING for '${name}' DONE'

		# rm $GEN/${name}.vcf*
		rm $SNV/${name}.indel.vcf*
		rm $SNV/${name}.snp.vcf*
		rm $SNV/${name}.filtered.indel.vcf*
		rm $SNV/${name}.filtered.snp.vcf*
	else
		printf "\nERROR: PROBLEMS WITH HARD FILTERING"
		exit 1
	fi


	end=`date +%s`
	runtime=$((end-start))
	printf '\nExecuting time: '$runtime 






fi



# Removing temporal forder

rm -r $TMP










# if [ "$cvcf" = "True" ]; then 
	
# 	printf '\n' $SNV/${name}_raw.vcf >> $SNV/filteredVCFs_to_combine_$run.list

# fi

