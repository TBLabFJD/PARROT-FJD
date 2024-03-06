#!/bin/bash 

# CNV results format change: tabulated files to VCF

MDAP=$1
run=$2
fasta=$3
nprograms=0
dt=`date +"%y%m%d"`

dir=$MDAP/cnvs/

for file in $dir/$run*combined*txt; do

	echo $file
	combined_file=$file
	combined_file_tmp=${combined_file}".tmp"
	combined_file_noheader=${combined_file}".noheader"

	outputpos_tmp=${combined_file}".pos"
	outputfile=`basename $file .txt`


	awk -F '\t' '{if($1!~/^#/ && $14>'$nprograms'){print $5"\t"$6-1"\t"$6}}'  ${combined_file} |sed 1d | sed 's/^/chr/g'  > ${combined_file_tmp}
	awk -F '\t' '{if($1!~/^#/ && $14>'$nprograms'){print $0}}'  ${combined_file} | sed 1d > ${combined_file_noheader}

	bedtools getfasta -fi $fasta -bed $combined_file_tmp  -tab -fo $outputpos_tmp



	echo '##fileformat=VCFv4.3' > $dir/${outputfile}.vcf
	echo "##fileDate=$dt"  >> $dir/${outputfile}.vcf
	echo "##reference=file:$fasta" >> $dir/${outputfile}.vcf
	echo '##INFO=<ID=BKPTID,Number=.,Type=String,Description="ID of the assembled alternate allele in the assembly file">' >> $dir/${outputfile}.vcf
	echo '##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">' >> $dir/${outputfile}.vcf
	echo '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">' >> $dir/${outputfile}.vcf
	echo '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record"> ' >> $dir/${outputfile}.vcf
	echo '##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">' >> $dir/${outputfile}.vcf
	echo '##INFO=<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints">' >> $dir/${outputfile}.vcf
	echo '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">' >> $dir/${outputfile}.vcf
	echo '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">' >> $dir/${outputfile}.vcf
	echo '##ALT=<ID=DEL,Description="Deletion">' >> $dir/${outputfile}.vcf
	echo '##ALT=<ID=DEL:ME:ALU,Description="Deletion of ALU element">' >> $dir/${outputfile}.vcf
	echo '##ALT=<ID=DEL:ME:L1,Description="Deletion of L1 element">' >> $dir/${outputfile}.vcf
	echo '##ALT=<ID=DUP,Description="Duplication">' >> $dir/${outputfile}.vcf
	echo '##ALT=<ID=DUP:TANDEM,Description="Tandem Duplication">' >> $dir/${outputfile}.vcf
	echo '##ALT=<ID=INS,Description="Insertion of novel sequence">' >> $dir/${outputfile}.vcf
	echo '##ALT=<ID=INS:ME:ALU,Description="Insertion of ALU element">' >> $dir/${outputfile}.vcf
	echo '##ALT=<ID=INS:ME:L1,Description="Insertion of L1 element">' >> $dir/${outputfile}.vcf
	echo '##ALT=<ID=INV,Description="Inversion">' >> $dir/${outputfile}.vcf
	echo '##ALT=<ID=CNV,Description="Copy number variable region">' >> $dir/${outputfile}.vcf
	echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' >> $dir/${outputfile}.vcf
	echo '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">' >> $dir/${outputfile}.vcf
	echo '##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">' >> $dir/${outputfile}.vcf
	echo '##FORMAT=<ID=CNQ,Number=1,Type=Float,Description="Copy number genotype quality for imprecise events">' >> $dir/${outputfile}.vcf
	echo '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT' >> $dir/${outputfile}.vcf

	awk -F '\t' 'NR==FNR { a[FNR]=$2;next } {if($4=="HET"){print $5"\t"$6"\t.\t"a[FNR]"\t<"$3">\t.\tPASS\tSVTYPE="$3";END="$7"\tGT\t0/1"}else{print $5"\t"$6"\t.\t"a[FNR]"\t<"$3">\t.\tPASS\tSVTYPE="$3";END="$7"\tGT\t1/1"}}' $outputpos_tmp $combined_file_noheader | sort | uniq >> $dir/${outputfile}.vcf

	rm ${combined_file_tmp}
	rm ${combined_file_noheader}
	rm ${outputpos_tmp}

done



	

