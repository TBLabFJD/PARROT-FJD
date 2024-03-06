module load bedtools

panel=$1
gff=$2


#awk '{if($3=="exon"){print $0}}' ${gff} > ${gff}_EXOMA.gtf
bedtools intersect -a ${panel} -b ${gff}_EXOMA.gtf -wao  > ${panel}_intersection.bed

awk '{print $1"\t"$2"\t"$3"\t"$NF"\t"$13}'  ${panel}_intersection.bed | sed 's/"//g' |  sed 's/;//g' | sort -k1,1 -k2,2n -k4,4nr | awk '!seen[$1, $2, $3]++' > ${panel}_intersection_uniq.bed

awk '{if($5!="."){print $1"\t"$2"\t"$3"\t"$5}}' ${panel}_intersection_uniq.bed > ${panel}_annotated.bed

rm ${panel}_intersection.bed
rm ${panel}_intersection_uniq.bed
