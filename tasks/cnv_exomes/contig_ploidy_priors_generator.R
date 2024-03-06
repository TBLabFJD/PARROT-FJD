

# cut -f1 hg38.fa.fai > /home/gonzalo/tblab/home/gonzalo/nextflowtest/NextVariantFJD/tasks/cnv_exomes/contig_ploidy_priors.hg38.tsv
# cut -f1 hg19.fa.fai > /home/gonzalo/tblab/home/gonzalo/nextflowtest/NextVariantFJD/tasks/cnv_exomes/contig_ploidy_priors.hg19.tsv

path_df="/home/gonzalo/tblab/home/gonzalo/nextflowtest/NextVariantFJD/tasks/cnv_exomes/contig_ploidy_priors.hg38.tsv"
path_df="/home/gonzalo/tblab/home/gonzalo/nextflowtest/NextVariantFJD/tasks/cnv_exomes/contig_ploidy_priors.hg19.tsv"


contig_df = read.delim(path_df, header = F, stringsAsFactors = F, quote = "")[,1,drop=F]

  
colnames(contig_df) = "CONTIG_NAME"
contig_df$PLOIDY_PRIOR_0 = "0.0"
contig_df$PLOIDY_PRIOR_1 = "0.01"
contig_df$PLOIDY_PRIOR_2 = "0.98"
contig_df$PLOIDY_PRIOR_3 = "0.01"

contig_df$PLOIDY_PRIOR_0[grepl("chrX", contig_df$CONTIG_NAME)] = "0.01"
contig_df$PLOIDY_PRIOR_1[grepl("chrX", contig_df$CONTIG_NAME)] = "0.49"
contig_df$PLOIDY_PRIOR_2[grepl("chrX", contig_df$CONTIG_NAME)] = "0.49"
contig_df$PLOIDY_PRIOR_3[grepl("chrX", contig_df$CONTIG_NAME)] = "0.01"

contig_df$PLOIDY_PRIOR_0[grepl("chrY", contig_df$CONTIG_NAME)] = "0.495"
contig_df$PLOIDY_PRIOR_1[grepl("chrY", contig_df$CONTIG_NAME)] = "0.495"
contig_df$PLOIDY_PRIOR_2[grepl("chrY", contig_df$CONTIG_NAME)] = "0.01"
contig_df$PLOIDY_PRIOR_3[grepl("chrY", contig_df$CONTIG_NAME)] = "0.0"


write.table(contig_df, path_df, col.names = T, row.names = F, quote = F, sep = "\t")



