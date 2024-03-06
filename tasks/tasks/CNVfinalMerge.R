#### FINAL CNV FILE ####

library(dplyr)
library(optparse)


option_list=list(
  make_option(c('-d','--resultsdir'),type="character", help="Results directory."),
  make_option(c('-r','--run'),type="character", help="Run name."))

opt <- parse_args(OptionParser(option_list = option_list))


cnv_dir=paste(opt$resultsdir, "/cnvs/", sep='')
mixer=paste(cnv_dir, opt$run, ".combined.txt", sep='')
annotsv=paste(cnv_dir, opt$run, ".combinedAnnotated.tsv", sep='')


df_mixer <- read.csv( mixer, fill=TRUE, sep= '\t', skip= 13)
print("MIXER... READ")

df_annot <- read.csv( annotsv, fill=TRUE, sep= '\t' )
print("ANNOTSV OUTPUT... READ")
 
# df_mixer_to_annot <- df_mixer[,c("SAMPLE", "HOMO_HETEROZYGOUS", "CHR", "START", "STOP", "NUM_OF_PROGRAMS" )]

drop_columns <- c("AnnotSV.ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "TADcoordinates", "ENCODEexperiments",  "GHid_not_elite", "ACMG" )

df_annot_corrected <-  df_annot[, ! names(df_annot) %in% drop_columns, drop = FALSE]
# , 

df_intermediate <- left_join( df_mixer, df_annot_corrected, 
	by=c("CHR"="SV.chrom","START"="SV.start", "STOP"="SV.end", "CNV_TYPE"="SV.type" )) %>% distinct()

df_final <- df_intermediate[, c(89, 1:7, 14, 8:13, 23, 25:88, 15:22)]

final_output= paste(cnv_dir,  opt$run, ".final.txt", sep="")

write.table(df_final, final_output, sep="\t", row.names=FALSE,  quote= FALSE )
print("CNV annotation is Finished")

