### CNV analysis 
### Author: Gonzalo Núñez Moreno
### Date: 22/03/19


rm(list=ls())


#**********************#
# package requirements #
#**********************#

library(optparse)
library(bedr)
library(plyr)

#==============#
#   Arguments  #
#==============#
option_list=list(
  make_option(c('-i','--inputdir'),type="character", help="Input directory."),
  make_option(c('-o','--outputfile'),type="character", help="Output directory."),
  make_option(c('-s','--samples'), type="character", default=NULL, help="Samples to keep")
)
  
opt_parser=OptionParser(option_list = option_list)
opt=parse_args(opt_parser) #list of the args


# opt=list()
# opt$inputdir <- "/home/gonzalo/Documents/prueba_merge_cnv_tmp"
# opt$outputfile <- "/home/gonzalo/Documents/prueba_merge_cnv_tmp/PRUEBA.CNV.merged.bed"
# opt$samples <- "/home/gonzalo/Documents/prueba_merge_cnv_tmp/samples2analyce.txt"

# opt=list()
# opt$inputdir <- "/home/gonzalo/tblab/mnt/genetica4/gonzalo/pruebas_merge_cnvs/work/49/7e9f779328cc16dff25a77f2a1c9e8/"
# opt$outputfile <- "/home/gonzalo/tblab/mnt/genetica4/gonzalo/pruebas_merge_cnvs/work/49/7e9f779328cc16dff25a77f2a1c9e8/PRUEBA.CNV.merged.bed"
# opt$samples <- NULL


inputdir <- opt$inputdir
outputfile <- opt$outputfile
samplesfilt <- opt$samples
  
setwd(inputdir)




#==================#
#   Input loading  #
#==================#
toAnnotate_files = list.files(inputdir, pattern = ".toAnnotate.txt")

programs = c()
extra_col_names = c()
extra_info_list = list()
all_progrmas_df = data.frame(stringsAsFactors = F)
for (toAnnotate_file in toAnnotate_files) {
  toAnnotate_df = read.delim(toAnnotate_file, header = T, stringsAsFactors = F, quote = "")
  cnv_program = gsub(".toAnnotate.txt","",toAnnotate_file)
  programs = c(programs,cnv_program)
  # Add "chr" prefix if it is not present
  toAnnotate_df[grep("chr", toAnnotate_df[,1],invert = T), 1] = paste0("chr", toAnnotate_df[grep("chr", toAnnotate_df[,1],invert = T),1])
  if (ncol(toAnnotate_df) == 5) { 
    colnames(toAnnotate_df) = c("CHR", "START", "END", "CNV_TYPE", "SAMPLE")
    toAnnotate_df$PROGRAM = cnv_program
    
    extra_info_df = data.frame(
      paste0(toAnnotate_df$CHR, ":", toAnnotate_df$START, "-", toAnnotate_df$END, "_", toAnnotate_df$CNV_TYPE),
      paste(toAnnotate_df$CHR, toAnnotate_df$START, toAnnotate_df$END, toAnnotate_df$CNV_TYPE, toAnnotate_df$SAMPLE, toAnnotate_df$PROGRAM, sep = "_")
    )
    colnames(extra_info_df) = c(cnv_program, "cnv_id")
    extra_info_list[[cnv_program]] = extra_info_df
    
  } else { 
    colnames(toAnnotate_df) = c("CHR", "START", "END", "CNV_TYPE", "SAMPLE", paste0(cnv_program, "_", colnames(toAnnotate_df)[6:ncol(toAnnotate_df)]) )
    toAnnotate_df$PROGRAM = cnv_program
    
    extra_col_names = c(extra_col_names, colnames(toAnnotate_df)[6:(ncol(toAnnotate_df)-1)])
    
    extra_info_df = toAnnotate_df[,6:(ncol(toAnnotate_df)-1), drop = F]
    extra_info_df[,cnv_program] = paste0(toAnnotate_df$CHR, ":", toAnnotate_df$START, "-", toAnnotate_df$END, "_", toAnnotate_df$CNV_TYPE)
    extra_info_df[,"cnv_id"] = paste(toAnnotate_df$CHR, toAnnotate_df$START, toAnnotate_df$END, toAnnotate_df$CNV_TYPE, toAnnotate_df$SAMPLE, toAnnotate_df$PROGRAM, sep = "_")
    extra_info_list[[cnv_program]] = extra_info_df
    
    toAnnotate_df = toAnnotate_df[,c("CHR", "START", "END", "CNV_TYPE", "SAMPLE", "PROGRAM")]
  }
  toAnnotate_df$COMBINED = paste(toAnnotate_df$CHR, toAnnotate_df$START, toAnnotate_df$END, toAnnotate_df$CNV_TYPE, toAnnotate_df$SAMPLE, toAnnotate_df$PROGRAM, sep = "_")
  toAnnotate_df$CNV_SIMPLE = gsub("HOM_", "", toupper(toAnnotate_df$CNV_TYPE))
  all_progrmas_df = rbind(all_progrmas_df, toAnnotate_df)
}

extra_info_df = rbind.fill(extra_info_list)
extra_info_df = extra_info_df[!duplicated(extra_info_df$cnv_id),]
rownames(extra_info_df) = extra_info_df$cnv_id



if (!is.null(samplesfilt)){
  samplesfilt_df = read.delim(samplesfilt, stringsAsFactors = F, header = F, quote = "")
  samplesfilt_df$V1 = gsub("_.*", "", samplesfilt_df$V1)
  all_progrmas_df = all_progrmas_df[grep(paste(samplesfilt_df$V1, collapse = "|"), all_progrmas_df$SAMPLE),]
}






all_progrmas_df_sampleSplit = split(all_progrmas_df, paste(all_progrmas_df$SAMPLE, all_progrmas_df$CNV_SIMPLE, sep = "__"))

merged_cnv_df = data.frame()

for (sample_cnvtype in all_progrmas_df_sampleSplit){
  sample_cnvtype = sample_cnvtype[order(sample_cnvtype$CHR, sample_cnvtype$START, sample_cnvtype$END),]
  merged_df <- bedr(
    input = list(i = sample_cnvtype), 
    method = "merge", 
    params = "-c 7,8,5 -o distinct,distinct,distinct -delim ';' ",
    verbose = F
  )
  
  merged_cnv_df = rbind(merged_cnv_df, merged_df)
}





extra_info_collapse_df = data.frame()
extra_info_collapse_list = list()
cnv_ids = strsplit(as.character(merged_cnv_df$V4), ";")
j = 1
for (i in cnv_ids){
  selected_rows = extra_info_df[i,]
  aa <- apply( selected_rows, 2 , function(x) {paste(x[!is.na(x)], collapse = ";")} )
  extra_info_collapse_list[[j]] = aa
  j = j + 1
}

extra_info_collapse_df = do.call("rbind", extra_info_collapse_list)

final_df = cbind(merged_cnv_df[,c(1:3,5,6)], extra_info_collapse_df[,c(programs, extra_col_names)])
colnames(final_df)[1:5] = c("#CHR", "START", "END", "CNV_TYPE", "SAMPLE")


final_df$N_PROGRAMS = unlist(apply(final_df[,programs,drop=F], 1, function(x) sum(!x == "")))

final_df = final_df[,c("#CHR", "START", "END", "CNV_TYPE", "SAMPLE", "N_PROGRAMS", programs, extra_col_names)]
final_df = final_df[order(final_df$N_PROGRAMS, decreasing = T),]
final_df$`#CHR` = gsub("chr", "", final_df$`#CHR`)


write.table(final_df, outputfile, col.names = T, row.names = F, sep = "\t", quote = F)

write.table(c(programs, extra_col_names), "colnames.txt", col.names = F, row.names = F, sep = "\t", quote = F)


