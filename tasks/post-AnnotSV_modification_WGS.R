



cnv_path = "/home/gonzalo/tblab/mnt/tblab/gonzalo/reanalysis/aDamian_151122_WGS2022_CNVs/cnvs/2022-11-16_11-04.CNV.annotated.tsv"
cnv = read.delim(cnv_path, header = TRUE, stringsAsFactors = F, quote = "", check.names=F, colClasses = "character")


samples_col = which(colnames(cnv) %in% c("FORMAT", "Annotation_mode")) + c(+1, -1)
samples_col = colnames(cnv)[samples_col[1]:samples_col[2]]

for (sample in samples_col){
  cnv[,paste0(sample, "_GT")] = gsub(":.*", "", cnv[,sample])
}
new_col = paste0(samples_col, "_GT")

cnv$N_samples = unlist(apply(cnv[,new_col],1, function(x){
  length(grep("*/1", x))
}))


cnv_filt = cnv[cnv$N_samples < 3,]

output="/home/gonzalo/tblab/mnt/tblab/gonzalo/reanalysis/aDamian_151122_WGS2022_CNVs/cnvs/2022-11-16_11-04.CNV.annotated.filt.tsv"
write.table(cnv_filt, output, sep = "\t", col.names = T, row.names = F, quote = F, na = "")

library(openxlsx)
write.xlsx(cnv_filt, paste0(output, ".xlsx"), colNames = TRUE)


