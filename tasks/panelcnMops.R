### CNV analysis 
### Author: Gonzalo Núñez Moreno
### Date: 20/03/19


rm(list=ls())

#start clock
ptm <- proc.time()
print(ptm)


#**********************#
# package requirements #
#**********************#


library(panelcn.mops)
library(optparse)


#**************#
#   Arguments  #
#**************#
option_list=list(
  make_option(c("-d","--dir"),type="character", help="Directory with input bam files to analyse."),
  make_option(c('-o','--outputdir'),type="character", help="Output directory."),
  make_option(c('-n','--name'),type="character", help="Project name."),
  make_option(c('-b','--bed'),type="character",default="genome", help="Bed file with genes found in panel"))

opt_parser=OptionParser(option_list = option_list)
opt=parse_args(opt_parser) #list of the args

dirPath <- opt$dir
# output_dir <- paste(opt$outputdir,"/cnvs/Panelcn.MOPS/",sep="")
bedFile <- opt$bed
projectname <- opt$name

# dirPath <- "/mnt/genetica3/Ionut/MCorton/SureSelect_Glaucoma_03/all_bams_new/Sure_123/some"
# output_dir <- "/mnt/genetica/gonzalo/panelcn.mops/"
# bedFile <- "/mnt/genetica/gonzalo/CoNIFER/muestras_Marta/SureSelect_Glaucoma_2018_real_target_manifest_v2_numbered.bed"
# projectname <- "some"

sink(paste("software_",opt$name,".txt",sep=""),append=TRUE)
print("R SESSION INFO (Panelcn.MOPS):")
sessionInfo()
sink()

# dir.create(output_dir)


#*****************************************#
# Getting read counts (RCs) from BAM file #
#*****************************************#
countWindows <- getWindows(bedFile) # Getting count windows from the BED file 
bamFile <- list.files(dirPath, pattern = '*.bam$')
setwd(dirPath)
test <- countBamListInGRanges(countWindows = countWindows,
                              bam.files = bamFile, 
                              read.width = 150)
# setwd(output_dir)
save.image(file = paste('panelcn.mops.','_counts_image.RData', sep=''))
write.table(data.frame(test), file = paste('panelcn.mops','_count_matrix.txt', sep=''), sep = "\t", quote = F, row.names = F, col.names = T)
# write.table(data.frame(test), file = '/home/gonzalo/Documents/cnvs/panelcn.mops_count_matrix.txt', sep = "\t", quote = F, row.names = F, col.names = T)

#***********************#
# Running the algorithm #
#***********************#
vec_pos <- 1:(ncol(elementMetadata(test)))
temporal_test <- test
finaltable <- data.frame()
for (i in 1:(ncol(elementMetadata(test)))){
  elementMetadata(temporal_test) <- cbind(elementMetadata(test)[vec_pos[i]],
                                          elementMetadata(test)[vec_pos[-i]])

  resultlist <- runPanelcnMops(XandCB = temporal_test,
                               testiv = 1,
                               countWindows = countWindows)
  
  sampleNames <- colnames(elementMetadata(temporal_test))
  resulttable <- createResultTable(resultlist = resultlist, 
                                   XandCB = temporal_test, 
                                   countWindows = countWindows, 
                                   sampleNames = sampleNames)
  resulttable <- resulttable[[1]]
  resulttable <- resulttable[!grepl("CN2",resulttable$CN),] # Delete rows that does not have a CNV
  finaltable <- rbind(finaltable,resulttable)
}

finaltable$Sample <- sub("_.*$","",finaltable$Sample) # Change samplename

finaltable$CN <- sub("CN0","HOM_DEL",finaltable$CN)
finaltable$CN <- sub("CN1","DEL",finaltable$CN)
finaltable$CN <- sub("CN3","DUP",finaltable$CN)
finaltable$CN <- sub("CN4","HOM_DUP",finaltable$CN)


finaltable <- finaltable[!grepl("lowQual",finaltable$lowQual),] # Delete rows with low quality

write.table(finaltable, file = 'panelcn.MOPS.results.txt', sep='\t', quote=F, row.names=F)

toAnnotateTable <- data.frame(finaltable$Chr, finaltable$Start, finaltable$End, finaltable$CN, finaltable$Sample, finaltable$RC.norm/finaltable$medRC.norm)
colnames(toAnnotateTable) = c("CHR", "START", "END", "CNV_TYPE", "SAMPLE","RATIO")
toAnnotateTable$SAMPLE = gsub("\\..*","", toAnnotateTable$SAMPLE, perl = TRUE)
toAnnotateTable$CHR[grep("chr", toAnnotateTable$CHR,invert = T)] = paste0("chr", toAnnotateTable$CHR[grep("chr", toAnnotateTable$CHR,invert = T)])
# toAnnotateTable = toAnnotateTable[toAnnotateTable$RATIO > 1.46 | toAnnotateTable$RATIO < 0.57, ]

write.table(toAnnotateTable, file ='panelcn.MOPS.toAnnotate.txt', sep='\t', quote=F, row.names=F, col.names = T)

#******************************************************#  
# Función que se le da directamete la matriz de conteo #  ¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡ CAMBIAR !!!!!!!!!!!!!!!!!!!!!!
#******************************************************# 
# matrix_count <- Y
# countWindows <- getWindows(bedFile) # Getting count windows from the BED file
# vec_pos <- 1:(ncol(matrix_count))
# temporal_test <- matrix()
# finaltable <- data.frame()
# for (i in 1: ncol(matrix_count)){
#   temporal_test <- cbind(matrix_count[,vec_pos[i]],
#                          matrix_count[,vec_pos[-i]])
#   resultlist <- panelcn.mops(input = matrix_count,
#                         testi = 1)
#   sampleNames <- colnames(matrix_count)
#   resulttable <- createResultTable(resultlist = resultlist,
#                                    XandCB = matrix_count,
#                                    countWindows = countWindows,
#                                    sampleNames = sampleNames)
#   resulttable <- resulttable[[1]]
#   resulttable <- resulttable[!grepl("CN2",resulttable$CN),] # Delete rows that does not have a CNV
#   finaltable <- rbind(finaltable,resulttable)
# }
# 
# finaltable$Sample <- sub("_.*$","",finaltable$Sample)
# 
# finaltable$CN <- sub("CN1|CN0","DEL",finaltable$CN)
# finaltable$CN <- sub("CN3|CN4","DUP",finaltable$CN)
# 
# write.table(finaltable, file = paste(projectname, '_panelcn.MOPS.txt', sep=''), sep='\t', quote=F, row.names=F)
# 
# toAnnotateTable <- data.frame(finaltable$Chr, finaltable$Start, finaltable$End, finaltable$CN)
# 
# write.table(toAnnotateTable, file = paste(projectname, '_panelcn.MOPS_toAnnotate.txt', sep=''), sep='\t', quote=F, row.names=F)


#stop clock
finish <- proc.time() - ptm
print(finish)
summary(finish)
