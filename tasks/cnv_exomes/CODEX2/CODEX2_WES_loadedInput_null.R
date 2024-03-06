### CNV analysis 
### Author: Lorena de la Fuente
### Date: 07/07/20


rm(list=ls())

#start clock
ptm <- proc.time()
print(ptm)


#**********************#
# package requirements #
#**********************#

library(CODEX2)
library(optparse)


#**************#
#   Arguments  #
#**************#
# ¡¡¡¡¡¡¡¡Samples from different batches are highly recommended to run separately!!!!!!!!!!!
# ¡¡¡¡¡¡¡¡Data sets of sample size ranging from 30 to 500!!!!!!!!!!!
option_list=list(
  make_option(c("-d","--dir"),type="character", help="Directory with input bam files to analyse."),
  make_option(c('-o','--outputdir'),type="character", help="Output directory."),
  make_option(c('-n','--name'),type="character", help="Project name."),
  make_option(c('-b','--bed'),type="character",default="genome", help="Bed file with genes found in panel"))

opt_parser=OptionParser(option_list = option_list)
opt=parse_args(opt_parser) #list of the args

bamdir <- opt$dir
output_dir <- paste(opt$outputdir,"/copy_number_variation_data/CODEX2/",sep="")
bedFile <- opt$bed
projectname <- opt$name

# bamdir <- "/mnt/genetica3/Ionut/MCorton/SureSelect_Glaucoma_03/all_bams_new/Sure_123/some"
# output_dir <- "/mnt/genetica/gonzalo/Codex2/"
# bedFile <- "/mnt/genetica/gonzalo/CoNIFER/muestras_Marta/SureSelect_Glaucoma_2018_real_target_manifest_v2_numbered.bed"
# projectname <- "some"

sink(paste(opt$outputdir,"/software_",opt$name,".txt",sep=""),append=TRUE)
print("R SESSION INFO (CODEX2):")
sessionInfo()
sink()

dir.create(output_dir)


#****************#
# Initialization #
#****************#
bamFile <- list.files(bamdir, pattern = '*.bam$')
bamdir <- file.path(bamdir, bamFile)
sampname <- sub("_.*$","",bamFile)

bambedObj <- getbambed(bamdir = bamdir, bedFile = bedFile, sampname = sampname, projectname = projectname)
bamdir <- bambedObj$bamdir; sampname <- bambedObj$sampname
ref <- bambedObj$ref; projectname <- bambedObj$projectname

setwd(output_dir)


#************************************#
# Getting GC content and mappability #
#************************************#
#gc <- getgc(ref)
#mapp <- getmapp(ref)
#values(ref) <- cbind(values(ref), DataFrame(gc, mapp))  



# #****************************************************#
# # Getting gene names, needed for targeted sequencing #
# #****************************************************#
# bedFile_data <- read.delim(bedFile, header = FALSE, stringsAsFactors = FALSE)
# bedFile_data[,1] <- gsub("chr", "", bedFile_data[,1])
# bedFile_data <- bedFile_data[order(bedFile_data[,1], bedFile_data[,2]),]
# gene <- bedFile_data[,4]
# values(ref) <- cbind(values(ref), DataFrame(gc, mapp, gene))




#***************************#
# Getting depth of coverage #
#***************************#

# Read depth matrix, as well as read lengths across all samples, will be returned. This will need to be generated for each chromosome.

#coverageObj <- getcoverage(bambedObj, mapqthres = 20)
#Y <- coverageObj$Y

load(file = paste('CODEX2',"_counts_image.RData", sep=''))

#load(paste('CODEX2',"_counts_image.RData", sep=''))




#*****************#
# Quality control #
#*****************#
# qcObj <- qc(Y, sampname, ref, cov_thresh = c(20, 10000),
#              length_thresh = c(20, 10000), mapp_thresh = 0.9,
#              gc_thresh = c(20, 80))
# Y_qc <- qcObj$Y_qc
# sampname_qc <- qcObj$sampname_qc
# ref_qc <- qcObj$ref_qc; 
# qcmat <- qcObj$qcmat
# gc_qc <- ref_qc$gc
# # write.table(qcmat, file = paste(projectname, '_qcmat', '.txt', sep=''),
# #             sep = '\t', quote = FALSE, row.names = FALSE)

# save(qcObj, file = paste('CODEX2',"_qc.RData", sep=''))
load(file = paste('CODEX2',"_norm.RData", sep=''))


#************************************************#
# Estimating library size factor for each sample #
#************************************************#
Y.nonzero <- Y_qc[apply(Y_qc, 1, function(x){!any(x==0)}),]
pseudo.sample <- apply(Y.nonzero,1,function(x){prod(x)^(1/length(x))})
N <- apply(apply(Y.nonzero, 2, function(x){x/pseudo.sample}), 2, median)





############# Running CODEX2 without specifying negative control samples/regions




#************************************************#
# Genome-wide normalization using normalize_null #
#************************************************#

# normObj.null <- normalize_null(Y_qc = Y_qc, gc_qc = gc_qc, K = 1:5, N=N)

Yhat.null <- normObj.null$Yhat
AIC.null <- normObj.null$AIC
BIC.null <- normObj.null$BIC
RSS.null <- normObj.null$RSS

# z.codex <- log(Y_qc/normObj.null$Yhat[[which.max(normObj.null$BIC)]])
# cnv_index1 <- which(apply(z.codex,1,sd)>=0.25) 

# print(cnv_index1)
save(normObj.null, file = paste('CODEX2',"_norm_null.RData", sep=''))


# normObj <- normalize_codex2_nr(Y_qc = Y_qc, gc_qc = gc_qc, cnv_index = cnv_index1, K = 1:5, N = N)
# Yhat <- normObj$Yhat
# AIC <- normObj$AIC; 
# BIC <- normObj$BIC
# RSS <- normObj$RSS

# save(normObj file = paste('CODEX2',"_norm_nr.RData", sep=''))




#**********************************************************#
# CBS segmentation per chromosome: optimal for WGS and WES #
#**********************************************************#

finalcall.CBS = NULL 

for (chr in unique(seqnames(ref_qc))){
  
  chr.index <- which(seqnames(ref_qc)==chr)
  Yhat.null.chr <- lapply(Yhat.null, function(x) x[chr.index,])  

  finalcall.CBS[[chr]] <- segmentCBS(Y_qc = Y_qc[chr.index,],  # recommended
                              Yhat = Yhat.null.chr, optK = which.max(BIC.null),
                              K = 1:5,
                              sampname_qc = colnames(Y_qc),
                              ref_qc = ranges(ref_qc)[chr.index],
                              chr=chr,
                              lmax = 400, mode = "integer")
}


save(finalcall.CBS, file = paste('CODEX2',"_CNV_null.RData", sep=''))









# optK=which.max(BIC)
# finalcall=matrix(ncol=14,nrow=0)
# colnames(finalcall)=c('sample_name','chr','gene','cnv',
#                       'st_bp','ed_bp','length_kb',
#                       'st_exon','ed_exon','raw_cov',
#                       'norm_cov','copy_no','lratio',
#                       'mBIC')
# for(genei in unique(ref_qc$gene)){
#   cat('Segmenting gene',genei,'\n')
#   geneindex=which(ref_qc$gene==genei)
#   yi=Y_qc[geneindex,, drop=FALSE]
#   yhati=Yhat[[optK]][geneindex,, drop=FALSE]
#   refi=ref_qc[geneindex]
#   finalcalli=segment_targeted(yi, yhati, sampname_qc, refi, genei, lmax=length(geneindex), mode='fraction') 
#   finalcall=rbind(finalcall,finalcalli)
# }
# cn=(as.numeric(as.matrix(finalcall[,'copy_no'])))
# cn.filter=(cn<=1.7)|(cn>=2.3) # removing calls with fractional copy numbers close to 2 (for heterogeneous cancer samples)
# finalcall=finalcall[cn.filter,]
# length_exon=as.numeric(finalcall[,'ed_exon'])-as.numeric(finalcall[,'st_exon'])+1
# finalcall=cbind(finalcall[,1:7],length_exon,finalcall[,10:14])

# finalcall[,"cnv"] <- toupper(finalcall[,"cnv"]) # Transform del/dup to DEL/DUP

# write.table(finalcall, file = paste( 'CODEX2', '_results.txt', sep=''), sep='\t', quote=F, row.names=F)

# toAnnotateTable <- data.frame(finalcall[,"chr"], finalcall[,"st_bp"], finalcall[,"ed_bp"], toupper(finalcall[,"cnv"]))

# write.table(toAnnotateTable, file = paste( 'CODEX2','_toAnnotate.txt', sep=''), sep='\t', quote=F, row.names=F, col.names = FALSE)


# #stop clock
# finish <- proc.time() - ptm
# print(finish)
# summary(finish)
