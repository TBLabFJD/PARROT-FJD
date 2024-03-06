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
output_dir <- paste(opt$outputdir,"/cnvs/CODEX2/",sep="")
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
sampname <- sub("*.bam$","",bamFile)

bambedObj <- getbambed(bamdir = bamdir, bedFile = bedFile, sampname = sampname, projectname = projectname)
bamdir <- bambedObj$bamdir; sampname <- bambedObj$sampname
ref <- bambedObj$ref; projectname <- bambedObj$projectname

setwd(output_dir)


#************************************#
# Getting GC content and mappability #
#************************************#
gc <- getgc(ref)
mapp <- getmapp(ref)


#****************************************************#
# Getting gene names, needed for targeted sequencing #
#****************************************************#
bedFile_data <- read.delim(bedFile, header = FALSE, stringsAsFactors = FALSE)
bedFile_data[,1] <- gsub("chr", "", bedFile_data[,1])
bedFile_data <- bedFile_data[order(bedFile_data[,1], bedFile_data[,2]),]
gene <- bedFile_data[,4]
values(ref) <- cbind(values(ref), DataFrame(gc, mapp, gene))


#***************************#
# Getting depth of coverage #
#***************************#
coverageObj <- getcoverage(bambedObj, mapqthres = 20)
Y <- coverageObj$Y

save.image(file = 'CODEX2_counts_image.RData')
write.table(Y, file = 'CODEX2_count_matrix.txt', sep = "\t", quote = F, row.names = T, col.names = T)
# write.table(Y, file = '/home/gonzalo/Documents/cnvs/CODEX2_count_matrix.txt', sep = "\t", quote = F, row.names = T, col.names = T)

#*****************#
# Quality control #
#*****************#
qcObj <- qc(Y, sampname, ref, cov_thresh = c(20, Inf),
            length_thresh = c(20, Inf), mapp_thresh = 0.9,
            gc_thresh = c(20, 80))
Y_qc <- qcObj$Y_qc; sampname_qc <- qcObj$sampname_qc
ref_qc <- qcObj$ref_qc; qcmat <- qcObj$qcmat; gc_qc <- ref_qc$gc


#************************************************#
# Estimating library size factor for each sample #
#************************************************#
Y.nonzero <- Y_qc[apply(Y_qc, 1, function(x){!any(x==0)}),]
pseudo.sample <- apply(Y.nonzero,1,function(x){prod(x)^(1/length(x))})
N <- apply(apply(Y.nonzero, 2, function(x){x/pseudo.sample}), 2, median)


#************************************************#
# Genome-wide normalization using normalize_null #
#************************************************#
normObj.null <- normalize_null(Y_qc = Y_qc, gc_qc = gc_qc, K = 1:4, N=N)
# Yhat <- normObj.null$Yhat
# AIC <- normObj.null$AIC; 
# BIC <- normObj.null$BIC
# RSS <- normObj.null$RSS
z.codex <- log(Y_qc/normObj.null$Yhat[[which.max(normObj.null$BIC)]])
cnv_index1 <- which(apply(z.codex,1,sd)>=0.25) 

normObj <- normalize_codex2_nr(Y_qc = Y_qc, gc_qc = gc_qc, cnv_index = cnv_index1, K = 1:4, N = N)
Yhat <- normObj$Yhat
AIC <- normObj$AIC; 
BIC <- normObj$BIC
RSS <- normObj$RSS


#**********************************************************#
# CBS segmentation per chromosome: optimal for WGS and WES #
#**********************************************************#
# chr='chr20'
# chr.index=which(seqnames(ref_qc)==chr)
# finalcall.CBS <- segmentCBS(Y_qc[chr.index,],
#                             Yhat, optK = which.max(BIC),
#                             K = 1:5,
#                             sampname_qc = sampname_qc,
#                             ref_qc = ranges(ref_qc)[chr.index],
#                             chr = chr, lmax = 400, mode = "integer")


#******************************************************#
# CBS segmentation per gene: optinmal for targeted seq #
#******************************************************#

  #=====================================================================#
  #                      Function for segmentation                      #
  #=====================================================================#

  segment_targeted <- function(yi, yhati, sampname_qc, refi, genei, lmax, mode) {
    chri <- as.matrix(seqnames(refi))[1]
    finalcall <- matrix(ncol = 10)
    lmax <- max(1, lmax - 1)
    if(is.vector(yi)){
      yi=t(yi)
      yhati=t(yhati)
    }
    for (sampno in 1:ncol(yi)) {
      #message("Segmenting sample ", sampno, ": ", sampname_qc[sampno], ".")
      y <- yi[, sampno]
      yhat <- yhati[, sampno]
      num <- length(y)
      y <- c(y, rep(0, lmax))
      yhat <- c(yhat, rep(0, lmax))
      i <- rep(1:num, rep(lmax, num))
      j <- rep(1:lmax, num) + i
      yact <- rep(0, length(i))
      lambda <- rep(0, length(i))
      for (k in 1:num) {
        yact[(lmax * k - (lmax - 1)):(lmax * k)] <- cumsum(y[k:(k + 
                                                                  lmax)])[-1]
        lambda[(lmax * k - (lmax - 1)):(lmax * k)] <- cumsum(yhat[k:(k + 
                                                                       lmax)])[-1]
      }
      i <- i[j <= num]
      yact <- yact[j <= num]
      lambda <- lambda[j <= num]
      j <- j[j <= num]
      yact[lambda<20]=20
      lambda[lambda<20]=20
      if (mode == "integer") {
        chat <- round(2 * (yact/lambda))
      } else if (mode == "fraction") {
        chat <- 2 * (yact/lambda)
      }
      lratio <- (1 - chat/2) * lambda + log((chat + 1e-04)/2.0001) * yact
      # chat[chat > 5] <- 5
      if (sum(lratio > 0) > 0) {
        if (sum(lratio > 0) >= 2) {
          finalmat <- (cbind(i, j, yact, lambda, chat, lratio))[lratio > 
                                                                  0, ]
          finalmat <- finalmat[order(-finalmat[, 6]), ]
          s <- 1
          while (s <= (nrow(finalmat))) {
            rowstart <- finalmat[s, 1]
            rowend <- finalmat[s, 2]
            rowsel <- (finalmat[, 1] <= rowend & finalmat[, 2] >= 
                         rowstart)
            rowsel[s] <- FALSE
            finalmat <- finalmat[!rowsel, ]
            if (is.vector(finalmat)) {
              finalmat <- t(as.matrix(finalmat))
            }
            s <- s + 1
          }
        }
        if (sum(lratio > 0) == 1) {
          finalmat <- (cbind(i, j, yact, lambda, chat, lratio))[lratio > 
                                                                  0, ]
          finalmat <- t(as.matrix(finalmat))
        }
        finalmat <- round(finalmat, digits = 3)
        loglikeij <- cumsum(finalmat[, 6])
        mBIC <- rep(NA, length(loglikeij))
        for (s in 1:nrow(finalmat)) {
          tau <- sort(unique(c(as.vector(finalmat[1:s, 1:2]), 1, num)))
          P <- length(tau) - 2
          mbic <- loglikeij[s]
          mbic <- mbic - 0.5 * sum(log(tau[2:length(tau)] - 
                                         tau[1:(length(tau) - 1)]))
          mbic <- mbic + (0.5 - P) * log(num)
          mBIC[s] <- mbic
        }
        mBIC <- round(mBIC, digits = 3)
        if (mBIC[1] > 0) {
          finalmat <- cbind(rep(sampname_qc[sampno], nrow(finalmat)), 
                            rep(chri, nrow(finalmat)),rep(genei, nrow(finalmat)), finalmat)
          finalmat <- (cbind(finalmat, mBIC)[1:which.max(mBIC), ])
          finalcall <- rbind(finalcall, finalmat)
        }
      }
    }
    finalcall <- finalcall[-1, ]
    if (is.vector(finalcall)) {
      finalcall <- t(as.matrix(finalcall))
    }
    st <- start(refi)[as.numeric(finalcall[, 4])]
    ed <- end(refi)[as.numeric(finalcall[, 5])]
    cnvtype <- rep(NA, length(st))
    cnvtype[as.numeric(finalcall[, 8]) < 2] <- "del"
    cnvtype[as.numeric(finalcall[, 8]) > 2] <- "dup"
    if (nrow(finalcall) == 1) {
      finalcall <- t(as.matrix(c(finalcall[, 1:3], cnvtype, st, ed, (ed - 
                                                                       st + 1)/1000, finalcall[, 4:10])))
    } else {
      finalcall <- cbind(finalcall[, 1:3], cnvtype, st, ed, (ed - st + 
                                                               1)/1000, finalcall[, 4:10])
    }
    colnames(finalcall) <- c("sample_name", "chr", "gene" ,"cnv", "st_bp", "ed_bp", 
                             "length_kb", "st_exon", "ed_exon", "raw_cov", 
                             "norm_cov", "copy_no", "lratio", "mBIC")
    rownames(finalcall) <- rep("", nrow(finalcall))
    finalcall
  }
  
  # source('segment_targeted.R')
  # Available at: https://github.com/yuchaojiang/CODEX2/blob/master/targeted_sequencing/segment_targeted.R
  
  #=====================================================================#
  #                  End of function for segmentation                   #
  #=====================================================================#

optK=which.max(BIC)
finalcall=matrix(ncol=14,nrow=0)
colnames(finalcall)=c('sample_name','chr','gene','cnv',
                      'st_bp','ed_bp','length_kb',
                      'st_exon','ed_exon','raw_cov',
                      'norm_cov','copy_no','lratio',
                      'mBIC')
for(genei in unique(ref_qc$gene)){
  cat('Segmenting gene',genei,'\n')
  geneindex=which(ref_qc$gene==genei)
  yi=Y_qc[geneindex,, drop=FALSE]
  yhati=Yhat[[optK]][geneindex,, drop=FALSE]
  refi=ref_qc[geneindex]
  finalcalli=segment_targeted(yi, yhati, sampname_qc, refi, genei, lmax=length(geneindex), mode='fraction') 
  finalcall=rbind(finalcall,finalcalli)
}
cn=(as.numeric(as.matrix(finalcall[,'copy_no'])))
cn.filter=(cn<=1.7)|(cn>=2.3) # removing calls with fractional copy numbers close to 2 (for heterogeneous cancer samples)
finalcall=finalcall[cn.filter,]
length_exon=as.numeric(finalcall[,'ed_exon'])-as.numeric(finalcall[,'st_exon'])+1
finalcall=cbind(finalcall[,1:7],length_exon,finalcall[,10:14])

finalcall[,"cnv"] <- toupper(finalcall[,"cnv"]) # Transform del/dup to DEL/DUP

write.table(finalcall, file = paste( 'CODEX2', '.results.txt', sep=''), sep='\t', quote=F, row.names=F)

toAnnotateTable <- data.frame(finalcall[,"chr"], finalcall[,"st_bp"], finalcall[,"ed_bp"], toupper(finalcall[,"cnv"]),finalcall[,"sample_name"],finalcall[,"copy_no"]/2)
colnames(toAnnotateTable) = c("CHR", "START", "END", "CNV_TYPE", "SAMPLE","RATIO")

write.table(toAnnotateTable, file = 'CODEX2.toAnnotate.txt', sep='\t', quote=F, row.names=F, col.names = T)


#stop clock
finish <- proc.time() - ptm
print(finish)
summary(finish)
