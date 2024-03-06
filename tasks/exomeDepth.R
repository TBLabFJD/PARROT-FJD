### CNV analysis 
### Author: Lorena de la Fuente
### Date: 14/01/19


rm(list=ls())

#start clock
ptm <- proc.time()


#**************************#
### package requirements ###
#**************************#

library(ExomeDepth)
library(GenomicRanges)
library(optparse)
library(stringr)
# data(exons.hg19)


#**************#
#   Arguments  #
#**************#

option_list=list(
  make_option(c("-d","--dir"),type="character", help="Directory with input bam files to analyse."),
  make_option(c('-o','--outputdir'),type="character", help="Output directory."),
  make_option(c('-n','--name'),type="character", help="Output directory."),
  make_option(c('-b','--bed'),type="character",default="genome", help="Bed file with genes found in panel"),
  make_option(c('-s','--samples'), type="character", help="File with negative samples", default="all"))


opt_parser=OptionParser(option_list = option_list)
opt=parse_args(opt_parser) #list of the args

# cnv_output=paste(opt$outputdir,"/cnvs/ExomeDepth/",sep="")
# dir.create(cnv_output)

sink(paste("software_",opt$name,".txt",sep=""),append=TRUE)
print("R SESSION INFO (ExomeDepth):")
sessionInfo()
sink() 



#**********************#
# Samples to call CNVs #
#**********************#

if(opt$samples!="all"){
  samples = read.delim2(opt$samples, sep=",", stringsAsFactors = F, header=F)[,1]
}




#*************#
#   Regions   #
#*************#

colnamesbed = c("chromosome", "start",  "end", "name")

# if(opt$bed=="genome"){
#   myexons = exons.hg19
#   myexons[,1] = paste("chr",myexons[,1] ,sep = "") 
# }else{
myexons = read.delim(file=opt$bed, header = F, stringsAsFactors = F)
if(length(colnames(myexons))<4){
  cat("ERROR: bed file without NAME field")
  quit(status = 1)
}else if (length(colnames(myexons))>4){
  myexons = myexons[,1:4]}
colnames(myexons) = colnamesbed
#myexons = myexons[!myexons$chromosome %in% c("chrX","chrY", "X","Y"),]
# }


cat('\nTotal number of testing regions: ', nrow(myexons), '\n', sep = "")


#reference.fasta <- '/home/lorena/floridaCluster/broad_bundle_hg19_v/ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.fasta'


#******************************#
#   List of sample bam files   #
#******************************#

setwd(opt$dir)
mybams <- Sys.glob('*.bam')
mybais <- Sys.glob('*.bai')



#***********************************#
#   counts dataframe for all BAMs   #
#***********************************#

cat('\nParsing bam files to get counts...\n')

my.counts <- getBamCounts(bed.frame = myexons, 
                          bam.files = mybams,
                          index.files = mybais,
                          min.mapq = 20,
                          include.chr = FALSE)
                          #referenceFasta = reference.fasta) # fasta if we want to compute GC %

# setwd(cnv_output)

save(my.counts, file = paste('exomedepth',"_counts_image.RData", sep=''))
write.table(my.counts, file = 'exomedepth_count_matrix.txt', sep = "\t", quote = F, row.names = F, col.names = T)
# write.table(my.counts, file = '/home/gonzalo/Documents/cnvs/exomedepth_count_matrix.txt', sep = "\t", quote = F, row.names = F, col.names = T)

ExomeCount.dafr <- my.counts



# ExomeCount.dafr <- as(my.counts[, colnames(my.counts)], 'data.frame')
# 
ExomeCount.dafr$chromosome <- gsub(as.character(ExomeCount.dafr$chromosome),pattern = 'chr',replacement = '')




#******************#
#   count matrix   #
#******************#


ExomeCount.mat <- as.matrix(ExomeCount.dafr[, grep(names(ExomeCount.dafr), pattern = '*.bam')])
#rownames(ExomeCount.mat) = ExomeCount.dafr$exon
nsamples <- ncol(ExomeCount.mat)
#colnames(ExomeCount.mat) = sapply(colnames(ExomeCount.mat), function(x) sub(".bam", "",strsplit(x, split = "_",  fixed = T)[[1]][1]))
colnames(ExomeCount.mat) = sapply(mybams, function(x) strsplit(x, split = ".",  fixed = T)[[1]][1])

if(opt$samples!="all" && !all(samples %in% colnames(ExomeCount.mat))){
  cat('\nSome of the selected samples can not be found...\n')
  
  }
if(opt$samples=="all"){
  samples=colnames(ExomeCount.mat) 
}

cat('\nSamples to analyse: ', paste(samples, collapse = ","), '\n', sep = "")
cat('\nAll detected samples: ',paste(colnames(ExomeCount.mat), collapse=','),'\n', sep = "")



#******************#
#   CNV calling    #
#******************#

CNVcalls_sorted_total.toAnnotate = NULL
CNVcalls_sorted_total = NULL
corDF = NULL

for (i in 1:nsamples) {
  
  samplename = colnames(ExomeCount.mat)[i]
  
  if(samplename %in% samples){
  
    cat("\nAnalysing Sample ",samplename,"...\n",sep = "")
    cat("..................................\n")
    

    #### Create the aggregate reference set for this sample
    
    cat('Read number distribution...\n')
    print(summary(ExomeCount.mat[,i]))
  
    my.choice <- select.reference.set (test.counts = ExomeCount.mat[,i],
                                       reference.counts = ExomeCount.mat[,-i,drop=F],
                                       bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000,
                                       n.bins.reduced = 10000)
  
    # PCA per sample
    # myfactors=data.frame(row.names = colnames(ExomeCount.mat), sample=colnames(ExomeCount.mat)=="17-2983")
    # myfactors[which(rownames(myfactors)%in%my.choice$reference.choice),"sample"]<-"reference"
    # myfactors[which(rownames(myfactors)==samplename),"sample"]<-"test"
    # mydata2 = readData(ExomeCount.mat,factors = myfactors)
    # myPCA = dat(mydata2, type = "PCA")
    # explo.plot(myPCA, factor = "sample")
    
    cat('Reference set choice: ', paste(my.choice$reference.choice, collapse = ","),"\n",sep = "")
    
    my.reference.selected <- apply(X = ExomeCount.mat[, my.choice$reference.choice, drop = FALSE],
                                      MAR = 1,
                                      FUN = sum)
    
    
    cat('\nComputing correlation between sample and references...\n')
    
    my.test = ExomeCount.mat[, i, drop = FALSE]
    my.ref.counts = ExomeCount.mat[, my.choice$reference.choice, drop = FALSE]
    correlation = cor(my.test[,1], apply(my.ref.counts,1,mean))
    
    corDF = rbind(corDF, c(samplename, round(correlation,4)))
    
    if (correlation > 0.97){
      
      cat(paste("Correlation between reference and tests count is",  round(correlation,4), ". Running CNV calling"))
      
      cat('Creating the ExomeDepth object...\n')
      
      all.exons <- new('ExomeDepth',
                       test = ExomeCount.mat[,i],
                       reference = my.reference.selected,
                       formula = 'cbind(test, reference) ~ 1')
  
      cat('Calling CNVs...\n')
      
      
      ################ Now call the CNVs
      
      all.exons <- CallCNVs(x = all.exons,
                            transition.probability = 10^-4,
                            chromosome = ExomeCount.dafr$chromosome,
                            start = ExomeCount.dafr$start,
                            end = ExomeCount.dafr$end,
                           # name = ExomeCount.dafr$names)
                            name = ExomeCount.dafr$exon)
      
      cat("Detected CNVs:",nrow(all.exons@CNV.calls), "\n", sep=" ")
      
      if (nrow(all.exons@CNV.calls) > 0) {
      
        # Adding name found in panel
        
        rangesMyAnnotation = GenomicRanges::GRanges(seqnames = all.exons@annotations$chromosome,
                                                    IRanges::IRanges(start=all.exons@annotations$start,end=all.exons@annotations$end),
                                                    names = all.exons@annotations$name)
        all.exons <- AnnotateExtra(x = all.exons,
                                   reference.annotation = rangesMyAnnotation,
                                   min.overlap = 0.0000001,
                                   column.name =
                                     'annotation')
        
        # Writting file for VEP annotation
        
        CNVcalls_sorted = all.exons@CNV.calls[ order ( all.exons@CNV.calls$BF, decreasing = TRUE),]

        CNVcalls_sorted[CNVcalls_sorted$type=="duplication","CNV_type"]<-"DUP"
        CNVcalls_sorted[CNVcalls_sorted$type=="deletion","CNV_type"]<-"DEL"

        CNVcalls_sorted_total = rbind( CNVcalls_sorted_total, data.frame(samplename,CNVcalls_sorted))
        # CNVcalls_sorted_total.toAnnotate = unique(rbind(CNVcalls_sorted_total.toAnnotate, CNVcalls_sorted[,c("chromosome","start","end","CNV_type"),]))
    
    
    }
    
    }else{cat(paste("WARNING: Correlation between reference and tests count is",  round(correlation,4), ". Skipped CNV calling.\n"))}
  }

}

toAnnotateTable = CNVcalls_sorted_total[,c("chromosome","start","end","CNV_type","samplename", "BF", "reads.ratio"), drop=F]
colnames(toAnnotateTable) = c("CHR", "START", "END", "CNV_TYPE", "SAMPLE","BF", "RATIO")
toAnnotateTable$SAMPLE = gsub("_.*", "", toAnnotateTable$SAMPLE)
toAnnotateTable$CHR[grep("chr", toAnnotateTable$CHR,invert = T)] = paste0("chr", toAnnotateTable$CHR[grep("chr", toAnnotateTable$CHR,invert = T)])

colnames(corDF) = c("sample", "cor")

output.file = paste('exomedepth','.results.txt', sep = '')
forvep.file = paste('exomedepth','.toAnnotate.txt', sep = '')
cor.file = paste('exomedepth',  ".correlation", '.txt', sep = '')

cat("\nOUTPUT FILES:\n")
cat(paste(output.file,"\n"))
cat(paste(forvep.file,"\n"))
cat(paste(cor.file,"\n"))

write.table(x = CNVcalls_sorted_total, file = output.file, row.names = FALSE, sep="\t", quote=F)
write.table(x = toAnnotateTable, file=forvep.file, quote=F, col.names = T, row.names=F, sep="\t")
write.table(x = corDF, file = cor.file, row.names = FALSE, sep="\t", quote=F)

cat("\nEXOMEDEPTH CNV DETECTION FINISHED\n")

#stop clock
finish <- proc.time() - ptm
print(finish)
summary(finish)

