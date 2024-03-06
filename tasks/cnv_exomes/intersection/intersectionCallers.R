#### INTERSECTION CALLERS AND BENCHMARK

library(readxl)
library(reshape2)
library(plyr)
library(vcfR)


####### LOAD AND PROCESS RESULTS

ED_raw = read.table("uamSSHFS/SVbenchmark/CODEX2/individual_copynumber/copy_number_variation_data/ExomeDepth/resultswoGC/exomedepth_results.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
ED = ED_raw[,c("samplename", "chromosome", "start", "end", "type")]
colnames(ED) = c("sample","chr", "start", "end", "cnv_type")
ED$chr = paste("chr",ED$chr, sep = "")
ED$cnv_type = factor(ED$cnv_type, levels = c("deletion", "duplication") , labels = c("DEL", "DUP"))

#ED_f = ED[ED$type%in%c("deletion") & ED$sample=="NA06985",2:4]
#write.table(ED, "ed_example.bed", quote=F, row.names = F, col.names = F, sep = "\t")

load("uamSSHFS/SVbenchmark/CODEX2/individual_copynumber/copy_number_variation_data/CODEX2/CODEX2_CNV_nr.RData")
finalcall.CBS.df <- do.call("rbind", finalcall.CBS)
C2nr = finalcall.CBS.df[,c("sample_name", "chr", "st_bp", "ed_bp", "cnv")]
colnames(C2nr) = c("sample","chr", "start", "end", "cnv_type")
C2nr$chr = as.character(C2nr$chr)
C2nr$cnv_type = factor(C2nr$cnv_type, levels = c("del", "dup") , labels = c("DEL", "DUP"))

#C2nr_f = C2nr[C2nr$type=="del" & C2nr$sample=="NA06985",2:4]
#write.table(C2nr, "c2_example.bed", quote=F, row.names = F, col.names = F, sep = "\t")

# CS_raw = read.table("uamSSHFS/SVbenchmark/CODEX2/individual_copynumber/copy_number_variation_data/CNMOPS/CNMOPS_results.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
# CS_raw$CN = gsub(pattern = "CN", replacement = "", CS_raw$CN)
# CS = CS_raw[,c("sampleName", "seqnames", "start", "end", "CN")]
# CS$CN = as.character(CS$CN>2)
# colnames(CS) = c("sample","chr", "start", "end", "cnv_type")
# CS$chr = as.character(CS$chr)
# CS$cnv_type = factor(CS$cnv_type, levels = c("FALSE", "TRUE") , labels = c("DEL", "DUP"))

#CS_f = CS[CS$type%in%c("CN0","CN1") & CS$sample=="NA06985",2:4]
#write.table(CS_f, "cs_example.bed", quote=F, row.names = F, col.names = F, sep = "\t")

XH_raw = read.table("uamSSHFS/SVbenchmark/CODEX2/individual_copynumber/copy_number_variation_data/XHMM/DATA.xcnv", header = T, sep = "\t", stringsAsFactors = FALSE)
XH_raw$start = as.integer(sapply(strsplit(XH_raw$INTERVAL, fixed =T, split = "-"), function(x) strsplit(x[1], split = ":", fixed = T)[[1]][2]))
XH_raw$end = as.integer(sapply(strsplit(XH_raw$INTERVAL, fixed =T, split = "-"), function(x) x[2]))
XH = XH_raw[,c("SAMPLE", "CHR","start","end", "CNV")]
colnames(XH) = c("sample","chr", "start", "end", "cnv_type")
XH$chr = as.character(XH$chr)
XH$cnv_type = factor(XH$cnv_type, levels = c("DEL", "DUP") , labels = c("DEL", "DUP"))



gatk_files = list.files("uamSSHFS/SVbenchmark/CODEX2/individual_copynumber/copy_number_variation_data/gatk/", pattern = "vcf.gz$", full.names = TRUE)
gatk_files = grep("segments", gatk_files, value = TRUE)

gatk_raw = NULL
for (file in gatk_files){
  
  a = read.vcfR(file)
  sample = colnames(a@gt)[2]
  aa <- data.frame(sample=sample, a@fix[which(a@fix[,"ALT"]!="NA"),c(1,2,8,5)], stringsAsFactors = F)
  colnames(aa) = c("sample","chr", "start", "end", "cnv_type")
  aa$end = as.integer(gsub("END=", "", aa$end))
  aa$start = as.integer(aa$start)
  #aa$QUAL = as.integer(aa$QUAL)
  aa$cnv_type = factor(aa$cnv_type, levels = c("<DEL>", "<DUP>") , labels = c("DEL", "DUP"))
  gatk_raw[[sample]] = aa
}

gatk = do.call("rbind", gatk_raw)



myProgramList = list(ED=ED, C2=C2nr, XH=XH, gatk=gatk)




######## RUN INTERSECION


source("uamSSHFS/SVbenchmark/CODEX2/multiProgramIntersection.R")

df_intersection = multiIntersectCNVs(list = myProgramList)
df_intersection$overlaplength = df_intersection$V3-df_intersection$V2
df_intersection_min1 = df_intersection[df_intersection$overlaplength>1,]

df_intersection_ok = df_intersection[df_intersection$n.overlaps>2,c(7,1,2,3,6)]
colnames(df_intersection_ok) = c("sample","chr", "start", "end", "cnv_type")

df_intersection_ok_min1 = df_intersection_min1[df_intersection_min1$n.overlaps>2,c(7,1,2,3,6)]
colnames(df_intersection_ok_min1) = c("sample","chr", "start", "end", "cnv_type")

write.table(df_intersection_ok_min1[,c("chr", "start", "end", "sample","cnv_type")], file="uamSSHFS/SVbenchmark/OpenEBench/results/structuralvariants-master/Results_pp/CNV_intersection_outputFJD_1bp.tsv",
            quote=F, row.names = F, sep = "\t")

# df_intersection_ok2 = df_intersection[df_intersection$n.overlaps>1,c(7,1,2,3,6)]
# colnames(df_intersection_ok2) = c("sample","chr", "start", "end", "cnv_type")


sample = "NA06985"
gs.bed.sorted <- bedr.sort.region(df_intersection_ok[df_intersection_ok$sample==sample,c(2,3,4)], check.chr =  TRUE, check.merge = TRUE)
is.merged <- is.merged.region(gs.bed.sorted);
print(is.merged)

gs.bed.sorted <- bedr.sort.region(df_intersection_ok_min1[df_intersection_ok_min1$sample==sample,c(2,3,4)], check.chr =  TRUE, check.merge = TRUE)
is.merged <- is.merged.region(gs.bed.sorted);
print(is.merged)




## Loading TEST RESULTS

prog.list = excel_sheets("uamSSHFS/SVbenchmark/CODEX2/tosend/CNVs_results.xlsx")
test = sapply( prog.list, function(x) data.frame( read_xlsx("uamSSHFS/SVbenchmark/CODEX2/tosend/CNVs_results.xlsx", sheet = x), stringsAsFactors = F))

# results as they would be provided by TransBioNET groups (sample  chr     start       end cnv_type).

test.formated = NULL
#test.formated$CODEX2 = test$CODEX2[, c(1,2,4,5,3)]  
#test.formated$CODEX = data.frame( sample=test$CODEX$sample_name, chr=paste("chr",as.character(test$CODEX$chr), sep = "" ), start=test$CODEX$st_bp, end=test$CODEX$ed_bp, cnv_type=test$CODEX$copy_no>2, stringsAsFactors = FALSE)
#test.formated$XHMM = data.frame( sample=test$XHMM$Sample, chr=paste("chr",test$XHMM$CHR, sep = "" ),  start=as.numeric(sapply(test$XHMM$INTERVAL, function(x) strsplit(x,split="[:,-]")[[1]][2])), end=as.numeric(sapply(test$XHMM$INTERVAL, function(x) strsplit(x,split="[:,-]")[[1]][3])),cnv_type=test$XHMM$CNV, stringsAsFactors = FALSE)
#test.formated$EXCAVATOR = data.frame( sample=test$EXCAVATOR$Sample, chr=paste("chr",test$EXCAVATOR$Chromosome, sep = "" ), start=test$EXCAVATOR$Start, end=test$EXCAVATOR$End, cnv_type=test$EXCAVATOR$CN>2, stringsAsFactors = FALSE)
#test.formated$CLAMMS = data.frame( sample=test$CLAMMS$sample, chr=paste("chr",test$CLAMMS$chr, sep = "" ), start=test$CLAMMS$start, end=test$CLAMMS$end, cnv_type=test$CLAMMS$cnv_type, stringsAsFactors = FALSE)
#test.formated$CLAMMS = test.formated$CLAMMS[test$CLAMMS$num_exons>1 & test$CLAMMS$num_exons<21,]
#test.formated$EXCAVATOR$cnv_type = factor(test.formated$EXCAVATOR$cnv_type, levels = c(FALSE, TRUE) , labels = c("DEL", "DUP"))
#test.formated$CODEX$cnv_type = factor(test.formated$CODEX$cnv_type, levels = c(FALSE, TRUE) , labels = c("DEL", "DUP"))


# Files

challenge.results = list.files(path = "Downloads/structuralvariants-master/Results/", "CNV")

JORGE_raw = read.table("Downloads/structuralvariants-master/Results/CNVs_FPGMX_IDIS.tab",
                       header = TRUE, sep = "\t",  stringsAsFactors = FALSE)
JORGE = JORGE_raw[,c("sample","chr", "start", "end", "cnv_type")]
colnames(JORGE) <- colnames(test.formated$CODEX2)
JORGE$chr = paste("chr",JORGE$chr, sep = "")
test.formated$JORGE <- JORGE


#JORGE AMIGO
JORGE_raw = read.table("Downloads/CNVs_FPGMX_IDIS.tab",
                       header = TRUE, sep = "\t",  stringsAsFactors = FALSE)
JORGE = JORGE_raw[JORGE_raw$bayes_factor>6,c("sample","chr", "start", "end", "cnv_type")]
colnames(JORGE) <- colnames(test.formated$CODEX2)
JORGE$chr = paste("chr",JORGE$chr, sep = "")
test.formated$JORGEfilt <- JORGE

#ANGELA
ANGELA_raw = read.table("Downloads/cnvs_test_bwa_3tools_restrictive.tsv",
                       header = TRUE, sep = "\t",  stringsAsFactors = FALSE)
ANGELA = ANGELA_raw[,c("sample","chrom", "start", "end")]
ANGELA$cnv_type = "DEL"
colnames(ANGELA) <- colnames(test.formated$CODEX2)
ANGELA$chr = paste("chr",ANGELA$chr, sep = "")
test.formated$ANGELA <- ANGELA
test.formated$ANGELA$sample <- sapply(test.formated$ANGELA$sample, function(x) strsplit(x, split = ".", fixed = T)[[1]][1])



#OTHER
OTHER_raw = read.table("Downloads/Results_CNVfilter_outputSJD.tsv",
                       header = TRUE, sep = "\t",  stringsAsFactors = FALSE)
OTHER = OTHER_raw[,c("sample","chr", "start", "end", "cnv_type")]
colnames(OTHER) <- colnames(test.formated$CODEX2)
OTHER$chr = paste("chr",OTHER$chr, sep = "")
test.formated$OTHER <- OTHER



test.formated$myED = ED
test.formated$myC2 = C2nr
test.formated$myXH = XH
test.formated$myINT2of4 = df_intersection_ok2
test.formated$myGATK = gatk
test.formated$myINT3of4 = df_intersection_ok
test.formated$ANGELA = ANGELA


test.formated2=NULL

test.formated2$myGATK = gatk
test.formated2$myINT4 = df_intersection_ok



######## RUN BENCHMARK


source("uamSSHFS/SVbenchmark/CODEX2/Benchmark_functions.R")


bench.list = benchmarkCNVs(test.formated = test.formated)
benchSTRINGENT.list = benchmarkCNVs(test.formated = test.formated)

bench.list2 = benchmarkCNVs(test.formated = test.formated2)
bench.list=bench.list2
#save(bench.list, bench.list2, file="uamSSHFS/SVbenchmark/CODEX2/individual_copynumber/results.RData")



#******** Plotting

  res.df = data.frame(rbind(bench.list$res.df, bench.list$res.df), stringsAsFactors = FALSE)
  res.df$Recall = as.numeric(res.df$Recall)
  res.df$Precision = as.numeric(res.df$Precision)
  a<-ggplot(data = res.df, aes(x=Precision, y=Recall, color=V3)) + geom_point(size=4) +
    facet_grid(cols = vars(V2)) +  theme_bw() + ylim(c(0,1)) +  xlim(c(0,1))
  
  res.rare.df = data.frame(rbind(bench.list$res.rare.df, bench.list$res.rare.df), stringsAsFactors = FALSE)
  res.rare.df$Recall = as.numeric(res.rare.df$Recall)
  res.rare.df$Precision = as.numeric(res.rare.df$Precision)
  b<-ggplot(data = res.rare.df, aes(x=Precision, y=Recall, color=V3)) + geom_point(size=4) +
    facet_grid(cols = vars(V2)) +  theme_bw() + ylim(c(0,1)) +  xlim(c(0,1))
  
  res.common.df = data.frame(rbind(bench.list$res.common.df, bench.list$res.common.df), stringsAsFactors = FALSE)
  res.common.df$Recall = as.numeric(res.common.df$Recall)
  res.common.df$Precision = as.numeric(res.common.df$Precision)
  c<-ggplot(data = res.common.df, aes(x=Precision, y=Recall, color=V3)) + geom_point(size=4) +
    facet_grid(cols = vars(V2)) +  theme_bw() + ylim(c(0,1)) +  xlim(c(0,1))
  
  aa = rbind(res.common.df, res.rare.df, res.df)
  
  ggplot(data = aa, aes(x=Precision, y=Recall, color=V3, shape=V1)) + geom_point(size=2) +
    facet_grid(cols = vars(V2)) +  theme_bw() + ylim(c(0,1)) +  xlim(c(0,1)) +
    scale_shape_manual(values = c(15,16,17))







    # 
# test.formated$CS_tuned <- CS_tuned
# 
# ED_raw = read.table("uamSSHFS2/SVbenchmark/CODEX2/individual_copynumber/copy_number_variation_data/ExomeDepth/resultswoGC/exomedepth_results.txt", header = T, sep = "\t", stringsAsFactors = FALSE)
# ED = ED_raw[,c("samplename", "chromosome", "start", "end", "type")]
# colnames(ED) <- colnames(test.formated$CODEX2)
# ED$chr = paste("chr",ED$chr, sep = "")
# test.formated$ED <- ED
# write.table(ED[,c(2:4)], file="uamSSHFS/SVbenchmark/CODEX2/individual_copynumber/copy_number_variation_data/ExomeDepth/exomedepth_results.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
# 
# load("uamSSHFS2/SVbenchmark/CODEX2/individual_copynumber/copy_number_variation_data/CODEX2/CODEX2_CNV_nr.RData")
# finalcall.CBS.df <- do.call("rbind", finalcall.CBS)
# C2nr = finalcall.CBS.df[,c("sample_name", "chr", "st_bp", "ed_bp", "cnv")]
# C2nr$chr = as.character(C2nr$chr)
# colnames(C2nr) <- colnames(test.formated$CODEX2)
# test.formated$C2nr <- C2nr
# write.table(C2nr[,c(2:4)], file="uamSSHFS/SVbenchmark/CODEX2/individual_copynumber/copy_number_variation_data/CODEX2/codex2_results.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
# 


# df_int=NULL
# 
# df_int = lapply(samples, function(x){
# 
#   
#   # duplications
#   
#   C2bed <- bedr.sort.region(C2nr[C2nr$sample==x & C2nr$cnv_type=="dup",c(2:4)], check.chr = FALSE)
#   EDbed <- bedr.sort.region(ED[ED$sample==x & ED$cnv_type=="duplication",c(2:4)], check.chr = FALSE)
#   
#   colnames(C2bed) <- c("chr", "start","end")
#   colnames(EDbed) <- c("chr", "start","end")
#   
# 
#   overlapping.exons.dup <- bedr.join.region(
#     C2bed,
#     EDbed,
#     report.n.overlap = TRUE,
#     check.chr = FALSE,
#     check.merge = FALSE,
#     check.sort = TRUE
#   )
#   
#   dup1 = overlapping.exons.dup[,1:3]
#   dup2 = overlapping.exons.dup[,4:6]
#   colnames(dup2) = c("chr","start","end")
#   dup2$start = as.integer(dup2$start)
#   dup2$end = as.integer(dup2$end)
#   dup = rbind(dup1, dup2)
#   dup <- bedr.sort.region(dup, check.chr = FALSE)
#   
#   merge.dup <- bedr.merge.region(dup)
#   
#   
#   # deletions
#   
#   
#   C2bed = NULL
#   EDbed = NULL
#   C2bed <- bedr.sort.region(C2nr[C2nr$sample==x & C2nr$cnv_type=="del",c(2:4)], check.chr = FALSE)
#   EDbed <- bedr.sort.region(ED[ED$sample==x & ED$cnv_type=="deletion",c(2:4)], check.chr = FALSE)
#   
#   colnames(C2bed) <- c("chr", "start","end")
#   colnames(EDbed) <- c("chr", "start","end")
#   
#   overlapping.exons.del <- bedr.join.region(
#     C2bed,
#     EDbed,
#     report.n.overlap = TRUE,
#     check.chr = FALSE,
#     check.merge = FALSE,
#     check.sort = TRUE
#   )
#   
#   del1 = overlapping.exons.del[,1:3]
#   del2 = overlapping.exons.del[,4:6]
#   colnames(del2) = c("chr","start","end")
#   del2$start = as.integer(del2$start)
#   del2$end = as.integer(del2$end)
#   del = rbind(del1, del2)
#   del <- bedr.sort.region(del, check.chr = FALSE)
# 
#   merge.del <- bedr.merge.region(del)
#   
#   
#   
#   
#   # joined 
#   
#   df_int = rbind(df_int,rbind(data.frame(merge.del, cnv_type="DEL", sample=x), data.frame(merge.dup, cnv_type="DUP", sample=x)))
#   
# 
#   return(df_int)
# 
#   }
# 
# )
#   
# 
# df_int_rbind = do.call("rbind", df_int)
#   
# 
# 
# 
# 
# 
# 
# 
