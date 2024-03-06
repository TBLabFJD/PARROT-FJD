



library(optparse)



#############
# Arguments # 
#############
option_list = list(
    make_option(c("-s", "--samplename"), type="character", default=NULL, help="sample", metavar="character"),
    make_option(c("-g", "--genomecov_path"), type="character", default=NULL,help="Genomic coverage", metavar="character"),
    make_option(c("-q", "--quality_path"), type="character", default=NULL, help="Basecall quality summary", metavar="character"),
    make_option(c("-f", "--samtools_flagstat_path"), type="character", default=NULL, help="Samtools flagstat output", metavar="character"),
    make_option(c("-l", "--read_lenghts_path"), type="character", default=NULL, help="Read lengths", metavar="character"),
    make_option(c("-c", "--cg_at_path"), type="character", default=NULL, help="GG and AT content", metavar="character"),
    make_option(c("-d", "--bedout"), type="character", default=NULL, help="Output bed file", metavar="character"),
    make_option(c("-p", "--panelcov_path"), type="character", default=NULL, help="Panel coverage", metavar="character"),
    make_option(c("-b", "--bed_path"), type="character", default=NULL, help="bed file", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, help="output file", metavar="character"),
    make_option(c("-u", "--nreads_nondup_uniq"), type="integer", default=NULL, help="Number of uniquely mapped, non-duplicated reads", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


sample = opt$samplename
genomecov_path = opt$genomecov_path
quality_path = opt$quality_path
samtools_flagstat_path = opt$samtools_flagstat_path
read_lenghts_path = opt$read_lenghts_path
cg_at_path = opt$cg_at_path
# n50_n90_path = opt$n50_n90_path
panelcov_path = opt$panelcov_path
bed_path = opt$bed_path
output = opt$output
nreads_nondup_uniq = opt$nreads_nondup_uniq
bedout = opt$bedout

# # 
# sample="22-0087-CES178"
# genomecov_path="/home/gonzalo/tblab//mnt/tblab/gonzalo/reanalysis/aAvila111122/work/c3/14488ff6c9220181a85dd7437a0356/22-1698.genomecov.bed"
# quality_path="/home/gonzalo/tblab/mnt/tblab/gonzalo/reanalysis/aAvila111122/work/88/858f07b9091e8a112b4798d40cd6ad/22-1698.quality.txt"
# samtools_flagstat_path="/home/gonzalo/tblab//mnt/tblab/gonzalo/reanalysis/aAvila111122/work/c9/15ae3d5149e86e4c1ccfb6d584f02c/22-1698.samtools_flagstat.txt"
# read_lenghts_path="/home/gonzalo/tblab//mnt/tblab/gonzalo/reanalysis/aAvila111122/work/c6/e3511ec10c91c341cd7d7e2c0a27d9/22-1698.read_lenghts.txt"
# cg_at_path="/home/gonzalo/tblab//mnt/tblab/gonzalo/reanalysis/aAvila111122/work/27/b09da2712fc6ac1cd590c5c4a0de75/22-1698.CG_AT.txt"
# # n50_n90_path="/home/gonzalo/tblab/"
# panelcov_path="/home/gonzalo/tblab//mnt/tblab/gonzalo/reanalysis/aAvila111122/work/c3/14488ff6c9220181a85dd7437a0356/22-1698.panelcov.bed"
# bed_path="/home/gonzalo/tblab//mnt/genetica7/beds/CES_v3_hg38_target.chr.formatted.sorted.annotated.bed"
# output="/home/gonzalo/tblab/Downloads/22-1698.quality.summary.txt"
# nreads_nondup_uniq="/home/gonzalo/tblab//mnt/tblab/gonzalo/reanalysis/aAvila111122/work/22/db6cbcc2a26d6fa102dcb1f57e6a84/22-1698.nreads_nondup_uniq.txt"
# bedout="/home/gonzalo/Downloads/22-1698.library.stats.txt"


################
# Data loading #
################
genomecov=read.delim(genomecov_path, header = F, stringsAsFactors = F)

quality=read.delim(quality_path, header = F, quote = "", sep = " ", stringsAsFactors = F)

samtools_flagstat=read.delim(samtools_flagstat_path, header = F, stringsAsFactors = F)
samtools_flagstat = data.frame(samtools_flagstat)
samtools_flagstat$V2 = as.numeric(gsub(" \\+.*", "", samtools_flagstat$V1, perl=T))
samtools_flagstat$V3 = gsub(".* \\+ [0-9]+ ", "", samtools_flagstat$V1, perl=T)
samtools_flagstat$V3 = gsub(" \\(mapQ", "_(mapQ", samtools_flagstat$V3, perl=T)
samtools_flagstat$V3 = gsub(" \\(.*", "", samtools_flagstat$V3, perl=T)
samtools_flagstat$V3 = gsub(" ", "_", samtools_flagstat$V3, perl=T)
row.names(samtools_flagstat) = samtools_flagstat$V3
samtools_flagstat = data.frame(t(samtools_flagstat[,2,drop=F]), stringsAsFactors = F, check.names = F)


read_lenghts=read.delim(read_lenghts_path, header = F, stringsAsFactors = F)

cg_at=read.delim(cg_at_path, header = F, quote = "", sep = " ", stringsAsFactors = F)

# n50_n90=read.delim(n50_n90_path, header = T, quote = "", stringsAsFactors = F)




df_out = data.frame(row.names = 1, stringsAsFactors = F)

save.image(file = "Rdata.RData")


#####################
# Basic Information #
#####################
df_out$Sample = sample
df_out$Raw_total_read_bases = sum(quality$V1)
df_out$Raw_total_read_number = samtools_flagstat$primary



# Quality
# chr <- function(n) { rawToChar(as.raw(n)) }
# aa = c()
# for (i in 33:73) aa = c(aa, chr(i))
qscore <- function(x) { strtoi(charToRaw(as.character(x)),16L)-33 }
quality$V3=0
for(i in 1:nrow(quality)) quality$V3[i] = qscore(quality$V2[i])
df_out$Q20=round(sum(quality$V1[quality$V3 >= 20])/sum(quality$V1)*100, 2) # Percentage of bases with a quality equal or higher than 20
df_out$Q30=round(sum(quality$V1[quality$V3 >= 30])/sum(quality$V1)*100, 2) # Percentage of bases with a quality equal or higher than 30



# cg_at
df_out$`GC(%)` = round( sum(cg_at$V1[cg_at$V2 %in% c("C", "c", "G", "g")]) / sum(cg_at$V1) * 100, 2)
df_out$`AT(%)` = round( sum(cg_at$V1[cg_at$V2 %in% c("A", "a", "T", "t")]) / sum(cg_at$V1) * 100, 2)
df_out$`N(%)` = 100 - (df_out$`GC(%)` + df_out$`AT(%)`)



# Read length info
length_info = data.frame(t(as.matrix(summary(read_lenghts$V1))))
colnames(length_info) = paste("Read_length", colnames(length_info), sep = "_")
df_out = cbind(df_out, length_info)
# df_out = cbind(df_out, n50_n90)

length_distribution=table(read_lenghts$V1)
length_distribution_mod=length_distribution*as.numeric(names(length_distribution))
n50_n90_calculator = function(selected_pos){
    readlengths=as.numeric(names(length_distribution_mod))
    for (i in 1:length(length_distribution_mod)){
        selected_pos = selected_pos - length_distribution_mod[i]
        if (selected_pos < 0) {return(readlengths[i])}
        else if (selected_pos == 0) {return((readlengths[i] + readlengths[i+1]) / 2)}
    }
}
df_out$N50 = n50_n90_calculator(sum(length_distribution_mod)/2)
df_out$N90 = n50_n90_calculator(sum(length_distribution_mod)/10)



# Mapping stats
df_out$mapped_reads = samtools_flagstat$mapped
df_out$`mapping_ratio(%)` = samtools_flagstat$mapped / samtools_flagstat$in_total * 100
df_out$primary_mapped_reads = samtools_flagstat$primary_mapped
df_out$`primary_mapping_ratio(%)` = samtools_flagstat$primary_mapped / samtools_flagstat$primary * 100
df_out$uniquemapped_nondup_reads = nreads_nondup_uniq
df_out$`uniquemapped_nondup_ratio(%)` = as.numeric(nreads_nondup_uniq) / samtools_flagstat$primary * 100

df_out$secondary_reads = samtools_flagstat$secondary
df_out$supplementary_reads = samtools_flagstat$supplementary
df_out$duplicated_reads = samtools_flagstat$duplicates

df_out$paired_in_sequencing = samtools_flagstat$paired_in_sequencing
df_out$read1 = samtools_flagstat$read1
df_out$read2 = samtools_flagstat$read2
df_out$properly_paired = samtools_flagstat$properly_paired
df_out$itself_and_mate_mapped = samtools_flagstat$with_itself_and_mate_mapped
df_out$singletons = samtools_flagstat$singletons
df_out$mate_mapped_to_a_different_chr = samtools_flagstat$with_mate_mapped_to_a_different_chr
df_out$`mate_mapped_to_a_different_chr(mapQ>=5)` = samtools_flagstat$`with_mate_mapped_to_a_different_chr_(mapQ>=5)`




# Coverage
genomecov$V5 = genomecov$V3 - genomecov$V2
df_out$`total_covered_bases_with(DP>1read)`= sum(genomecov[genomecov$V4 >= 1,"V5"])
df_out$`total_covered_bases_with(DP>10reads)`= sum(genomecov[genomecov$V4 >= 10,"V5"])
df_out$`total_covered_bases_with(DP>20reads)`= sum(genomecov[genomecov$V4 >= 20,"V5"])
df_out$`total_covered_bases_with(DP>30reads)`= sum(genomecov[genomecov$V4 >= 30,"V5"])
df_out$`total_covered_bases_with(DP>40reads)`= sum(genomecov[genomecov$V4 >= 40,"V5"])
df_out$`total_covered_bases_with(DP>50reads)`= sum(genomecov[genomecov$V4 >= 50,"V5"])
df_out$`total_covered_bases_with(DP>75reads)`= sum(genomecov[genomecov$V4 >= 75,"V5"])
df_out$`total_covered_bases_with(DP>100reads)`= sum(genomecov[genomecov$V4 >= 100,"V5"])




# Coverage (panel)
if (!is.null(panelcov_path)){
    panelcov=read.delim(panelcov_path, header = F, stringsAsFactors = F)
    bed=read.delim(bed_path, header = F, stringsAsFactors = F)
    
    panelcov$fragment_size = panelcov$V3 - panelcov$V2
    panel_bases = sum(panelcov$fragment_size)
    
    cov_vector = rep.int(panelcov$V4, panelcov$fragment_size)

    cov_info = data.frame(t(as.matrix(summary(cov_vector))))
    colnames(cov_info) = paste("Panel_cov_dp", colnames(cov_info), sep = "_")
    df_out = cbind(df_out, cov_info)

    df_out$`panel_covered_ratio_with(DP>1read)`= sum(panelcov[panelcov$V4 >= 1,"fragment_size"]) / panel_bases * 100
    df_out$`panel_covered_bases_with(DP>10reads)`= sum(panelcov[panelcov$V4 >= 10,"fragment_size"]) / panel_bases * 100
    df_out$`panel_covered_bases_with(DP>20reads)`= sum(panelcov[panelcov$V4 >= 20,"fragment_size"]) / panel_bases * 100
    df_out$`panel_covered_bases_with(DP>30reads)`= sum(panelcov[panelcov$V4 >= 30,"fragment_size"]) / panel_bases * 100
    df_out$`panel_covered_bases_with(DP>40reads)`= sum(panelcov[panelcov$V4 >= 40,"fragment_size"]) / panel_bases * 100
    df_out$`panel_covered_bases_with(DP>50reads)`= sum(panelcov[panelcov$V4 >= 50,"fragment_size"]) / panel_bases * 100
    df_out$`panel_covered_bases_with(DP>75reads)`= sum(panelcov[panelcov$V4 >= 75,"fragment_size"]) / panel_bases * 100
    df_out$`panel_covered_bases_with(DP>100reads)`= sum(panelcov[panelcov$V4 >= 100,"fragment_size"]) / panel_bases * 100
    
    # library(ggplot2)
    # coverage_df_plot = data.frame(table(cov_vector))
    # coverage_df_plot$cov_vector = as.numeric(as.character(coverage_df_plot$cov_vector))
    # 
    # ggplot(data = coverage_df_plot, aes(x=cov_vector, y=Freq)) +
    #     geom_line() +
    #     theme_bw()
    
    
    panelcov_split = split(panelcov, f=paste(panelcov$V5, panelcov$V6, panelcov$V7, panelcov$V8, sep = "___"))
    bed_stats = do.call("rbind", lapply(panelcov_split, function(x){
        dpvalues = rep.int(x$V4, x$V3-x$V2)
        p10x = round(table(dpvalues>=10)["TRUE"]/length(dpvalues)*100, 2)
        return(c(summary(dpvalues), p10x = p10x))
    }))
    bed_stats=cbind(do.call("rbind", strsplit(rownames(bed_stats),split = "___")), bed_stats)
    colnames(bed_stats)=c("#chr", "start", "end", "gene", colnames(bed_stats)[5:10], "ratio_DP>10x(%)")
    bed_stats = data.frame(bed_stats, stringsAsFactors = F, check.names=F)
    bed_stats$start = as.numeric(bed_stats$start)
    bed_stats$end = as.numeric(bed_stats$end)
    bed_stats$Mean = round(as.numeric(bed_stats$Mean), 2)
    bed_stats = bed_stats[order(bed_stats$`#chr`, bed_stats$start, bed_stats$end),]
    bed_stats$size = bed_stats$end - bed_stats$start
    write.table(bed_stats, bedout, row.names = F, col.names = T, sep = "\t", quote = F)
}else{
    print("No panel used")
}


write.table(df_out, output, row.names = F, col.names = T, sep = "\t", quote = F)






