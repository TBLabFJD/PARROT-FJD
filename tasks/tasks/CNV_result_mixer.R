### CNV analysis 
### Author: Gonzalo Núñez Moreno
### Date: 22/03/19


rm(list=ls())


#**********************#
# package requirements #
#**********************#

library(optparse)


#**************#
#   Arguments  #
#**************#
option_list=list(
  make_option(c('-o','--outputdir'),type="character", help="Output directory."),
  make_option(c('-n','--name'),type="character", help="Project name."),
  make_option(c('-b','--bed'),type="character",default="genome", help="Bed file with genes found in panel"),
  make_option(c('-s','--samples'), type="character", help="Samples to keep"),
  make_option(c("-d","--dir"),type="character", help="Directory with input bam files to analyse."))


opt_parser=OptionParser(option_list = option_list)
opt=parse_args(opt_parser) #list of the args


# opt=list()
# opt$name <- "M96_2019_03_29_10_09_14"
# opt$outputdir <- "/mnt/genetica/gonzalo/Prueba_pipeline/output/"
# opt$bed <- "/mnt/genetica/CIBERERBioinfo/CNVs/ICR96_dataset/bed_file/TruSightCancer_hg19_sorted_final.bed"
# opt$samples <- "all"
# opt$dir <- "/mnt/genetica/CIBERERBioinfo/CNVs/ICR96_dataset/results/applied_bqsr_data/"



output_dir <- paste(opt$outputdir,"/cnvs/",sep="")
bedFile <- opt$bed
projectname <- opt$name

setwd(output_dir)


#=====================#
#                     #
#   #=============#   #
#   # IMPORT DATA #   # 
#   #=============#   #
#                     #
#=====================#

#*****************#
# Import Bed File #
#*****************#
bedFile_data <- read.delim(bedFile, header = FALSE, stringsAsFactors = FALSE)

if(length(colnames(bedFile_data))<4){
  cat("ERROR: bed file without NAME field")
  quit(status = 1)} else if (length(colnames(bedFile_data))>4){
  bedFile_data = bedFile_data[,1:4]}
  
colnames(bedFile_data) <- c("CHR", "START", "STOP", "GENE_NAME")
bedFile_data$AVG_POS <- (bedFile_data$START + bedFile_data$STOP) / 2 # Create a column with the medium point of the exons
bedFile_data$CHR <- gsub("chr", "", bedFile_data$CHR)
#bedFile_data <- bedFile_data[!grepl("^rs.*$",bedFile_data$GENE_NAME),] # Delete SNP rows
#bedFile_data <- bedFile_data[!grepl("^.*_SNP$",bedFile_data$GENE_NAME),] # Delete SNP rows


# bedFile_data <- bedFile_data[order(bedFile_data$CHR, bedFile_data$START, decreasing = T),]
# zero_rows = which(bedFile_data$GENE_NAME == "0")
# 
# for (i in zero_rows){
# 
#   if ( i == 1 ){ # If it is the first
#     if ( bedFile_data[i+1, 4] != "0" ) { bedFile_data[i, 4] = bedFile_data[i+1, 4] }
#     else { bedFile_data[i, 4] = paste(bedFile_data[i, 1], bedFile_data[i, 2], sep = "_") }
#   }
#   else if (i == nrow(bedFile_data)) { # If it is the last
#     if ( bedFile_data[i-1, 4] != "0" ) { bedFile_data[i, 4] = bedFile_data[i-1, 4] }
#     else { bedFile_data[i, 4] = paste(bedFile_data[i, 1], bedFile_data[i, 2], sep = "_") }
#   }
#   else { # Beginning
#     if (bedFile_data[i, 1] != bedFile_data[i-1, 1]){
#       if (bedFile_data[i+1, 4] != "0") { bedFile_data[i, 4] = bedFile_data[i+1, 4] }
#       #else { bedFile_data[i, 4] = paste(bedFile_data[i, 1], bedFile_data[i, 2], sep = "_") }
#     }
#     else if (bedFile_data[i, 1] != bedFile_data[i+1, 1]){
#       if (bedFile_data[i-1, 4] != "0") { bedFile_data[i, 4] = bedFile_data[i-1, 4] }
#       #else { bedFile_data[i, 4] = paste(bedFile_data[i, 1], bedFile_data[i, 2], sep = "_") }
#     }
#     else {
#       if ( which.min(c(abs(bedFile_data[i-1, 3] - bedFile_data[i, 2]), abs(bedFile_data[i+1, 2] - bedFile_data[i, 3])) ) == 1){
#         bedFile_data[i, 4] = bedFile_data[i-1, 4]
#       } else {
#           if ( bedFile_data[i+1, 4] != "0" ) { bedFile_data[i, 4] = bedFile_data[i+1, 4] }
#           else { bedFile_data[i, 4] = bedFile_data[i-1, 4] }
#       }
#     }
#   }
# }


# Create exon names. It changes the order in case there are overlapping genes.
bedFile_data <- bedFile_data[order(bedFile_data$CHR, bedFile_data$GENE_NAME, bedFile_data$START),]
j <- 1
for (i in 1:nrow(bedFile_data)){
  bedFile_data$TARGET[i] <- paste(bedFile_data$GENE_NAME[i], j, sep = "_")
  j <- j + 1
  if (bedFile_data$GENE_NAME[i] != bedFile_data$GENE_NAME[i+1] & i != nrow(bedFile_data)){j <- 1}
}

bedFile_data <- bedFile_data[order(bedFile_data$CHR, bedFile_data$START),]




# Create total_data. It contains the output of the programs
total_data <- data.frame(stringsAsFactors = FALSE)

#************************#
# Import ExomeDepth data #
#************************#
exomedepth_output <- paste(output_dir, "/ExomeDepth/", 'exomedepth', '.results.txt', sep = '')
tryCatch(
  {
    exomedepth_data <- read.delim(exomedepth_output, stringsAsFactors = FALSE) # Get the information from ExomeDepth's output
    exomedepth_data$PROGRAM <- "ED"
    
    df_parcial <- data.frame(exomedepth_data$samplename, exomedepth_data$chromosome, exomedepth_data$start, exomedepth_data$end, exomedepth_data$CNV_type, 
                             exomedepth_data$PROGRAM, exomedepth_data$reads.ratio, NA, NA, NA, exomedepth_data$BF, NA, NA, NA, stringsAsFactors = FALSE)
    colnames(df_parcial) <- c("SAMPLE_NAME", "CHR", "START", "STOP", "CNV_TYPE", "PROGRAM", "ED_ratio", "CN_ratio", "C2_ratio", "PM_ratio", "ED_BF", "CN_ZSCORE", "CN_SHAPIRO.WILK", "C2_mBIC")
    
    total_data <- rbind(total_data, df_parcial)
    
  },
  error=function(e) print("There is no ExomeDepth output"), 
  warning=function(e) print("There is no ExomeDepth output")
)
  

#***********************#
# Import Convading data #
#***********************#
convading_output <- paste(output_dir, "/CoNVaDING/", 'CoNVaDING', '.results.txt', sep = '')
tryCatch(
  {
    convading_data <- read.delim(convading_output, stringsAsFactors = FALSE) # Get the information from CoNVaDING's output
    convading_data$PROGRAM <- "CN"
    
    df_parcial <- data.frame(convading_data$SAMPLE_NAME, convading_data$CHR,convading_data$START, convading_data$STOP, convading_data$ABBERATION, 
                             convading_data$PROGRAM, NA, convading_data$AUTO_RATIO, NA, NA, NA, convading_data$AUTO_ZSCORE, convading_data$SHAPIRO.WILK, NA, stringsAsFactors = FALSE)
    colnames(df_parcial) <- c("SAMPLE_NAME", "CHR", "START", "STOP", "CNV_TYPE", "PROGRAM", "ED_ratio", "CN_ratio", "C2_ratio", "PM_ratio", "ED_BF", "CN_ZSCORE", "CN_SHAPIRO.WILK", "C2_mBIC")

    total_data <- rbind(total_data, df_parcial)
    
  },
  error=function(e) print("There is no CoNVaDING output"), 
  warning=function(e) print("There is no CoNVaDING output")
)


#********************#
# Import Codex2 data #
#********************#
codex2_output <- paste(output_dir, "/CODEX2/", 'CODEX2', '.results.txt', sep = '')
tryCatch(
  {
    codex2_data <- read.delim(codex2_output, stringsAsFactors = FALSE)
    codex2_data$chr <- gsub("chr", "", codex2_data$chr)
    codex2_data$cnv <- toupper(codex2_data$cnv)
    codex2_data$PROGRAM <- "C2"
    
    df_parcial <- data.frame(codex2_data$sample_name, codex2_data$chr, codex2_data$st_bp, codex2_data$ed_bp,
                             codex2_data$cnv, codex2_data$PROGRAM, NA, NA, codex2_data$copy_no/2, NA, NA, NA, NA, codex2_data$mBIC, stringsAsFactors = FALSE) # mBIC: modified bayesian information criterion
    colnames(df_parcial) <- c("SAMPLE_NAME", "CHR", "START", "STOP", "CNV_TYPE", "PROGRAM", "ED_ratio", "CN_ratio", "C2_ratio", "PM_ratio", "ED_BF", "CN_ZSCORE", "CN_SHAPIRO.WILK", "C2_mBIC")
    
    total_data <- rbind(total_data, df_parcial)
    
  },
  error=function(e) print("There is no CODEX2 output"), 
  warning=function(e) print("There is no CODEX2 output")
)


#**************************#
# Import Panelcn.MOPS data #
#**************************#
panelMops_output <- paste(output_dir, "/Panelcn.MOPS/", 'panelcn.MOPS', '.results.txt', sep = '')
tryCatch(
  {
    panelMops_data <- read.delim(panelMops_output, stringsAsFactors = FALSE)
    panelMops_data <- panelMops_data[!grepl("lowQual",panelMops_data$lowQual),]
    panelMops_data$PROGRAM <- "PM"
    
    df_parcial <- data.frame(panelMops_data$Sample, panelMops_data$Chr, panelMops_data$Start, panelMops_data$End,
                             panelMops_data$CN, panelMops_data$PROGRAM, NA, NA, NA, panelMops_data$RC.norm/panelMops_data$medRC.norm, NA, NA, NA, NA, stringsAsFactors = FALSE)
    colnames(df_parcial) <- c("SAMPLE_NAME", "CHR", "START", "STOP", "CNV_TYPE", "PROGRAM", "ED_ratio", "CN_ratio", "C2_ratio", "PM_ratio", "ED_BF", "CN_ZSCORE", "CN_SHAPIRO.WILK", "C2_mBIC")
    
    total_data <- rbind(total_data, df_parcial)

  },
  error=function(e) print("There is no Panelcn.MOPS output"), 
  warning=function(e) print("There is no Panelcn.MOPS output")
)

total_data <- data.frame(lapply(total_data, function(x) if (is.factor(x)) as.character(x) else {x}), stringsAsFactors=FALSE)

total_data$SAMPLE_NAME <- gsub("\\..*","", total_data$SAMPLE_NAME, perl = TRUE)
total_data$SAMPLE_NAME <- gsub("_.*","", total_data$SAMPLE_NAME, perl = TRUE)


#===================#
#                   #
#   #===========#   #
#   # PROCESING #   # 
#   #===========#   #
#                   #
#===================#

#*****************#
# Exon assignment #
#*****************#
total_exons <- do.call("rbind", apply(bedFile_data, 1, function(x){
  
  positionRef = as.integer(x[5])
  chrRef = x[1]
  exonName = x[6]
  
  total_data_WE = total_data[total_data$CHR == chrRef & positionRef > total_data$START & positionRef < total_data$STOP,]
  
  if (nrow(total_data_WE) != 0){
    total_data_WE$EXON <- exonName
  }
  return(total_data_WE)
}))

total_exons$HOMOZYGOUS <- NA
total_exons$HOMOZYGOUS[grep("HOM", total_exons$CNV_TYPE)] <- "HOM" # Mark homozygous CNVs (CN and PM)
total_exons$CNV_TYPE <- gsub("HOM_", "", total_exons$CNV_TYPE) # Convert HOM-DUP/DEL to DUP/DEL
total_exons$EXON_MIX <- paste(total_exons$SAMPLE_NAME, total_exons$EXON, total_exons$CNV_TYPE, sep = "_") # Create a key for comparison
total_exons$GENE <- gsub("_[0-9]*$", "", total_exons$EXON)
total_exons$GENE_MIX <- paste(total_exons$SAMPLE_NAME, total_exons$GENE, total_exons$CNV_TYPE, sep = "_")


total_exons$PROGRAM_MIX <- paste(total_exons$PROGRAM, total_exons$EXON_MIX, sep = "_")
total_exons <- total_exons[!duplicated(total_exons$PROGRAM_MIX),] # Delete rows with duplicated $PROGRAM_MIX
total_exons$EXON_POS <- as.numeric(gsub("^.*_", "", total_exons$EXON))


#*****************#
# Output creation #
#*****************#
interval.maker <- function(numbers, intput_type){
  # Transform a list of numbers to intervals
  if (intput_type == "str"){num_vector <- as.numeric(strsplit(numbers, split = ", ")[[1]])
  }else if (intput_type == "vec") {num_vector <- numbers
  }else{return(print("intput_type not detected"))}
  
  num_vector <- unique(num_vector)
  
  if (length(num_vector) <= 1){return(numbers)}
  num_vector <- sort(num_vector)
  output_interval <- c("")
  for (i in 2:length(num_vector)){
    if (num_vector[i] - 1 == num_vector[i - 1]){
      if (tail(output_interval, 1) == "-"){
        if (i == length(num_vector)){
          output_interval <- c(output_interval, num_vector[i])
        }
      }else{
        if (i == length(num_vector)){
          output_interval <- c(output_interval, num_vector[i - 1], "-", num_vector[i])
        }else{
          output_interval <- c(output_interval, num_vector[i - 1], "-")
        }
      }
    }else{
      if (i == length(num_vector)){
        output_interval <- c(output_interval, num_vector[i - 1], ", ", num_vector[i])
      }else{
        output_interval <- c(output_interval, num_vector[i - 1], ", ")
      }
    }
  }
  output_interval <- paste(output_interval, collapse = "")
  if (intput_type == "vec") {output_interval <- strsplit(output_interval, split = ", ")[[1]]}
  
  return(output_interval)
}

total_list <- split(total_exons, total_exons$GENE_MIX)

data_out <- data.frame(stringsAsFactors = FALSE)

row_idx <- 1
for (i in total_list){
  
  intervals <- interval.maker(unique(i$EXON_POS), "vec")
  
  for (j in intervals){
    min <- as.numeric(gsub("-.*$", "", j))
    max <- as.numeric(gsub("^.*-", "", j))
    
    CNV = i[i$EXON_POS >= min & i$EXON_POS <= max,]
    
    data_out[row_idx, "SAMPLE"] <- CNV$SAMPLE_NAME[1]
    data_out[row_idx, "GENE"] <- CNV$GENE[1]
    data_out[row_idx, "CNV_TYPE"] <- CNV$CNV_TYPE[1]
    
    homo_exons_ED<- gsub("^.*_", "", CNV[CNV$ED_ratio < 0.1 | CNV$ED_ratio > 1.75, "EXON"])
    homo_exons_ED <- paste(homo_exons_ED[!is.na(homo_exons_ED)], collapse = ",")
    if (homo_exons_ED != ""){homo_exons_ED <- paste("ED(", homo_exons_ED,") ", sep = "")}
    homo_exons_C2<- gsub("^.*_", "", CNV[CNV$C2_ratio < 0.1 | CNV$C2_ratio > 1.75, "EXON"])
    homo_exons_C2 <- paste(homo_exons_C2[!is.na(homo_exons_C2)], collapse = ",")
    if (homo_exons_C2 != ""){homo_exons_C2 <- paste("C2(", homo_exons_C2,") ", sep = "")}
    homo_exons_CN<- gsub("^.*_", "", CNV[CNV$PROGRAM == "CN" & CNV$HOMOZYGOUS == "HOM", "EXON"])
    homo_exons_CN <- paste(homo_exons_CN[!is.na(homo_exons_CN)], collapse = ",")
    if (homo_exons_CN != ""){homo_exons_CN <- paste("CN(", homo_exons_CN,") ", sep = "")}
    homo_exons_PM<- gsub("^.*_", "", CNV[CNV$PROGRAM == "PM" & CNV$HOMOZYGOUS == "HOM", "EXON"])
    homo_exons_PM <- paste(homo_exons_PM[!is.na(homo_exons_PM)], collapse = ",")
    if (homo_exons_PM != ""){homo_exons_PM <- paste("PM(", homo_exons_PM,") ", sep = "")}
    homo_exons = paste("HOM: ", homo_exons_ED, homo_exons_CN, homo_exons_C2, homo_exons_PM, sep = "")
    if (homo_exons == "HOM: "){data_out[row_idx, "HOMO_HETEROZYGOUS"] <- "HET"
    }else{data_out[row_idx, "HOMO_HETEROZYGOUS"] <- homo_exons}
    
    data_out[row_idx, "CHR"] <- CNV$CHR[1]
    data_out[row_idx, "START"] <- CNV$START[which.min(CNV$START)]
    data_out[row_idx, "STOP"] <- CNV$STOP[which.max(CNV$STOP)]
    
    data_out[row_idx, "NUM_OF_TARGETS"] <- length(unique(gsub("^.*_","",CNV$EXON)))
    data_out[row_idx, "TARGETS"] <- paste(unique(gsub("^.*_","",CNV$EXON)), collapse = ", ")
    data_out[row_idx, "ED_TARGETS"] <- paste(gsub("^.*_", "", CNV[CNV$PROGRAM == "ED",]$EXON), collapse = ", ")
    data_out[row_idx, "CN_TARGETS"] <- paste(gsub("^.*_", "", CNV[CNV$PROGRAM == "CN",]$EXON), collapse = ", ")
    data_out[row_idx, "C2_TARGETS"] <- paste(gsub("^.*_", "", CNV[CNV$PROGRAM == "C2",]$EXON), collapse = ", ")
    data_out[row_idx, "PM_TARGETS"] <- paste(gsub("^.*_", "", CNV[CNV$PROGRAM == "PM",]$EXON), collapse = ", ")
    data_out[row_idx, "NUM_OF_PROGRAMS"] <- max(as.vector(table(CNV$EXON)))
  
    data_out[row_idx, "ED_ratio"] <- mean(CNV[CNV$PROGRAM == "ED",]$ED_ratio)
    data_out[row_idx, "CN_ratio"] <- mean(CNV[CNV$PROGRAM == "CN",]$CN_ratio)
    data_out[row_idx, "C2_ratio"] <- mean(CNV[CNV$PROGRAM == "C2",]$C2_ratio)
    data_out[row_idx, "PM_ratio"] <- mean(CNV[CNV$PROGRAM == "PM",]$PM_ratio)
    
    data_out[row_idx, "ED_BF"] <- mean(CNV[CNV$PROGRAM == "ED",]$ED_BF)
    data_out[row_idx, "CN_ZSCORE"] <- mean(CNV[CNV$PROGRAM == "CN",]$CN_ZSCORE)
    data_out[row_idx, "CN_SHAPIRO.WILK"] <- mean(CNV[CNV$PROGRAM == "CN",]$CN_SHAPIRO.WILK)
    data_out[row_idx, "C2_mBIC"] <- mean(CNV[CNV$PROGRAM == "C2",]$C2_mBIC)
    
    row_idx <- row_idx + 1
  }
}



data_out2 <- data_out

data_out2$TARGETS <- sapply(data_out2$TARGETS, function(x) interval.maker(x, "str"))
data_out2$ED_TARGETS <- sapply(data_out2$ED_TARGETS, function(x) interval.maker(x, "str"))
data_out2$CN_TARGETS <- sapply(data_out2$CN_TARGETS, function(x) interval.maker(x, "str")) 
data_out2$C2_TARGETS <- sapply(data_out2$C2_TARGETS, function(x) interval.maker(x, "str"))
data_out2$PM_TARGETS <- sapply(data_out2$PM_TARGETS, function(x) interval.maker(x, "str"))


# Remove lines that are not in the sample list
total_samples <- unique(unlist(lapply(list.files(opt$dir, pattern = ".bam"), function(x) strsplit(x, split = ".", fixed = T)[[1]][1])))
if(opt$samples!="all"){
  analyced_samples = intersect(total_samples, strsplit(opt$samples, split=",")[[1]])
  analyced_samples = gsub("_.*", "", analyced_samples)
  data_out2 = data_out2[data_out2$SAMPLE %in% analyced_samples,]
  # data_out3 <- data.frame(stringsAsFactors = FALSE)
  # for (i in analyced_samples){
  #   data_out_parcial <- data_out2[grepl(i,data_out2[,1]),]
  #   if (nrow(data_out_parcial) >= 1){
  #     data_out3 <- rbind(data_out3, data_out_parcial)
  #   }
  # }
  # data_out2 <- data_out3
  mysamples = analyced_samples
}else{mysamples = total_samples}

cat("#Analyzed samples: ", paste(mysamples, collapse = ", "), "\n", sep = "", file = paste(projectname, ".combined.txt", sep = ""))

# Order the output
if(nrow(data_out2)>0){
  data_out2 <- data_out2[order(data_out2$NUM_OF_PROGRAMS, decreasing = TRUE),]
  data_out2[data_out2 == ""] <- NA
}


cat("#HOMO_HETEROZYGOUS: HOMOZYGOUS/HETEROZYGOUS CNV. Only the homozygous CNVs are specyfied like this: PROGRAM(TARGET_DETECTED_BY_THE_PROGRAM), otherwise it was detected as heterozygous.
#NUM_OF_TARGETS: Total number of targets detected
#TARGETS: All targets detected
#ED_TARGETS: Targets detected by ExomeDepth
#CN_TARGETS: Targets detected by CoNVaDING
#C2_TARGETS: Targets detected by CODEX2
#PM_TARGETS: Targets detected by Panelcn.MOPS
#NUM_OF_PROGRAMS: Number of programs that have detected at least the same target. The table is sorted descendingly following this column
#ED_BF: Bayes Factor (ExomeDepth)
#CN_ZSORE: Z-score (CoNVaDING)
#CN_SHAPIRO.WILK: Shapiro-Wilk test (CoNVaDING)
#C2_mBIC: Modified Bayesian Information Criterion (CODEX2)\n",
file = paste(projectname, ".combined.txt", sep = ""), append = TRUE)


# Final table
write.table(data_out2, file = paste(projectname, ".combined.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE)

# To Annotate table
#data_out_toAnnotate <- data.frame(data_out2$CHR, data_out2$START, data_out2$STOP, data_out2$CNV_TYPE, stringsAsFactors = FALSE)
#data_out_toAnnotate <-unique(data_out_toAnnotate)
#write.table(data_out_toAnnotate, file = paste(projectname, "_combined_toAnnotate.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Extended bed File
bedout <- bedFile_data
bedout$AVG_POS <- NULL
write.table(bedout, file = paste(projectname, ".extended.bed", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)



