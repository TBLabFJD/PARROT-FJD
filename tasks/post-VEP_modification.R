# Author: Gonzalo Nunez Moreno
# Date: 16/03/2021

#rm(list=ls()) 


library(optparse)




#############
# Arguments # 
#############
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="\t\t input TSV file (output from VEP)", metavar="character"),
  
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="\t\t Output file", metavar="character"),
  
  make_option(c("-d", "--dbNSFPgene"), type="character", default=NULL, 
              help="\t\tdbNSFP_gene file", metavar="character"),
  
  make_option(c("-r", "--regiondict"), type="character", default=NULL,
              help="\t\tRegion dictionary", metavar="character"),
  
  make_option(c("-m", "--omim"), type="character", default=NULL,
              help="\t\tOMIM information", metavar="character"),
  
  make_option(c("-n", "--numheader"), type="integer", default=1,
              help="\t\tNumber of the row where the header is", metavar="character"),
  
  make_option(c("-a", "--automap"), type="character", default=NULL,
              help="\t\tAutomap output directory", metavar="character"),
  
  make_option(c("-f", "--maf"), type="double", default=0.1,
              help="\t\tMinimum allele frequency to filter", metavar="character"),
  
  make_option(c("-g", "--genefilter"), type="character", default=NULL,
              help="\t\tGene list to filter the resutls", metavar="character"),
  
  make_option(c("-w", "--glowgenes"), type="character", default=NULL,
              help="\t\tGLOWgenes output file to annotate and srt the results", metavar="character")
 )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


input = opt$input
output = opt$output
dbNSFPgenepath <- opt$dbNSFPgene
dict_region_path <- opt$regiondict
omim_path = opt$omim
skip = opt$numheader
automap_path = opt$automap
maf = opt$maf
genefilter_path = opt$genefilter
glowgenes_path = opt$glowgenes

# 
# input = "/home/gonzalo/tblab/mnt/genetica4/gonzalo/trio_externo_12oct/trio.annotated.1.tsv"
# output = "/home/gonzalo/tblab/mnt/genetica4/gonzalo/trio_externo_12oct/trio.annotated.pvm.tsv"
# dbNSFPgenepath = "/home/gonzalo/tblab/mnt/genetica7/vep_plugins_gene/dbNSFP4.3_gene.complete.pvm.txt"
# dict_region_path = "/home/gonzalo/tblab/mnt/genetica7/vep_plugins_gene/dict_region.csv"
# omim_path = "/home/gonzalo/tblab/mnt/genetica7/vep_plugins_gene/omim_genemap2.txt"
# skip = 166
# automap_path = "/home/gonzalo/tblab/mnt/genetica4/gonzalo/trio_externo_12oct/"
# maf = 0.1
# genefilter_path = NULL
# glowgenes_path = NULL
# 
# input = "/home/gonzalo/tblab/mnt/tblab/raquel/reanalysis/mCorton/results/sis-2355498526-49-21855776-full.vep.tsv"
# output = "/home/gonzalo/tblab/mnt/tblab/raquel/reanalysis/mCorton/results/sis-2355498526-49-21855776-full.pvm.tsv"
# dbNSFPgenepath = "/home/gonzalo/tblab/mnt/genetica7/vep_plugins_gene/dbNSFP4.3_gene.complete.pvm.txt"
# dict_region_path = "/home/gonzalo/tblab/mnt/genetica7/vep_plugins_gene/dict_region.csv"
# omim_path = "/home/gonzalo/tblab/mnt/genetica7/vep_plugins_gene/omim_genemap2.txt"
# skip = 162
# automap_path = "/home/gonzalo/tblab/mnt/genetica6/reanotacion/results/"
# maf = 0.1
# genefilter_path = "/home/gonzalo/tblab/mnt/genetica7/distrofias_retina_sindr_nosindr.txt"
# glowgenes_path = "/home/gonzalo/tblab/mnt/genetica7/GLOWgenes_prioritization_Random.txt"


################
# Data loading # 
################

# VEP
vep = read.delim(input, header = TRUE, skip = skip-1, stringsAsFactors = F, quote = "", check.names=F, colClasses = "character")
#copy_vep=vep
#vep=copy_vep
# dbNSFP gene
dbNSFP_gene = read.delim(dbNSFPgenepath, header = TRUE, stringsAsFactors = F, quote = "")
vep = merge(vep, dbNSFP_gene, by.x = "SYMBOL", by.y = "Gene_name", all.x = T)


# OMIM
if (!is.null(omim_path)){
  omim = read.delim(omim_path, header = F, stringsAsFactors = F, comment.char = "#", quote = "", check.names=F)
  colnames(omim) = c("Chromosome", "Genomic_Position_Start", "Genomic Position End", "Cyto_Location", "Computed_Cyto_Location", "MIM_Number",
    "Gene_Symbols", "Gene_Name",	"Approved_Gene_Symbol", "Entrez_Gene_ID", "Ensembl_Gene_ID", "Comments", "Phenotypes", "Mouse_Gene_Symbol-ID")
  vep = merge(vep, omim, by.x = "SYMBOL", by.y = "Approved_Gene_Symbol", all.x = T)
}


# Region dictionary
if (!is.null(dict_region_path)){
  dict_region = read.csv(dict_region_path, header = F, sep = ",", stringsAsFactors = F)
  priority_list = c("SPLICING", "5UTR", "3UTR", "ncRNA", "regulatory", "UPSTREAM", "DOWNSTREAM", "EXONIC", "INTRONIC", "INTERGENIC", "-")
  dict_region$V2[is.na(dict_region$V2)] = "-"
  dict_region$V2 = factor(dict_region$V2, priority_list)
  rownames(dict_region) = dict_region$V1
}


# genefilter
if (!is.null(genefilter_path)){
  genefilter = read.delim(genefilter_path, header = F, stringsAsFactors = F, quote = "", check.names=F)
}


# glowgenes
if (!is.null(glowgenes_path)){
  
  glowgenes = read.delim(glowgenes_path, header = F, stringsAsFactors = F, quote = "", check.names=F)
  colnames(glowgenes) = c("SYMBOL", "score", "GLOWgenes")
  
  # if (!is.null(genefilter_path)){
  #   genefilter$score = NA
  #   genefilter$GLOWgenes = 0
  #   colnames(genefilter) = c("SYMBOL", "score", "GLOWgenes")
  # 
  #   glowgenes = rbind(genefilter, glowgenes)
  # }
  
  vep = merge(vep, glowgenes, by.x = "SYMBOL", by.y = "SYMBOL", all.x = T)
}


# copy_vep = vep

#===========#
# Filtering #
#===========#
# AF filtering
print(nrow(vep))
vep$gnomADe_AF_popmax = as.numeric(unlist(lapply(vep$gnomADe_AF_popmax, function(x) strsplit(x, ",")[[1]][1])))
vep$gnomADg_AF_popmax = as.numeric(unlist(lapply(vep$gnomADg_AF_popmax, function(x) strsplit(x, ",")[[1]][1])))

vep = vep[is.na(vep$gnomADe_AF_popmax) | as.numeric(vep$gnomADe_AF_popmax) < as.numeric(maf) | vep$gnomADe_filt != "PASS",]
print(nrow(vep))
vep = vep[is.na(vep$gnomADg_AF_popmax) | as.numeric(vep$gnomADg_AF_popmax) < as.numeric(maf) | vep$gnomADg_filt != "PASS",]
print(nrow(vep))

# Gene Filter
if ((!is.null(genefilter_path)) & (is.null(glowgenes_path))){
  vep = vep[vep$SYMBOL %in% genefilter$V1,]
}
# print("HOLA")
# save.image(file = 'PVM.RData')

if (nrow(vep) == 0) {
  stop("There are no remaining variants")
}

# [1] 54242
# [1] 54349

# vep[is.na(vep)] = "-"
# vep[vep=="-"] = ""
df_out  = data.frame(row.names = 1:nrow(vep), stringsAsFactors = F)
  




#==================================#
# Basic information of the variant #
#==================================#
print("Basic information of the variant")

df_out$CHROM = unlist(lapply(vep$Location, function(x) strsplit(x, ":")[[1]][1]))
df_out$POS = as.numeric(unlist(lapply(vep$`#Uploaded_variation`, function(x) rev(strsplit(x, "_")[[1]])[2])))
df_out$REF = vep$USED_REF
df_out$ALT = vep$Allele
df_out$Location = vep$Location
df_out$SYMBOL = vep$SYMBOL
df_out$Gene_full_name = vep$Gene_full_name
if (!is.null(glowgenes_path)) df_out$GLOWgenes = vep$GLOWgenes
if ((!is.null(genefilter_path)) & (!is.null(glowgenes_path))) df_out$GLOWgenes[df_out$SYMBOL %in% genefilter$V1] = 0 # We assume the genes of the list are the genes from the panel
df_out$VARIANT_CLASS = vep$VARIANT_CLASS





#====================#
# Sample information #
#====================#
print("Sample information")
# 
# samples = gsub("_GT$$", "", gsub("^SAMPLE_", "", colnames(vep)[grepl(".*_GT$", colnames(vep), perl = T)])) # Get
# alleles = unlist(lapply(vep$`#Uploaded_variation`, function(x) strsplit(x, "_")[[1]][3])) # Get all posible alleles (including reference)
# alleles_df = data.frame(alleles, vep$Allele, row.names = 1:nrow(vep), stringsAsFactors = F)
# vep$alleles_pos = as.numeric(unlist(apply(alleles_df,1, function(x) which(strsplit(x[1], "/")[[1]] == x[2]))))
# for (sample in samples){
#   
#   tryCatch(
#     {
#       df_out[,paste0(sample,"_GT")] = vep[,paste0("SAMPLE_", sample,"_GT")]
#     },
#     error=function(e) print(paste0("There is no GT information of the sample ", sample)), 
#     warning=function(e) print(paste0("There is no GT information of the sample ", sample))
#   )
#   
#   tryCatch(
#     {
#       AD = as.numeric(apply(vep, 1, function(x) strsplit(x[paste0("SAMPLE_", sample,"_AD")], split = "_")[[1]][as.numeric(x["alleles_pos"])]))
#     },
#     error=function(e) print(paste0("There is no AD information of the sample ", sample)), 
#     warning=function(e) print(paste0("There is no AD information of the sample ", sample))
#   )
#   
#   tryCatch(
#     {
#       df_out[,paste0(sample,"_VF")] = round(AD/as.numeric(vep[,paste0("SAMPLE_", sample,"_DP")]),2)
#       
#     },
#     error=function(e) print(paste0("There is no AD or DP information of the sample ", sample)), 
#     warning=function(e) print(paste0("There is no AD or DP information of the sample ", sample))
#   )
#   
#   tryCatch(
#     {
#       df_out[,paste0(sample,"_AD")] = vep[,paste0("SAMPLE_", sample,"_AD")] 
#     },
#     error=function(e) print(paste0("There is no AD information of the sample ", sample)), 
#     warning=function(e) print(paste0("There is no AD information of the sample ", sample))
#   )
#   
#   tryCatch(
#     {
#       df_out[,paste0(sample,"_DP")] = as.numeric(vep[,paste0("SAMPLE_", sample,"_DP")])
#     },
#     error=function(e) print(paste0("There is no DP information of the sample ", sample)), 
#     warning=function(e) print(paste0("There is no DP information of the sample ", sample))
#   )
#     
#   tryCatch(
#     {
#       df_out[,paste0(sample,"_GQ")] = as.numeric(vep[,paste0("SAMPLE_", sample,"_GQ")])
#     },
#     error=function(e) print(paste0("There is no GQ information of the sample ", sample)), 
#     warning=function(e) print(paste0("There is no GQ information of the sample ", sample))
#   )
#   
#   # Sacar del output de autopmap
#   tryCatch(
#     {
#       # automap = read.delim(paste0(automap_path,"/", sample, "/", sample, ".HomRegions.tsv"), header = F, comment.char = "#", stringsAsFactors = F)
#       automap = read.delim(paste0(automap_path, "/", sample, ".HomRegions.tsv"), header = F, comment.char = "#", stringsAsFactors = F)
#       df_out[,paste0(sample,"_ROH")] = "False"
#       for (i in 1:nrow(automap)) df_out[df_out$POS >= automap$V2[i] & df_out$POS <= automap$V3[i] & gsub("chr","",df_out$CHROM) == gsub("chr","",automap$V1[i]), paste0(sample,"_ROH")] = "True"
#     },
#     error=function(e) print(paste0("There is no AutoMap information of the sample ", sample)), 
#     warning=function(e) print(paste0("There is no AutoMap information of the sample ", sample))
#   )
# }

samples = unique(gsub("_.*$", "", gsub("^SAMPLE_", "", colnames(vep)[grepl(".*_GT$", colnames(vep), perl = T)])))
for (sample in samples){
  for (field in c("GT", "VAF", "AD", "DP", "SF", "GD")){
    tryCatch(
      {
        print(paste0(sample, "_", field))
        df_out[,paste0(sample, "_", field)] = vep[,paste0("SAMPLE_", sample, "_", field)]
      },
      error=function(e) print(paste0("There is no ", field, " information of the sample ", sample)),
      warning=function(e) print(paste0("There is no ", field, " information of the sample ", sample))
      )
  }
  
  # Sacar del output de autopmap
  tryCatch(
    {
      # automap = read.delim(paste0(automap_path,"/", sample, "/", sample, ".HomRegions.tsv"), header = F, comment.char = "#", stringsAsFactors = F)
      automap = read.delim(paste0(automap_path, "/", sample, ".HomRegions.tsv"), header = F, comment.char = "#", stringsAsFactors = F)
      df_out[,paste0(sample,"_ROH")] = "False"
      for (i in 1:nrow(automap)) df_out[df_out$POS >= automap$V2[i] & df_out$POS <= automap$V3[i] & gsub("chr","",df_out$CHROM) == gsub("chr","",automap$V1[i]), paste0(sample,"_ROH")] = "True"
    },
    error=function(e) print(paste0("There is no AutoMap information of the sample ", sample)), 
    warning=function(e) print(paste0("There is no AutoMap information of the sample ", sample))
  )
}
df_out$hiConfDeNovo = vep$SAMPLE_hiConfDeNovo
df_out$loConfDeNovo = vep$SAMPLE_loConfDeNovo



df_out$hiConfDeNovo = vep$SAMPLE_hiConfDeNovo
df_out$loConfDeNovo = vep$SAMPLE_loConfDeNovo





#=====================#
# Feature information #
#=====================#
print("Feature information")

df_out$Existing_variation = vep$Existing_variation
if (!is.null(dict_region_path)) df_out$Genomic_region = unlist(lapply(vep$Consequence, function(x) as.character(dict_region[strsplit(x, ",")[[1]],2])[which.min(dict_region[strsplit(x, ",")[[1]],2])]))
df_out$CANONICAL = vep$CANONICAL
df_out$Feature = vep$Feature
df_out$Feature_type = vep$Feature_type
df_out$BIOTYPE = vep$BIOTYPE
df_out$Consequence = vep$Consequence
df_out$INTRON = vep$INTRON
df_out$EXON = vep$EXON
df_out$HGVSc = vep$HGVSc
df_out$HGVSp = vep$HGVSp
df_out$DISTANCE = as.numeric(vep$DISTANCE)
df_out$STRAND = vep$STRAND
df_out$Interpro_domain = vep$Interpro_domain




#===============#
# Pathogenicity #
#===============#
print("Pathogenicity")

df_out$CLNSIG = vep$ClinVar_CLNSIG
df_out$CLNREVSTAT = vep$ClinVar_CLNREVSTAT
df_out$CLNDN = vep$ClinVar_CLNDN
df_out$OMIM_phenotype = vep$Phenotypes
df_out$Orphanet_disorder = vep$Orphanet_disorder
df_out$Orphanet_association_type = vep$Orphanet_association_type
df_out$HPO_name = vep$HPO_name
df_out$PUBMED = vep$PUBMED





#=============#
# Frequencies #
#=============#
print("Frequencies")

df_out$gnomADg_AF =  as.numeric(unlist(lapply(vep$gnomADg_AF, function(x) strsplit(x, ",")[[1]][1])))
df_out$gnomADg_AC = as.numeric(unlist(lapply(vep$gnomADg_AC, function(x) strsplit(x, ",")[[1]][1]))) 
df_out$gnomADg_AN = as.numeric(unlist(lapply(vep$gnomADg_AN, function(x) strsplit(x, ",")[[1]][1])))
df_out$gnomADg_nhomalt = as.numeric(unlist(lapply(vep$gnomADg_nhomalt, function(x) strsplit(x, ",")[[1]][1])))
df_out$gnomADg_cov_median = round(unlist(lapply(vep$gnomADg_cov_median, function(x) mean(as.numeric(strsplit(gsub("(?<![eE])-","0",x, perl = T), ",")[[1]])))))
df_out$gnomADg_cov_perc_20x = round(unlist(lapply(vep$gnomADg_cov_perc_20x, function(x) mean(as.numeric(strsplit(gsub("(?<![eE])-","0",x, perl = T), ",")[[1]])))),2)
df_out$gnomADg_filter = vep$gnomADg_filt
df_out$gnomADg_popmax = vep$gnomADg_popmax
df_out$gnomADg_AF_popmax = vep$gnomADg_AF_popmax
df_out$gnomADg_AC_popmax = as.numeric(unlist(lapply(vep$gnomADg_AC_popmax, function(x) strsplit(x, ",")[[1]][1])))
df_out$gnomADg_AF_nfe = as.numeric(unlist(lapply(vep$gnomADg_AF_nfe, function(x) strsplit(x, ",")[[1]][1])))
df_out$gnomADg_AC_nfe = as.numeric(unlist(lapply(vep$gnomADg_AC_nfe, function(x) strsplit(x, ",")[[1]][1])))

df_out$gnomADe_AF = as.numeric(unlist(lapply(vep$gnomADe_AF, function(x) strsplit(x, ",")[[1]][1])))
df_out$gnomADe_AC = as.numeric(unlist(lapply(vep$gnomADe_AC, function(x) strsplit(x, ",")[[1]][1])))
df_out$gnomADe_AN = as.numeric(unlist(lapply(vep$gnomADe_AN, function(x) strsplit(x, ",")[[1]][1])))
df_out$gnomADe_nhomalt = as.numeric(unlist(lapply(vep$gnomADe_nhomalt, function(x) strsplit(x, ",")[[1]][1])))
df_out$gnomADe_cov_median = round(unlist(lapply(vep$gnomADe_cov_median, function(x) mean(as.numeric(strsplit(gsub("(?<![eE])-","0",x, perl = T), ",")[[1]])))))
df_out$gnomADe_cov_perc_20x = round(unlist(lapply(vep$gnomADe_cov_perc_20x, function(x) mean(as.numeric(strsplit(gsub("(?<![eE])-","0",x, perl = T), ",")[[1]])))),2)
df_out$gnomADe_filter = vep$gnomADe_filt
df_out$gnomADe_popmax = vep$gnomADe_popmax
df_out$gnomADe_AF_popmax = vep$gnomADe_AF_popmax
df_out$gnomADe_AC_popmax = as.numeric(unlist(lapply(vep$gnomADe_AC_popmax, function(x) strsplit(x, ",")[[1]][1])))
df_out$gnomADe_AF_nfe = as.numeric(unlist(lapply(vep$gnomADe_AF_nfe, function(x) strsplit(x, ",")[[1]][1])))
df_out$gnomADe_AC_nfe = as.numeric(unlist(lapply(vep$gnomADe_AC_nfe, function(x) strsplit(x, ",")[[1]][1])))

df_out$kaviar_AF = vep$kaviar_AF
df_out$kaviar_AC = vep$kaviar_AC
df_out$CSVS_AF = as.numeric(unlist(lapply(vep$CSVS_AF, function(x) strsplit(x, ",")[[1]][1])))
df_out$CSVS_AC = as.numeric(unlist(lapply(vep$CSVS_AC, function(x) strsplit(x, ",")[[1]][1])))
df_out$FJD_MAF_AF = as.numeric(unlist(lapply(vep$FJD_MAF_AF, function(x) strsplit(x, ",")[[1]][1])))
df_out$FJD_MAF_AC = as.numeric(unlist(lapply(vep$FJD_MAF_AC, function(x) strsplit(x, ",")[[1]][1])))
df_out$denovoVariants_SAMPLE_CT = vep$denovoVariants_SAMPLE_CT





#==========================#
# Pathogenicity prediction #
#==========================#
print("Pathogenicity prediction")

df_out$CADD_PHRED = as.numeric(vep$CADD_PHRED)
df_out$CADD_RAW = as.numeric(vep$CADD_RAW)
df_out$MutScore = as.numeric(vep$Mut_Score)

patho_norm_func = function(predictions){
  predictions = gsub(";", ",", predictions)
  predictions = tolower(predictions)
  unlist(lapply(predictions, function(x) {
    y = strsplit(x,",")[[1]]
    y = y[!y %in% c("-", ".", "U", "")]
    y[y %in% c("tolerated", "tolerated_low_confidence", "benign", "n", "l", "p", "t")] = "T"
    y[y %in% c("deleterious", "deleterious_low_confidence", "probably_damaging", 
               "possibly_damaging", "a", "m", "h", "d", "Dominant", "Recessive")] = "D"
    if ("D" %in% y) { return("D") }
    else if ("T" %in% y) { return("T") }
    else {return("")}
  }))
}  

df_pathogenic_predictors = data.frame(row.names = 1:nrow(vep))
df_pathogenic_predictors$SIFT = patho_norm_func(vep$SIFT)
df_pathogenic_predictors$PolyPhen = patho_norm_func(vep$PolyPhen)
df_pathogenic_predictors$Polyphen2_HDIV_pred = patho_norm_func(vep$Polyphen2_HDIV_pred)
df_pathogenic_predictors$Polyphen2_HVAR_pred = patho_norm_func(vep$Polyphen2_HVAR_pred)
df_pathogenic_predictors$LRT_pred = patho_norm_func(vep$LRT_pred)
df_pathogenic_predictors$`M-CAP_pred` = patho_norm_func(vep$`M-CAP_pred`)
df_pathogenic_predictors$MetaLR_pred = patho_norm_func(vep$MetaLR_pred)
df_pathogenic_predictors$MetaSVM_pred = patho_norm_func(vep$MetaSVM_pred)
df_pathogenic_predictors$MutationAssessor_pred = patho_norm_func(vep$MutationAssessor_pred)
df_pathogenic_predictors$MutationTaster_pred = patho_norm_func(vep$MutationTaster_pred)
df_pathogenic_predictors$PROVEAN_pred = patho_norm_func(vep$PROVEAN_pred)
df_pathogenic_predictors$FATHMM_pred = patho_norm_func(vep$FATHMM_pred)
df_pathogenic_predictors$MetaRNN_pred = patho_norm_func(vep$MetaRNN_pred)
df_pathogenic_predictors$PrimateAI_pred = patho_norm_func(vep$PrimateAI_pred)
df_pathogenic_predictors$DEOGEN2_pred = patho_norm_func(vep$DEOGEN2_pred)
df_pathogenic_predictors$BayesDel_addAF_pred = patho_norm_func(vep$BayesDel_addAF_pred)
df_pathogenic_predictors$BayesDel_noAF_pred = patho_norm_func(vep$BayesDel_noAF_pred)
df_pathogenic_predictors$ClinPred_pred = patho_norm_func(vep$ClinPred_pred)
df_pathogenic_predictors$`LIST-S2_pred` = patho_norm_func(vep$`LIST-S2_pred`)
df_pathogenic_predictors$Aloft_pred = patho_norm_func(vep$Aloft_pred)
df_pathogenic_predictors$`fathmm-MKL_coding_pred` = patho_norm_func(vep$`fathmm-MKL_coding_pred`)
df_pathogenic_predictors$`fathmm-XF_coding_pred` = patho_norm_func(vep$`fathmm-XF_coding_pred`)
  

df_out$N_Pathogenic_pred = apply(df_pathogenic_predictors, 1, function(x) table(x)["D"])
df_out$N_Pathogenic_pred[is.na(df_out$N_Pathogenic_pred)] = 0
df_out$N_Benign_pred = apply(df_pathogenic_predictors, 1, function(x) table(x)["T"])
df_out$N_Benign_pred[is.na(df_out$N_Benign_pred)] = 0
df_out$N_predictions = df_out$N_Pathogenic_pred + df_out$N_Benign_pred
df_out$Pathogenic_pred = apply(df_pathogenic_predictors, 1, function(x) paste(names(x)[which(x == "D")], collapse = ","))
df_out$Benign_pred = apply(df_pathogenic_predictors, 1, function(x) paste(names(x)[which(x == "T")], collapse = ","))





#=====================#
# Splicing predictors #
#=====================#
print("Splicing predictors")

# Select one splice prediction per row
for (j in c("SpliceAI_SNV_SpliceAI", "SpliceAI_INDEL_SpliceAI")){
  multi_gene_sites = grep(",",vep[,j])
  for (i in multi_gene_sites){
    splice_predictions = do.call("rbind",strsplit(strsplit(vep[i,j], ",")[[1]], "|",fixed = T))
    if (vep$SYMBOL[i] %in% splice_predictions[,2]) {
      vep[i,j] = paste(splice_predictions[splice_predictions[,2] == vep$SYMBOL[i],,drop = F][1,],collapse = "|")
    } else {
      max_value_row = which(splice_predictions[,3:6] == max(splice_predictions[,3:6]), arr.ind = TRUE)[1,1]
      vep[i,j] = paste(splice_predictions[max_value_row,],collapse = "|")
    }
  }
}
# Merge SpliceAI predictions for INDELs and SNVs 
vep$SpliceAI_INDEL_SpliceAI[vep$SpliceAI_INDEL_SpliceAI == "-"] = vep$SpliceAI_SNV_SpliceAI[vep$SpliceAI_INDEL_SpliceAI == "-"]
# Create data.frame with separated SpliceAI values
SpliceAI = data.frame(do.call("rbind", strsplit(vep$SpliceAI_INDEL_SpliceAI, "|", fixed = T)), stringsAsFactors = F)
colnames(SpliceAI) = c("ALLELE", "SYMBOL", "DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL")
df_out$SpliceAI_SYMBOL = SpliceAI$SYMBOL
df_out$SpliceAI_DS_AG = as.numeric(SpliceAI$DS_AG)
df_out$SpliceAI_DS_AL = as.numeric(SpliceAI$DS_AL)
df_out$SpliceAI_DS_DG = as.numeric(SpliceAI$DS_DG)
df_out$SpliceAI_DS_DL = as.numeric(SpliceAI$DS_DL)
df_out$SpliceAI_DP_AG = as.numeric(SpliceAI$DP_AG)
df_out$SpliceAI_DP_AL = as.numeric(SpliceAI$DP_AL)
df_out$SpliceAI_DP_DG = as.numeric(SpliceAI$DP_DG)
df_out$SpliceAI_DP_DL = as.numeric(SpliceAI$DP_DL)

df_out$ada_score = as.numeric(vep$ada_score)
df_out$rf_score = as.numeric(vep$rf_score)
df_out$MaxEntScan_alt = as.numeric(vep$MaxEntScan_alt)
df_out$MaxEntScan_diff = as.numeric(vep$MaxEntScan_diff)
df_out$MaxEntScan_ref = as.numeric(vep$MaxEntScan_ref)





#============================#
# Conservation and phylogeny #
#============================#
print("Conservation and phylogeny")

df_out$LoFtool = as.numeric(vep$LoFtool)
df_out$ExACpLI = as.numeric(vep$ExACpLI)
df_out$gnomAD_exomes_CCR = vep$gnomAD_exomes_CCR
df_out$phastCons30way_mammalian = as.numeric(vep$phastCons30way_mammalian)
df_out$phyloP30way_mammalian = as.numeric(vep$phyloP30way_mammalian)
df_out$MGI_mouse_phenotype = vep$MGI_mouse_phenotype_filt





#=====================================================#
#Expression, process, route, function and interaction #
#=====================================================#
print("Conservation and phylogeny")

df_out$GTEx_V8_gene = vep$GTEx_V8_gene
df_out$GTEx_V8_tissue = vep$GTEx_V8_tissue
df_out$`Expression_GNF-Atlas` = vep$Expression.GNF.Atlas.
df_out$Pathway_KEGG = vep$Pathway.KEGG._full
df_out$GO_biological_process = vep$GO_biological_process
df_out$GO_cellular_component = vep$GO_cellular_component
df_out$GO_molecular_function = vep$GO_molecular_function
df_out$Interactions_IntAct = vep$Interactions.IntAct.


df_out$Original_pos = vep$SAMPLE_Original_pos
df_out$variant_id = vep$SAMPLE_variant_id





#==============================================#
#Extra sample information (individual callers) #
#==============================================#
for (sample in samples){
  for (field in c("DV_GT", "DV_DP", "DV_VD", "DR_GT", "DR_DP", "DR_VD", "GK_GT", "GK_DP", "GK_VD")){
    tryCatch(
      {
        print(paste0(sample, "_", field))
        df_out[,paste0(sample, "_", field)] = vep[,paste0("SAMPLE_", sample, "_", field)]
      },
      error=function(e) print(paste0("There is no ", field, " information of the sample ", sample)),
      warning=function(e) print(paste0("There is no ", field, " information of the sample ", sample))
    )
  }
} 




#======#
# Sort #
#======#
if (!is.null(glowgenes_path)){
  df_out = df_out[order(as.numeric(df_out$GLOWgenes), df_out$CHROM, df_out$POS),]
}else{
  df_out = df_out[order(df_out$CHROM, df_out$POS),]
}



#==============#
# Write output #
#==============#
df_out[df_out=="-"] = NA
write.table(df_out, output, sep = "\t", col.names = T, row.names = F, quote = F, na = "")

library(openxlsx)
write.xlsx(df_out, paste0(output, ".xlsx"), colNames = TRUE)

