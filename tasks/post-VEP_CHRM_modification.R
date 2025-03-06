# Author: Yolanda Benítez Quesada (based on tasks/post-AnnotSV_modification.R)
# Date: 24/02/2025

#rm(list=ls()) 


library(optparse)
library(data.table)
library(dplyr)
library(tidyr)



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
  
  make_option(c("-s", "--domino"), type="character", default=NULL, 
              help="\t\tdomino file", metavar="character"),

  make_option(c("-e", "--expression"), type="character", default=NULL, 
            help="\t\ttissue expression file", metavar="character"),
  
  make_option(c("-n", "--numheader"), type="integer", default=1,
              help="\t\tNumber of the row where the header is", metavar="character"),
  
  make_option(c("-a", "--automap"), type="character", default=NULL,
              help="\t\tAutomap output directory", metavar="character"),
  
  make_option(c("-g", "--genefilter"), type="character", default=NULL,
              help="\t\tGene list to filter the resutls", metavar="character"),
  
  make_option(c("-w", "--glowgenes"), type="character", default=NULL,
              help="\t\tGLOWgenes output file to annotate and srt the results", metavar="character"),
  
  make_option(c("-p", "--panels"), type="character", default=NULL,
              help="\t\tGene-Panel file to annotate", metavar="character")
 )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


input = opt$input
output = opt$output
dbNSFPgenepath <- opt$dbNSFPgene
dominopath <- opt$domino
expressionpath <- opt$expression
dict_region_path <- opt$regiondict
omim_path = opt$omim
skip = opt$numheader
automap_path = opt$automap
genefilter_path = opt$genefilter
glowgenes_path = opt$glowgenes
panels_path = opt$panels

# 
# input = "/home/gonzalo/tblab/mnt/genetica4/gonzalo/trio_externo_12oct/trio.annotated.1.tsv"
# output = "/home/gonzalo/tblab/mnt/genetica4/gonzalo/trio_externo_12oct/trio.annotated.pvm.tsv"
# dbNSFPgenepath = "/home/gonzalo/tblab/mnt/genetica7/vep_plugins_gene/dbNSFP4.3_gene.complete.pvm.txt"
# dict_region_path = "/home/gonzalo/tblab/mnt/genetica7/vep_plugins_gene/dict_region.csv"
# omim_path = "/home/gonzalo/tblab/mnt/genetica7/vep_plugins_gene/omim_genemap2.txt"
# skip = 166
# automap_path = "/home/gonzalo/tblab/mnt/genetica4/gonzalo/trio_externo_12oct/"
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
# genefilter_path = "/home/gonzalo/tblab/mnt/genetica7/distrofias_retina_sindr_nosindr.txt"
# glowgenes_path = "/home/gonzalo/tblab/mnt/genetica7/GLOWgenes_prioritization_Random.txt"


################
# Data loading # 
################

# VEP
vep = read.delim(input, header = TRUE, skip = skip-1, stringsAsFactors = F, quote = "", check.names=F, colClasses = "character")

# Split the rows
vep_split <- vep %>% tidyr::separate_rows(Mitomap_genomeloci, sep = ",")
#print(nrow(vep))
#head(vep_split$Mitomap_genomeloci)
#print(nrow(vep))

vep <- vep_split
#copy_vep=vep
#vep=copy_vep
# dbNSFP gene
dbNSFP_gene = read.delim(dbNSFPgenepath, header = TRUE, stringsAsFactors = F, quote = "")
#vep = merge(vep, dbNSFP_gene, by.x = "SYMBOL", by.y = "Gene_name", all.x = T)
vep = merge(vep, dbNSFP_gene, by.x = "Mitomap_genomeloci", by.y = "Gene_name", all.x = T)

# domino: add Graci
domino = read.delim(dominopath, header = TRUE, stringsAsFactors = F, quote = "")
vep = merge(vep, domino, by.x = "Mitomap_genomeloci", by.y = "Gene_name", all.x = T)

# tissue expression: add Yoli
expression = read.delim(expressionpath, header = TRUE, stringsAsFactors = F, quote = "")
vep = merge(vep, expression, by.x = "Mitomap_genomeloci", by.y = "Gene.name", all.x = T)

# OMIM
if (!is.null(omim_path)){
  omim = read.delim(omim_path, header = F, stringsAsFactors = F, comment.char = "#", quote = "", check.names=F)
  colnames(omim) = c("Chromosome", "Genomic_Position_Start", "Genomic Position End", "Cyto_Location", "Computed_Cyto_Location", "MIM_Number",
    "Gene_Symbols", "Gene_Name",	"Approved_Gene_Symbol", "Entrez_Gene_ID", "Ensembl_Gene_ID", "Comments", "Phenotypes", "Mouse_Gene_Symbol-ID")
  vep = merge(vep, omim, by.x = "Mitomap_genomeloci", by.y = "Approved_Gene_Symbol", all.x = T)
}


# Gene-Panel 
if (!is.null(panels_path)){
  gene_panel = read.delim(panels_path, header = TRUE, stringsAsFactors = F, quote = "")
  vep = merge(vep, gene_panel, by.x = "Mitomap_genomeloci", by.y = "gene", all.x = T)
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


df_out  = data.frame(row.names = 1:nrow(vep), stringsAsFactors = F)
  

# Function to extract REF and ALT without underscores
extract_alleles <- function(variant) {
  # Extract the portion with REF/ALT using regex
  matches <- regmatches(variant, regexec("_[A-Za-z-]+/[A-Za-z-]+$", variant))
  if (length(matches[[1]]) > 0) {
    # Remove the leading underscore and split by "/"
    ref_alt <- unlist(strsplit(sub("_", "", matches[[1]]), "/"))
    return(ref_alt)
  }
  return(c(NA, NA)) # Return NA if no match
}

alt_ref_alleles <- t(sapply(vep$`#Uploaded_variation`, extract_alleles))



#==================================#
# Basic information of the variant #
#==================================#
print("Basic information of the variant")

df_out$CHROM = unlist(lapply(vep$Location, function(x) strsplit(x, ":")[[1]][1]))
df_out$POS = as.numeric(unlist(lapply(vep$`#Uploaded_variation`, function(x) rev(strsplit(x, "_")[[1]])[2])))
df_out$REF = alt_ref_alleles[, 1]
df_out$ALT = alt_ref_alleles[, 2]
# YBQ: en vez de coger el ref y el alt del USED_REF y el Allele que nos da VEP, parseamos el #Uploaded_variation, porque VEP no funciona bien en algunos genes del chrM.  
#df_out$REF = vep$USED_REF
#df_out$ALT = vep$Allele

df_out$REF_COUNT = vep$SAMPLE_AD_REF
df_out$ALT_COUNT = vep$SAMPLE_AD_ALT
df_out$AF = vep$SAMPLE_AF
df_out <- df_out %>% dplyr::mutate(GT = ifelse(as.numeric(AF) >= 0.95, "1/1", "0/1"))

df_out$Location = vep$Location
df_out$SYMBOL = vep$SYMBOL
df_out$GeneLoci = vep$Mitomap_genomeloci
df_out$Gene_full_name = vep$Gene_full_name
df_out$FILTER = vep$SAMPLE_filter
if (!is.null(glowgenes_path)) df_out$GLOWgenes = vep$GLOWgenes
if ((!is.null(genefilter_path)) & (!is.null(glowgenes_path))) df_out$GLOWgenes[df_out$SYMBOL %in% genefilter$V1] = 0 # We assume the genes of the list are the genes from the panel
df_out$VARIANT_CLASS = vep$VARIANT_CLASS
#df_out$Panels_name = vep$panels


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
df_out$Functional_effect_general = vep$MitImpact_Functional_effect_general
df_out$Functional_effect_detailed = vep$MitImpact_Functional_effect_detailed
#df_out$INTRON = vep$INTRON
#df_out$EXON = vep$EXON
df_out$HGVSc = vep$HGVSc
df_out$HGVSp = vep$HGVSp
df_out$HGVS = vep$MitImpact_HGVS 


df_out$MitoMap_aachange = vep$Mitomap_disease_aachange
df_out$MitoMap_disease = vep$Mitomap_disease_Disease
df_out$MitoMap_disease_status = vep$Mitomap_disease_DiseaseStatus
df_out$MitoMap_disease_HGFL = vep$Mitomap_disease_HGFL
df_out$MitoTip_score = vep$MitoTip_MitoTIP_Score
df_out$MitoTip_quartile = vep$MitoTip_MitoTIP_Quartile
df_out$MitoTip_count = vep$MitoTip_MitoTIP_Count
df_out$MitoTip_percentage = vep$MitoTip_MitoTIP_Percentage
df_out$MitoTip_status = vep$MitoTip_MitoTIP_Status


#df_out$DISTANCE = as.numeric(vep$DISTANCE)
#df_out$STRAND = vep$STRAND
df_out$Interpro_domain = vep$Interpro_domain
df_out$Interpro_domain = vep$Interpro_domain
#df_out$Domino_Score = vep$Domino_Score


#===============#
# Pathogenicity #
#===============#
print("Pathogenicity")

df_out$CLNSIG = vep$ClinVar_CLNSIG
df_out$CLNREVSTAT = vep$ClinVar_CLNREVSTAT
df_out$CLNDN = vep$ClinVar_CLNDN
df_out$CLNSIGCONF = vep$ClinVar_CLNSIGCONF                                                                      
df_out$OMIM_phenotype = vep$Phenotypes
df_out$Orphanet_disorder = vep$Orphanet_disorder
df_out$Orphanet_association_type = vep$Orphanet_association_type
df_out$HPO_name = vep$HPO_name
df_out$PUBMED = vep$PUBMED


#=============#
# Frequencies #
#=============#
print("Frequencies")

df_out$GenBank_AC = vep$Mitomap_GenBank_AC
df_out$GenBank_AF = vep$Mitomap_GenBank_AF

df_out$gnomAD_AC_hom = vep$gnomAD_AC_hom
df_out$gnomAD_AC_het = vep$gnomAD_AC_het
df_out$gnomAD_AF_hom = vep$gnomAD_AF_hom
df_out$gnomAD_AF_het = vep$gnomAD_AF_het
df_out$gnomAD_AN = vep$gnomAD_AN
df_out$gnomAD_max_observed_heteroplasmy = vep$gnomAD_max_observed_heteroplasmy
df_out$gnomAD_hap_defining_variant = vep$gnomAD_full_hap_defining_variant
df_out$gnomAD_hap_defining_variant[df_out$gnomAD_hap_defining_variant == 1] <- TRUE
df_out$gnomAD_filter = vep$gnomAD_FILTER

df_out$kaviar_AC = vep$kaviar_AC
df_out$kaviar_AF = vep$kaviar_AF

df_out$ToMMo_54KJPN_AC = vep$MitImpact_ToMMo_54KJPN_AC
df_out$ToMMo_54KJPN_AF = vep$MitImpact_ToMMo_54KJPN_AF
df_out$ToMMo_54KJPN_AN = vep$MitImpact_ToMMo_54KJPN_AN

df_out$HelixMTdb_AC_hom = vep$MitImpact_HelixMTdb_AC_hom
df_out$HelixMTdb_AF_hom = vep$MitImpact_HelixMTdb_AF_hom
df_out$HelixMTdb_AC_het = vep$MitImpact_HelixMTdb_AC_het
df_out$HelixMTdb_AF_het = vep$MitImpact_HelixMTdb_AF_het
df_out$HelixMTdb_mean_ARF = vep$MitImpact_HelixMTdb_mean_ARF
df_out$HelixMTdb_max_ARF = vep$MitImpact_HelixMTdb_max_ARF

#df_out$CSVS_AF = as.numeric(unlist(lapply(vep$CSVS_AF, function(x) strsplit(x, ",")[[1]][1])))
#df_out$CSVS_AC = as.numeric(unlist(lapply(vep$CSVS_AC, function(x) strsplit(x, ",")[[1]][1])))
#df_out$FJD_MAF_AF = as.numeric(unlist(lapply(vep$FJD_MAF_AF, function(x) strsplit(x, ",")[[1]][1])))
#df_out$FJD_MAF_AC = as.numeric(unlist(lapply(vep$FJD_MAF_AC, function(x) strsplit(x, ",")[[1]][1])))
##add new columns del MAF_FJD de DHR vs pseudocontroles (SON DE OJO LOS PSEUDOCONTROLES)
#df_out$FJD_MAF_AF_DS_IRD = as.numeric(unlist(lapply(vep$FJD_MAF_AF_DS_irdt, function(x) strsplit(x, ",")[[1]][1])))
#df_out$FJD_MAF_AC_DS_IRD = as.numeric(unlist(lapply(vep$FJD_MAF_AC_DS_irdt, function(x) strsplit(x, ",")[[1]][1])))
#df_out$FJD_MAF_AF_P_IRD = as.numeric(unlist(lapply(vep$FJD_MAF_AF_P_eyeg, function(x) strsplit(x, ",")[[1]][1])))
#df_out$FJD_MAF_AC_P_IRD = as.numeric(unlist(lapply(vep$FJD_MAF_AC_P_eyeg, function(x) strsplit(x, ",")[[1]][1])))
                                             
#df_out$denovoVariants_SAMPLE_CT = vep$denovoVariants_SAMPLE_CT



#==========================#
# Pathogenicity prediction #
#==========================#
print("Pathogenicity prediction")

df_out$CADD_PHRED = as.numeric(vep$MitImpact_CADD_phred_score)
df_out$CADD_RAW = as.numeric(vep$MitImpact_CADD_score)
df_out$MutScore = as.numeric(vep$Mut_Score)
df_out$REVELScore = as.numeric(vep$REVEL_Score)
df_out$tRNA_APOGEE_Score = as.numeric(vep$`tRNA_APOGEE_t-APOGEE_unbiased_score`)

patho_norm_func = function(predictions){
  predictions = gsub(";", ",", predictions)
  predictions = tolower(predictions)
  unlist(lapply(predictions, function(x) {
    y = strsplit(x,",")[[1]]
    y = y[!y %in% c("-", ".", "U", "")]
    y[y %in% c("tolerated", "tolerated_low_confidence", "benign", "n", "l", "p", "t","neutral","likely-benign",
              "likely_benign")] = "T"
    y[y %in% c("deleterious", "deleterious_low_confidence", "probably_damaging","disease", 
               "possibly_damaging", "a", "m", "h", "d", "Dominant", "Recessive","pathogenic")] = "D"
    if ("D" %in% y) { return("D") }
    else if ("T" %in% y) { return("T") }
    else {return("")}
  }))
}  

df_pathogenic_predictors = data.frame(row.names = 1:nrow(vep))
#df_pathogenic_predictors$SIFT = patho_norm_func(vep$SIFT)
#df_pathogenic_predictors$PolyPhen = patho_norm_func(vep$PolyPhen)
df_pathogenic_predictors$Polyphen2_HDIV_pred = patho_norm_func(vep$Polyphen2_HDIV_pred)
df_pathogenic_predictors$Polyphen2_HVAR_pred = patho_norm_func(vep$Polyphen2_HVAR_pred)
df_pathogenic_predictors$LRT_pred = patho_norm_func(vep$LRT_pred)
df_pathogenic_predictors$`M-CAP_pred` = patho_norm_func(vep$`M-CAP_pred`)
df_pathogenic_predictors$MetaLR_pred = patho_norm_func(vep$MetaLR_pred)
df_pathogenic_predictors$MetaSVM_pred = patho_norm_func(vep$MetaSVM_pred)
df_pathogenic_predictors$MutationAssessor_pred = patho_norm_func(vep$MutationAssessor_pred)
df_pathogenic_predictors$MutationTaster_pred = patho_norm_func(vep$MutationTaster_pred)
df_pathogenic_predictors$PROVEAN_pred = patho_norm_func(vep$PROVEAN_pred)
#df_pathogenic_predictors$FATHMM_pred = patho_norm_func(vep$FATHMM_pred)
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
df_pathogenic_predictors$`fathmm-XF_coding_pred` = patho_norm_func(vep$`fathmm-XF_coding_pred`)

# los predictores sacadois de MitImpact
df_pathogenic_predictors$PolyPhen = patho_norm_func(vep$MitImpact_PolyPhen2)
df_pathogenic_predictors$SIFT = patho_norm_func(vep$MitImpact_SIFT)
df_pathogenic_predictors$SIFT4G = patho_norm_func(vep$MitImpact_SIFT4G)
df_pathogenic_predictors$VEST = patho_norm_func(vep$MitImpact_VEST)
df_pathogenic_predictors$MITOCLASS = patho_norm_func(vep$MitImpact_Mitoclass1)
df_pathogenic_predictors$SNPDryad = patho_norm_func(vep$MitImpact_SNPDryad)
df_pathogenic_predictors$FATHMM = patho_norm_func(vep$MitImpact_FATHMM)
df_pathogenic_predictors$AlphaMissense = patho_norm_func(vep$MitImpact_AlphaMissense)
df_pathogenic_predictors$PROVEAN = patho_norm_func(vep$MitImpact_PROVEAN)
df_pathogenic_predictors$EFIN_SP = patho_norm_func(vep$MitImpact_EFIN_SP)
df_pathogenic_predictors$EFIN_HD = patho_norm_func(vep$MitImpact_EFIN_HD)
df_pathogenic_predictors$MLC = patho_norm_func(vep$MitImpact_MLC)
df_pathogenic_predictors$APOGEE1 = patho_norm_func(vep$MitImpact_APOGEE1)
df_pathogenic_predictors$APOGEE2 = patho_norm_func(vep$MitImpact_APOGEE2)
df_pathogenic_predictors$CAROL = patho_norm_func(vep$MitImpact_CAROL)
df_pathogenic_predictors$Condel = patho_norm_func(vep$MitImpact_Condel)
df_pathogenic_predictors$PANTHER = patho_norm_func(vep$MitImpact_PANTHER)


df_out$N_Pathogenic_pred = apply(df_pathogenic_predictors, 1, function(x) table(x)["D"])
df_out$N_Pathogenic_pred[is.na(df_out$N_Pathogenic_pred)] = 0
df_out$N_Benign_pred = apply(df_pathogenic_predictors, 1, function(x) table(x)["T"])
df_out$N_Benign_pred[is.na(df_out$N_Benign_pred)] = 0
df_out$N_predictions = df_out$N_Pathogenic_pred + df_out$N_Benign_pred
df_out$Pathogenic_pred = apply(df_pathogenic_predictors, 1, function(x) paste(names(x)[which(x == "D")], collapse = ","))
df_out$Benign_pred = apply(df_pathogenic_predictors, 1, function(x) paste(names(x)[which(x == "T")], collapse = ","))



#============================#
# Conservation and phylogeny #
#============================#
print("Conservation and phylogeny")

df_out$LoFtool = as.numeric(vep$LoFtool)
#antiguo "ExACpLI ahora es pLI_gene_value
df_out$ExACpLI = as.numeric(vep$pLI_gene_value)
#df_out$gnomAD_exomes_CCR = vep$gnomAD_exomes_CCR
df_out$phastCons30way_mammalian = as.numeric(vep$phastCons470way_mammalian)
df_out$phyloP30way_mammalian = as.numeric(vep$phyloP470way_mammalian)
df_out$MGI_mouse_phenotype = vep$MGI_mouse_phenotype_filt





#=====================================================#
#Expression, process, route, function and interaction #
#=====================================================#
print("Expression, process, route, function and interaction")

df_out$GTEx_V8_gene = vep$GTEx_V8_eQTL_gene
df_out$GTEx_V8_tissue = vep$GTEx_V8_eQTL_tissue
df_out$`Expression_GNF-Atlas` = vep$Expression.GNF.Atlas.
df_out$Pathway_KEGG = vep$Pathway.KEGG._full
df_out$GO_biological_process = vep$GO_biological_process
df_out$GO_cellular_component = vep$GO_cellular_component
df_out$GO_molecular_function = vep$GO_molecular_function
df_out$Interactions_IntAct = vep$Interactions.IntAct.

df_out$retina_RNA_tissue_consensus = round(vep$retina,2)
df_out$testis_RNA_tissue_consensus = round(vep$testis,2)
df_out$kidney_RNA_tissue_consensus = round(vep$kidney,2)
df_out$brain_max_RNA_tissue_consensus = round(vep$brain_max,2)
df_out$glands_max_RNA_tissue_consensus = round(vep$glands_max,2)
df_out$digestive_max_RNA_tissue_consensus = round(vep$digestive_max,2)
df_out$heart_RNA_tissue_consensus = round(vep$heart.muscle,2)
df_out$liver_RNA_tissue_consensus = round(vep$liver,2)
df_out$lung_RNA_tissue_consensus = round(vep$lung,2)
df_out$pancreas_RNA_tissue_consensus = round(vep$pancreas,2)
df_out$skel_muscle_RNA_tissue_consensus = round(vep$skeletal.muscle,2)
df_out$skin_RNA_tissue_consensus = round(vep$skin,2)
df_out$mean_expression_RNA_tissue_consensus = round(vep$mean_exp,2)
df_out$retina_ratio_exp_RNA_tissue_consensus = round(vep$retina_ratio, 2)
                           
df_out$Original_pos = vep$SAMPLE_Original_pos
df_out$variant_id = vep$SAMPLE_variant_id





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
