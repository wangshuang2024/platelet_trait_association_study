##################################################################################
# Gene-centric analysis for coding rare variants in long masks using STAARpipeline
# Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 03/10/2022
##################################################################################
rm(list=ls())
gc()

## load required packages
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)

###########################################################
#           User Input
###########################################################
## aGDS directory
agds_dir <- get(load("/GWAS/data/nov_gwas_staar_update/step0_preSTAAR/agds_dir.Rdata"))
## Null model
obj_nullmodel <- get(load("/GWAS/data/nov_gwas_staar_update/plt/snv/step1_fitnullmodel/plt_snv_nullmodel_GENESIS.Rdata"))

## QC_label
QC_label <- "annotation/filter"
## variant_type
variant_type <- "SNV"
## geno_missing_imputation
geno_missing_imputation <- "mean"

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load("/GWAS/data/nov_gwas_staar_update/step0_preSTAAR/Annotation_name_catalog.Rdata"))
# Or equivalently
# Annotation_name_catalog <- read.csv("/path_to_the_file/Annotation_name_catalog.csv")
## Use_annotation_weights
Use_annotation_weights <- TRUE
## Annotation name
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")

## output path
output_path <- "/GWAS/data/nov_gwas_staar_update/plt/snv/snv_5mask/step3_gene_centric_coding_maf5percent/"
## output file name
output_file_name <- "plt_snv_coding"
## input array id from batch file (Harvard FAS RC cluster)
arrayid_longmask <- as.numeric(commandArgs(TRUE)[1])

###########################################################
#           Main Function 
###########################################################
## gene number in job
gene_num_in_array <- 50 
group.num.allchr <- ceiling(table(genes_info[,2])/gene_num_in_array)
sum(group.num.allchr)

### exclude large genes
arrayid <- c(57,112,112,113,113,113,113,113,113,113)
sub_seq_id <- c(840,543,544,575,576,577,578,579,580,582)

region_spec <- data.frame(arrayid,sub_seq_id) 
sub_seq_id <- ((arrayid_longmask-1)*5+1):min(arrayid_longmask*5,length(arrayid))

### gds file
genes <- genes_info

results_coding <- c()
for(kk in sub_seq_id)
{
	print(kk)
	arrayid <- region_spec$arrayid[kk]
	sub_id <- region_spec$sub_seq_id[kk]
	
	chr <- which.max(arrayid <= cumsum(group.num.allchr))
	gds.path <- agds_dir[chr]
	genofile <- seqOpen(gds.path)
	
	genes_info_chr <- genes_info[genes_info[,2]==chr,]
	gene_name <- genes_info_chr[sub_id,1]

	results <- Gene_Centric_Coding(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
	                               rare_maf_cutoff=0.05,rv_num_cutoff=2,
	                               QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
	                               Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
	                               Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
	
	results_coding <- append(results_coding,results)

	seqClose(genofile)
}

save(results_coding,file=paste0(output_path,output_file_name,"_",arrayid_longmask+379,".Rdata"))

