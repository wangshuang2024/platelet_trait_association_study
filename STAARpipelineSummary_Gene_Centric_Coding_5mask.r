##########################################################
# Summarization and visualization of gene-centric 
# coding analysis results using STAARpipelineSummary
# Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 03/10/2022
##########################################################
rm(list=ls())
gc()

## load required packages
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)
library(STAARpipelineSummary)

###########################################################
#           User Input
###########################################################
## aGDS directory
agds_dir <- get(load("/GWAS/data/nov_gwas_staar_update/step0_preSTAAR/agds_dir.Rdata"))
## Known loci
known_loci <- get(load("/GWAS/data/nov_gwas_staar_update/plt/plt_known_loci_pruning/plt_known_loci_genome_LD_pruning.Rdata"))
## Null model
obj_nullmodel <- get(load("/GWAS/data/nov_gwas_staar_update/plt/snv/step1_fitnullmodel/plt_snv_nullmodel_GENESIS.Rdata"))

## results path
input_path <- "/GWAS/data/nov_gwas_staar_update/plt/snv/snv_5mask/step3_gene_centric_coding_maf5percent/"
output_path <- "/GWAS/data/nov_gwas_staar_update/plt/snv/snv_5mask/summary/step3_gene_centric_coding_maf5percent/"
## number of jobs
gene_centric_coding_jobs_num <- 381
## results name
gene_centric_results_name <- "plt_snv_coding"

## QC_label
QC_label <- "annotation/filter"
## geno_missing_imputation
geno_missing_imputation <- "mean"
## variant_type
variant_type <- "SNV"
## method_cond
method_cond <- "optimal"
## alpha level
alpha <- 1E-5

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

cMAC_cutoff <- 5

###########################################################
#           Main Function 
###########################################################
Gene_Centric_Coding_Results_Summary(agds_dir=agds_dir,gene_centric_coding_jobs_num=gene_centric_coding_jobs_num,
                                    input_path=input_path,output_path=output_path,gene_centric_results_name=gene_centric_results_name,
                                    obj_nullmodel=obj_nullmodel,known_loci=known_loci,cMAC_cutoff=cMAC_cutoff,
                                    method_cond=method_cond,rare_maf_cutoff=0.05,
                                    QC_label=QC_label,geno_missing_imputation=geno_missing_imputation,variant_type=variant_type,
                                    Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                    Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                    alpha=alpha,manhattan_plot=TRUE,QQ_plot=TRUE)

