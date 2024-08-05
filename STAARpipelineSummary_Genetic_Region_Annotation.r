###########################################################
# Annotate rare variants in genetic regions
# using STAARpipelineSummary
# Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 03/10/2022
###########################################################
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

## output path
output_path <- "/GWAS/data/nov_gwas_staar_update/plt/snv/snv_5mask/summary/step6_sliding_window_maf5percent_1/"

## QC_label
QC_label <- "annotation/filter"
## geno_missing_imputation
geno_missing_imputation <- "mean"
## variant_type
variant_type <- "SNV"
# method_cond
method_cond <- "optimal"

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load("/GWAS/data/nov_gwas_staar_update/step0_preSTAAR/Annotation_name_catalog.Rdata"))
# Or equivalently
# Annotation_name_catalog <- read.csv("/path_to_the_file/Annotation_name_catalog.csv")
## Annotation name
Annotation_name <- c("rs_num","GENCODE.Category","GENCODE.Info","GENCODE.EXONIC.Category","MetaSVM","GeneHancer","CAGE","DHS",
                     "CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")

## Chr
chr_seq <- c(8,13,20,20)
## Start location
start_loc_seq <- c(34543807,83284912,38380454,38382454)
## End location
end_loc_seq <- c(34547806,83288911,38384453,38386453)

###########################################################
#           Main Function 
###########################################################
for(kk in 1:length(chr_seq))
{
	chr <- chr_seq[kk]
	start_loc <- start_loc_seq[kk]
	end_loc <- end_loc_seq[kk]
	
	print(paste0(chr,"_",start_loc,"_",end_loc))

	### gds file
	gds.path <- agds_dir[chr]
	genofile <- seqOpen(gds.path)

	results_info <- Sliding_Window_Info(chr=chr,genofile=genofile,obj_nullmodel=obj_nullmodel,start_loc=start_loc,end_loc=end_loc,known_loci=known_loci,
	                                    QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
	                                    Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,Annotation_name=Annotation_name)

	seqClose(genofile)

	save(results_info,file=paste0(output_path,"window_",chr,"_",start_loc,"_",end_loc,".Rdata"))
	write.csv(results_info,paste0(output_path,"window_",chr,"_",start_loc,"_",end_loc,".csv"),row.names=FALSE)
}

