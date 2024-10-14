####################################################################################################
# Function to calculate the tmp
# A dataframe with the read count per gene for the patient_j	
# A vector with gene length in the same order of read_counts_per_patient_j
tpm <- function(counts, lengths) 
{
	rate <- counts / lengths
	rate / sum(rate) * 1e6  
}
####################################################################################################
# Set name of the TCGA project
TCGA_project <- "TCGA-STAD"

# Execute TCGA coomand
t1 <- try(system(paste("mkdir /home/felipe/Documents/Cancer_staging/outputfolder/"), intern = TRUE))

# Set outfolder
output_dir=paste("/home/felipe/Documents/Cancer_staging/outputfolder/",TCGA_project,"/",sep="")

# Execute TCGA coomand
t1 <- try(system(paste("mkdir",output_dir), intern = TRUE))

# FDR threshold for the comparisson all tumor samples / all control samples
threshold_FDR<-0.05 

# log2foldchange threshold for the comparisson all tumor samples / all control samples
threshold_tumor<-0.0 

# rpkm=0,fpkm=0,tmm=0, tpm=0
list_threshold_filters<-list(raw=4,rpkm=4,fpkm=4,tmm=4, tpm=4,tpm_calc=4 )

# log2foldchange threshold for the comparisson all tumor samples / all control samples
threshold_stage<-1.0 

# File to save results
results_files <- paste(output_dir,"/outputfile_FDR_005_TUMOR_00_FILTER_0_STAGE_10.txt",sep="")

print("\nCancerStaging_SetupAllParamters.R")
####################################################################################################
# Write to to file the number of cases
cat(paste("results_files","       : ", results_files, sep=""),file=results_files,sep="\n", append=FALSE)
cat(paste("threshold_FDR","       : ", threshold_FDR, sep=""),file=results_files,sep="\n", append=TRUE)
cat(paste("threshold_tumor","     : ", threshold_tumor, sep=""),file=results_files,sep="\n", append=TRUE)
cat(paste("threshold_stage","     : ", threshold_stage, sep=""),file=results_files,sep="\n", append=TRUE)
cat(paste("threshold_filters raw  : ", list_threshold_filters[["raw"]], sep=""),file=results_files,sep="\n", append=TRUE)
cat(paste("threshold_filters rpkm : ", list_threshold_filters[["rpkm"]], sep=""),file=results_files,sep="\n", append=TRUE)
cat(paste("threshold_filters fpkm : ", list_threshold_filters[["fpkm"]], sep=""),file=results_files,sep="\n", append=TRUE)
cat(paste("threshold_filters tmm : ", list_threshold_filters[["tmm"]], sep=""),file=results_files,sep="\n", append=TRUE)
cat(paste("threshold_filters tpm : ", list_threshold_filters[["tpm"]], sep=""),file=results_files,sep="\n", append=TRUE)
####################################################################################################
