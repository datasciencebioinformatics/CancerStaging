# Set name of the TCGA project
TCGA_project <- "TCGA-LUSC"

# Execute TCGA coomand
t1 <- try(system(paste("mkdir /home/felipe/Documents/Cancer_staging/outputfolder/"), intern = TRUE))

# Set outfolder
output_dir=paste("/home/felipe/Documents/Cancer_staging/outputfolder/",TCGA_project,"/",sep="")

# Execute TCGA coomand
t1 <- try(system(paste("mkdir",output_dir), intern = TRUE))

# FDR threshold for the comparisson all tumor samples / all control samples
threshold_FDR<-0.05 

# log2foldchange threshold for the comparisson all tumor samples / all control samples
threshold_tumor<-1.0 

# rpkm=0,fpkm=0,tmm=0, tpm=0
list_threshold_filters<-list(raw=0,rpkm=10,fpkm=4,tmm=4, tpm=10)

# log2foldchange threshold for the comparisson all tumor samples / all control samples
threshold_stage<-1.0 

print("\nCancerStaging_SetupAllParamters.R")
