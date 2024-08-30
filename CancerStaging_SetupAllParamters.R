# Set outfolder
output_dir="/home/felipe/Documents/Cancer_staging/outputfolder/"

# FDR threshold for the comparisson all tumor samples / all control samples
threshold_FDR<-0.05 

# log2foldchange threshold for the comparisson all tumor samples / all control samples
threshold_tumor<-0.0 

# rpkm=0,fpkm=0,tmm=0, tpm=0
list_threshold_filters<-list(rpkm=4,fpkm=4,tmm=4, tpm=4)

# log2foldchange threshold for the comparisson all tumor samples / all control samples
threshold_stage<-1.0 
