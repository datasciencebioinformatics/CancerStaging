#####################################################################################################################
# This script will take the TSV file (metadata), unstranded.rna_seq.augmented_star_gene_counts (rna-seq count data), 
# the name of genes and samples already processed for the primary_diagnosis=Squamous cell carcinoma, NOS
#####################################################################################################################
# Reading the contents of TSV file using read_tsv() method
merged_data_patient_info_file<- "/home/felipe/googledrive/Cancer_staging/merged_data_patient_info.tsv"

# Load metadata table
merged_data_patient_info     <-read.table(file = merged_data_patient_info_file, sep = '\t', header = TRUE,fill=TRUE)   
#####################################################################################################################
