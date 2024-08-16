#####################################################################################################################
# This script will take the TSV file (metadata), unstranded.rna_seq.augmented_star_gene_counts (rna-seq count data), 
#####################################################################################################################
# Reading the contents of TSV file using read_tsv() method
merged_data_patient_info_file<- "/home/felipe/googledrive/Cancer_staging/merged_data_patient_info.tsv"

# Load metadata table
merged_data_patient_info     <-read.table(file = merged_data_patient_info_file, sep = '\t', header = TRUE,fill=TRUE)   
#####################################################################################################################
# A script to load cancer data base in R
unstranded_file       <- "/home/felipe/googledrive/Cancer_staging/tables/unstranded.rna_seq.augmented_star_gene_counts.tsv"
colnames_file         <- "/home/felipe/googledrive/Cancer_staging/tables/header.txt"
rownames_file         <- "/home/felipe/googledrive/Cancer_staging/tables/gene_ids.txt"

# Load data
unstranded_data<-read.table(file = unstranded_file, sep = '\t', header = FALSE,fill=TRUE)    
colnames_data<-read.table(file = colnames_file, sep = '\t', header = FALSE,fill=TRUE)                                    
rownames_data<-read.table(file = rownames_file, sep = '\t', header = FALSE,fill=TRUE)     

# Set colnames and rownames
rownames(unstranded_data)<-rownames_data[,1]
colnames(unstranded_data)<-colnames_data[,1]
#####################################################################################################################
# A list to store the datasets
reads_count_per_project<-list()

# for each project 
for (project in rownames(table_cases_per_stage))
{
    # Set project ids
    project_ids <-merged_data_patient_info[merged_data_patient_info$project_id==project,"sample_id"]  

    # Set project data
    project_data<-unstranded_data[,which(colnames(unstranded_data) %in% project_ids)]

    # Store dataset
    reads_count_per_project[[project]]<-project_data    
}
