#####################################################################################################################
# This script will take the TSV file (metadata), unstranded.rna_seq.augmented_star_gene_counts (rna-seq count data), 
#####################################################################################################################
# Reading the contents of TSV file using read_tsv() method
merged_data_patient_info_file<- "/home/felipe/Documents/Cancer_staging/merged_data_patient_info.tsv"

# Load metadata table
merged_data_patient_info     <-read.table(file = merged_data_patient_info_file, sep = '\t', header = TRUE,fill=TRUE)   
#####################################################################################################################
# A script to load cancer data base in R
unstranded_file       <- "/home/felipe/Documents/Cancer_staging/tables/unstranded.rna_seq.augmented_star_gene_counts.tsv"
colnames_file         <- "/home/felipe/Documents/Cancer_staging/tables/header.txt"
rownames_file         <- "/home/felipe/Documents/Cancer_staging/tables/gene_ids.txt"

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
#####################################################################################################################
# Count the number of reads per project
dim(reads_count_per_project[["TCGA-BRCA"]])
dim(reads_count_per_project[["TCGA-LIHC"]])
dim(reads_count_per_project[["TCGA-LUAD"]])
dim(reads_count_per_project[["TCGA-LUSC"]])
dim(reads_count_per_project[["TCGA-PRAD"]])
dim(reads_count_per_project[["TCGA-READ"]])
dim(reads_count_per_project[["TCGA-SKCM"]])
dim(reads_count_per_project[["TCGA-STAD"]])

# Count total number of samples
dim(reads_count_per_project[["TCGA-BRCA"]])[2]+dim(reads_count_per_project[["TCGA-LIHC"]])[2]+dim(reads_count_per_project[["TCGA-LUAD"]])[2]+dim(reads_count_per_project[["TCGA-LUSC"]])[2]+dim(reads_count_per_project[["TCGA-PRAD"]])[2]+dim(reads_count_per_project[["TCGA-READ"]])[2]+dim(reads_count_per_project[["TCGA-SKCM"]])[2]+dim(reads_count_per_project[["TCGA-STAD"]])[2]

merged_data_patient_info_BRCA<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-BRCA",]
merged_data_patient_info_LIHC<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-LIHC",]
merged_data_patient_info_LUAD<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-LUAD",]
merged_data_patient_info_LUSC<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-LUSC",]
merged_data_patient_info_PRAD<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-PRAD",]
merged_data_patient_info_READ<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-READ",]
merged_data_patient_info_SKCM<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-SKCM",]
merged_data_patient_info_STAD<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-STAD",]

all_samples<-c(colnames(reads_count_per_project[["TCGA-BRCA"]]),
colnames(reads_count_per_project[["TCGA-LIHC"]]),
colnames(reads_count_per_project[["TCGA-LUAD"]]),
colnames(reads_count_per_project[["TCGA-LUSC"]]),
colnames(reads_count_per_project[["TCGA-PRAD"]]),
colnames(reads_count_per_project[["TCGA-READ"]]),
colnames(reads_count_per_project[["TCGA-SKCM"]]),
colnames(reads_count_per_project[["TCGA-STAD"]]))

# I am checckin here, number of samples - 19-08-2024
# Save table
merged_data_patient_info_merged<-rbind(merged_data_patient_info_BRCA,merged_data_patient_info_LIHC,merged_data_patient_info_LUAD,merged_data_patient_info_LUSC,merged_data_patient_info_PRAD,merged_data_patient_info_READ,merged_data_patient_info_SKCM, merged_data_patient_info_STAD)

# Use only data that has read counts
merged_data_patient_info_count<-unique(merged_data_patient_info_merged[which( merged_data_patient_info_merged$sample_id %in%  all_samples),c("project_id","stages", "Sample.ID")])   

# Cases per stage
table_cases_per_stage<-table(merged_data_patient_info_count$project_id)

#####################################################################################################################
# Organize how to send to Carles
write_tsv(merged_data_patient_info_merged, "/home/felipe/Documents/Cancer_staging/merged_data_patient_info.tsv")

# Total number of samples
#TCGA-BRCA TCGA-LIHC TCGA-LUAD TCGA-LUSC TCGA-PRAD TCGA-READ TCGA-SKCM TCGA-STAD 
#     1180       411       579       536       545       167       104       440 

# Count the number of reads per project
reads_count_all_projects<-cbind(reads_count_per_project[["TCGA-BRCA"]],
reads_count_per_project[["TCGA-LIHC"]],
reads_count_per_project[["TCGA-LUAD"]],
reads_count_per_project[["TCGA-LUSC"]],
reads_count_per_project[["TCGA-PRAD"]],
reads_count_per_project[["TCGA-READ"]],
reads_count_per_project[["TCGA-SKCM"]],
reads_count_per_project[["TCGA-STAD"]])
#####################################################################################################################
write_tsv(reads_count_all_projects, "/home/felipe/Documents/Cancer_staging/reads_count_all_projects.tsv")
#####################################################################################################################
