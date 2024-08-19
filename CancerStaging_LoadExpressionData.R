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

merged_data_patient_info_BRCA<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-BRCA",]
merged_data_patient_info_LIHC<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-LIHC",]
merged_data_patient_info_LUAD<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-LUAD",]
merged_data_patient_info_LUSC<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-LUSC",]
merged_data_patient_info_PRAD<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-PRAD",]
merged_data_patient_info_READ<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-READ",]
merged_data_patient_info_SKCM<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-SKCM",]
merged_data_patient_info_STAD<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-STAD",]

merged_data_patient_info_BRCA<-unique(merged_data_patient_info_BRCA[which(merged_data_patient_info_BRCA$sample_id %in%  colnames(reads_count_per_project[["TCGA-BRCA"]])),c("project_id","stages", "Sample.ID")])   
merged_data_patient_info_LIHC<-unique(merged_data_patient_info_LIHC[which(merged_data_patient_info_LIHC$sample_id %in%  colnames(reads_count_per_project[["TCGA-LIHC"]])),c("project_id","stages", "Sample.ID")])    
merged_data_patient_info_LUAD<-unique(merged_data_patient_info_LUAD[which(merged_data_patient_info_LUAD$sample_id %in%  colnames(reads_count_per_project[["TCGA-LUAD"]])),c("project_id","stages", "Sample.ID")])    
merged_data_patient_info_LUSC<-unique(merged_data_patient_info_LUSC[which(merged_data_patient_info_LUSC$sample_id %in%  colnames(reads_count_per_project[["TCGA-LUSC"]])),c("project_id","stages", "Sample.ID")])    
merged_data_patient_info_PRAD<-unique(merged_data_patient_info_PRAD[which(merged_data_patient_info_PRAD$sample_id %in%  colnames(reads_count_per_project[["TCGA-PRAD"]])),c("project_id","stages", "Sample.ID")])    
merged_data_patient_info_READ<-unique(merged_data_patient_info_READ[which(merged_data_patient_info_READ$sample_id %in%  colnames(reads_count_per_project[["TCGA-READ"]])),c("project_id","stages", "Sample.ID")])    
merged_data_patient_info_SKCM<-unique(merged_data_patient_info_SKCM[which(merged_data_patient_info_SKCM$sample_id %in%  colnames(reads_count_per_project[["TCGA-SKCM"]])),c("project_id","stages", "Sample.ID")])    

# Save table
merged_data_patient_info_merged<-rbind(merged_data_patient_info_BRCA,merged_data_patient_info_LIHC,merged_data_patient_info_LUAD,merged_data_patient_info_LUSC,merged_data_patient_info_PRAD,merged_data_patient_info_READ,merged_data_patient_info_SKCM)


# Check the number of samples per stage
# Cases per stage
table_cases_per_stage<-table(merged_data_patient_info_merged$project_id, merged_data_patient_info_merged$stages)

# Cases per stage
table_cases_per_stage<-table(merged_data_patient_info_merged$project_id, merged_data_patient_info_merged$stages)

# Cases per stage
table_cases_per_stage<-table_cases_per_stage[,c("Stage I","Stage II","Stage III","Stage IV")]
#####################################################################################################################
# Organize how to send to Carles
write_tsv(merged_data_patient_info, "/home/felipe/googledrive/Cancer_staging/merged_data_patient_info.tsv")
