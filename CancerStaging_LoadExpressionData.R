#####################################################################################################################
# This script will take the TSV file (metadata), unstranded.rna_seq.augmented_star_gene_counts (rna-seq count data), 
#####################################################################################################################
# Reading the contents of TSV file using read_tsv() method
merged_data_patient_info_file<- "/home/felipe/Documents/Cancer_staging/merged_data_patient_info.tsv"

# Load metadata table
merged_data_patient_info     <-read.table(file = merged_data_patient_info_file, sep = '\t', header = TRUE,fill=TRUE)   
#####################################################################################################################
# A script to load cancer data base in R
unstranded_raw_file      <- "/home/felipe/Documents/Cancer_staging/tables/unstranded.rna_seq.augmented_star_gene_counts.tsv"
unstranded_fpkm_file     <- "/home/felipe/Documents/Cancer_staging/tables/fpkm_unstranded.rna_seq.augmented_star_gene_counts.tsv"
unstranded_tpm_file      <- "/home/felipe/Documents/Cancer_staging/tables/tpm_unstranded.rna_seq.augmented_star_gene_counts.tsv"

colnames_file         <- "/home/felipe/Documents/Cancer_staging/tables/header.txt"
rownames_file         <- "/home/felipe/Documents/Cancer_staging/tables/gene_ids.txt"

# Load data
unstranded_raw_data <-read.table(file = unstranded_raw_file, sep = '\t', header = FALSE,fill=TRUE)    
unstranded_fpkm_data<-read.table(file = unstranded_fpkm_file, sep = '\t', header = FALSE,fill=TRUE)    
unstranded_tpm_data <-read.table(file = unstranded_tpm_file, sep = '\t', header = FALSE,fill=TRUE)    

colnames_data<-read.table(file = colnames_file, sep = '\t', header = FALSE,fill=TRUE)                                    
rownames_data<-read.table(file = rownames_file, sep = '\t', header = FALSE,fill=TRUE)     

# Set colnames and rownames
rownames(unstranded_raw_data)<-rownames_data[,1]
colnames(unstranded_raw_data)<-colnames_data[,1]
rownames(unstranded_tpm_data)<-rownames_data[,1]
colnames(unstranded_tpm_data)<-colnames_data[,1]
rownames(unstranded_fpkm_data)<-rownames_data[,1]
colnames(unstranded_fpkm_data)<-colnames_data[,1]
#####################################################################################################################
# A list to store the datasets
reads_count_per_project_raw  <-list()
reads_count_per_project_tpm  <-list()
reads_count_per_project_fpkm <-list()

# for each project 
for (project in rownames(table_cases_per_stage))
{
    # Set project ids
    sample_ids <-merged_data_patient_info[merged_data_patient_info$project_id==project,"sample_id"]  

    # Set project data
    project_data_raw<-unstranded_raw_data[,which(colnames(unstranded_raw_data) %in% sample_ids)]
    project_data_tpm<-unstranded_tpm_data[,which(colnames(unstranded_tpm_data) %in% sample_ids)]
    project_data_fpkm<-unstranded_fpkm_data[,which(colnames(unstranded_fpkm_data) %in% sample_ids)]

    # Store dataset
    reads_count_per_project_raw[[project]]<-project_data_raw    
    reads_count_per_project_tpm[[project]]<-project_data_tpm
    reads_count_per_project_fpkm[[project]]<-project_data_fpkm
}
#####################################################################################################################
# Count the number of reads per project
dim(reads_count_per_project_raw[["TCGA-BRCA"]])
dim(reads_count_per_project_raw[["TCGA-LIHC"]])
dim(reads_count_per_project_raw[["TCGA-LUAD"]])
dim(reads_count_per_project_raw[["TCGA-LUSC"]])
dim(reads_count_per_project_raw[["TCGA-PRAD"]])
dim(reads_count_per_project_raw[["TCGA-READ"]])
dim(reads_count_per_project_raw[["TCGA-SKCM"]])
dim(reads_count_per_project_raw[["TCGA-STAD"]])

# Count the number of reads per project
dim(reads_count_per_project_tpm[["TCGA-BRCA"]])
dim(reads_count_per_project_tpm[["TCGA-LIHC"]])
dim(reads_count_per_project_tpm[["TCGA-LUAD"]])
dim(reads_count_per_project_tpm[["TCGA-LUSC"]])
dim(reads_count_per_project_tpm[["TCGA-PRAD"]])
dim(reads_count_per_project_tpm[["TCGA-READ"]])
dim(reads_count_per_project_tpm[["TCGA-SKCM"]])
dim(reads_count_per_project_tpm[["TCGA-STAD"]])

# Count the number of reads per project
dim(reads_count_per_project_fpkm[["TCGA-BRCA"]])
dim(reads_count_per_project_fpkm[["TCGA-LIHC"]])
dim(reads_count_per_project_fpkm[["TCGA-LUAD"]])
dim(reads_count_per_project_fpkm[["TCGA-LUSC"]])
dim(reads_count_per_project_fpkm[["TCGA-PRAD"]])
dim(reads_count_per_project_fpkm[["TCGA-READ"]])
dim(reads_count_per_project_fpkm[["TCGA-SKCM"]])
dim(reads_count_per_project_fpkm[["TCGA-STAD"]])

# Count total number of samples
dim(reads_count_per_project_raw[["TCGA-BRCA"]])[2]+dim(reads_count_per_project_raw[["TCGA-LIHC"]])[2]+dim(reads_count_per_project_raw[["TCGA-LUAD"]])[2]+dim(reads_count_per_project_raw[["TCGA-LUSC"]])[2]+dim(reads_count_per_project_raw[["TCGA-PRAD"]])[2]+dim(reads_count_per_project_tpm[["TCGA-READ"]])[2]+dim(reads_count_per_project_raw[["TCGA-SKCM"]])[2]+dim(reads_count_per_project_raw[["TCGA-STAD"]])[2]
dim(reads_count_per_project_fpkm[["TCGA-BRCA"]])[2]+dim(reads_count_per_project_fpkm[["TCGA-LIHC"]])[2]+dim(reads_count_per_project_fpkm[["TCGA-LUAD"]])[2]+dim(reads_count_per_project_fpkm[["TCGA-LUSC"]])[2]+dim(reads_count_per_project_fpkm[["TCGA-PRAD"]])[2]+dim(reads_count_per_project_fpkm[["TCGA-READ"]])[2]+dim(reads_count_per_project_fpkm[["TCGA-SKCM"]])[2]+dim(reads_count_per_project_fpkm[["TCGA-STAD"]])[2]
dim(reads_count_per_project_tpm[["TCGA-BRCA"]])[2]+dim(reads_count_per_project_tpm[["TCGA-LIHC"]])[2]+dim(reads_count_per_project_tpm[["TCGA-LUAD"]])[2]+dim(reads_count_per_project_tpm[["TCGA-LUSC"]])[2]+dim(reads_count_per_project_tpm[["TCGA-PRAD"]])[2]+dim(reads_count_per_project_tpm[["TCGA-READ"]])[2]+dim(reads_count_per_project_tpm[["TCGA-SKCM"]])[2]+dim(reads_count_per_project_tpm[["TCGA-STAD"]])[2]

merged_data_patient_info_BRCA<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-BRCA",]
merged_data_patient_info_LIHC<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-LIHC",]
merged_data_patient_info_LUAD<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-LUAD",]
merged_data_patient_info_LUSC<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-LUSC",]
merged_data_patient_info_PRAD<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-PRAD",]
merged_data_patient_info_READ<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-READ",]
merged_data_patient_info_SKCM<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-SKCM",]
merged_data_patient_info_STAD<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-STAD",]

all_samples<-c(colnames(reads_count_per_project_raw[["TCGA-BRCA"]]),
colnames(reads_count_per_project_raw[["TCGA-LIHC"]]),
colnames(reads_count_per_project_raw[["TCGA-LUAD"]]),
colnames(reads_count_per_project_raw[["TCGA-LUSC"]]),
colnames(reads_count_per_project_raw[["TCGA-PRAD"]]),
colnames(reads_count_per_project_raw[["TCGA-READ"]]),
colnames(reads_count_per_project_raw[["TCGA-SKCM"]]),
colnames(reads_count_per_project_raw[["TCGA-STAD"]]))

# Use only data that has read counts
merged_data_patient_info      <-unique(merged_data_patient_info[which( merged_data_patient_info$sample_id %in%  all_samples),])   
merged_data_patient_info_count<-unique(merged_data_patient_info[which( merged_data_patient_info$sample_id %in%  all_samples),c("project_id","stages", "Sample.ID","sample_id", "Sample.Type")])   

# Cases per stage
table_cases_per_stage<-table(merged_data_patient_info_count$project_id)
#####################################################################################################################
# Organize how to send to Carles
write_tsv(merged_data_patient_info, "/home/felipe/Documents/Cancer_staging/merged_data_patient_info.tsv")
#####################################################################################################################
# Total number of samples
#TCGA-BRCA TCGA-LIHC TCGA-LUAD TCGA-LUSC TCGA-PRAD TCGA-READ TCGA-SKCM TCGA-STAD 
#     1180       411       579       536       545       167       104       440 
# Count the number of reads per project
df_reads_count_all_projects_raw<-cbind(reads_count_per_project_raw[["TCGA-BRCA"]],
reads_count_per_project_raw[["TCGA-LIHC"]],
reads_count_per_project_raw[["TCGA-LUAD"]],
reads_count_per_project_raw[["TCGA-LUSC"]],
reads_count_per_project_raw[["TCGA-PRAD"]],
reads_count_per_project_raw[["TCGA-READ"]],
reads_count_per_project_raw[["TCGA-SKCM"]],
reads_count_per_project_raw[["TCGA-STAD"]])

# Count the number of reads per project
df_reads_count_all_projects_fpkm<-cbind(reads_count_per_project_fpkm[["TCGA-BRCA"]],
reads_count_per_project_fpkm[["TCGA-LIHC"]],
reads_count_per_project_fpkm[["TCGA-LUAD"]],
reads_count_per_project_fpkm[["TCGA-LUSC"]],
reads_count_per_project_fpkm[["TCGA-PRAD"]],
reads_count_per_project_fpkm[["TCGA-READ"]],
reads_count_per_project_fpkm[["TCGA-SKCM"]],
reads_count_per_project_fpkm[["TCGA-STAD"]])

# Count the number of reads per project
df_reads_count_all_projects_tpm<-cbind(reads_count_per_project_tpm[["TCGA-BRCA"]],
reads_count_per_project_tpm[["TCGA-LIHC"]],
reads_count_per_project_tpm[["TCGA-LUAD"]],
reads_count_per_project_tpm[["TCGA-LUSC"]],
reads_count_per_project_tpm[["TCGA-PRAD"]],
reads_count_per_project_tpm[["TCGA-READ"]],
reads_count_per_project_tpm[["TCGA-SKCM"]],
reads_count_per_project_tpm[["TCGA-STAD"]])
###########################################################################################################################
#write.table(df_reads_count_all_projects_raw,  paste(output_dir,"df_reads_count_all_projects_raw.tsv",sep="/"), na = "NA", append = TRUE, col.names = TRUE, row.names = TRUE, sep = "\t", quote = TRUE)
#write.table(df_reads_count_all_projects_fpkm, paste(output_dir,"df_reads_count_all_projects_fpkm.tsv",sep="/"), na = "NA", append = TRUE, col.names = TRUE, row.names = TRUE, sep = "\t", quote = TRUE)
#write.table(df_reads_count_all_projects_tpm,  paste(output_dir,"df_reads_count_all_projects_tpm.tsv",sep="/"), na = "NA", append = TRUE, col.names = TRUE, row.names = TRUE, sep = "\t", quote = TRUE)
###########################################################################################################################




