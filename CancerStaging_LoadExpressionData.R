#####################################################################################################################
# This script will take the TSV file (metadata), unstranded.rna_seq.augmented_star_gene_counts (rna-seq count data), 
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
# Remove genes whose counts are zero in at least one patient
unstranded_raw_data<-unstranded_raw_data[rowSums(unstranded_raw_data[])>0,]

# Genes whose counts are zero in one patient
nonzero_genes<-rownames(unstranded_raw_data)
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
    project_data_raw<-unstranded_raw_data[nonzero_genes,which(colnames(unstranded_raw_data) %in% sample_ids)]
    project_data_tpm<-unstranded_tpm_data[nonzero_genes,which(colnames(unstranded_tpm_data) %in% sample_ids)]
    project_data_fpkm<-unstranded_fpkm_data[nonzero_genes,which(colnames(unstranded_fpkm_data) %in% sample_ids)]

    # Store dataset
    reads_count_per_project_raw[[project]]<-project_data_raw
    reads_count_per_project_tpm[[project]]<-project_data_tpm
    reads_count_per_project_fpkm[[project]]<-project_data_fpkm
}
#####################################################################################################################
# If TCGA_project is not all
if (TCGA_project != "ALL")
{
    # Count the number of reads per project, raw, tpm, tpkm
    dim(reads_count_per_project_raw[[TCGA_project]])
    dim(reads_count_per_project_tpm[[TCGA_project]])
    dim(reads_count_per_project_fpkm[[TCGA_project]])

    reads_count_per_project_raw<-reads_count_per_project_raw[[TCGA_project]]
    reads_count_per_project_tpm<-reads_count_per_project_tpm[[TCGA_project]]
    reads_count_per_project_fpkm<-reads_count_per_project_fpkm[[TCGA_project]]

    # Count the number of reads per project
    df_reads_count_all_projects_raw<-reads_count_per_project_raw
    
    # Count the number of reads per project
    df_reads_count_all_projects_fpkm<-reads_count_per_project_fpkm
    
    # Count the number of reads per project
    df_reads_count_all_projects_tpm<-reads_count_per_project_tpm

    # Take name of all samples
    all_samples<-c(colnames(reads_count_per_project_raw))        
}else
{    
    # Take name of all samples
    all_samples<-c(colnames(reads_count_per_project_raw[["TCGA-BRCA"]]),
    colnames(reads_count_per_project_raw[["TCGA-LIHC"]]),
    colnames(reads_count_per_project_raw[["TCGA-LUAD"]]),
    colnames(reads_count_per_project_raw[["TCGA-LUSC"]]),
    colnames(reads_count_per_project_raw[["TCGA-PRAD"]]),
    colnames(reads_count_per_project_raw[["TCGA-READ"]]),
    colnames(reads_count_per_project_raw[["TCGA-SKCM"]]),
    colnames(reads_count_per_project_raw[["TCGA-STAD"]]),
    colnames(reads_count_per_project_raw[["TCGA-THCA"]]))
                   
    # Split metadata per project
    merged_data_patient_info_BRCA<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-BRCA",]
    merged_data_patient_info_LIHC<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-LIHC",]
    merged_data_patient_info_LUAD<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-LUAD",]
    merged_data_patient_info_LUSC<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-LUSC",]
    merged_data_patient_info_PRAD<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-PRAD",]
    merged_data_patient_info_READ<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-READ",]
    merged_data_patient_info_SKCM<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-SKCM",]
    merged_data_patient_info_STAD<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-STAD",]    
    merged_data_patient_info_STAD<-merged_data_patient_info[merged_data_patient_info$project_id == "TCGA-THCA",]    

    # Count the number of reads per project
    df_reads_count_all_projects_raw<-cbind(reads_count_per_project_raw[["TCGA-BRCA"]],
    reads_count_per_project_raw[["TCGA-LIHC"]],
    reads_count_per_project_raw[["TCGA-LUAD"]],
    reads_count_per_project_raw[["TCGA-LUSC"]],
    reads_count_per_project_raw[["TCGA-PRAD"]],
    reads_count_per_project_raw[["TCGA-READ"]],
    reads_count_per_project_raw[["TCGA-SKCM"]],
    reads_count_per_project_raw[["TCGA-STAD"]],
    reads_count_per_project_raw[["TCGA-THCA"]])
    
    # Count the number of reads per project
    df_reads_count_all_projects_fpkm<-cbind(reads_count_per_project_fpkm[["TCGA-BRCA"]],
    reads_count_per_project_fpkm[["TCGA-LIHC"]],
    reads_count_per_project_fpkm[["TCGA-LUAD"]],
    reads_count_per_project_fpkm[["TCGA-LUSC"]],
    reads_count_per_project_fpkm[["TCGA-PRAD"]],
    reads_count_per_project_fpkm[["TCGA-READ"]],
    reads_count_per_project_fpkm[["TCGA-SKCM"]],
    reads_count_per_project_fpkm[["TCGA-STAD"]],
    reads_count_per_project_raw[["TCGA-THCA"]])
    
    # Count the number of reads per project
    df_reads_count_all_projects_tpm<-cbind(reads_count_per_project_tpm[["TCGA-BRCA"]],
    reads_count_per_project_tpm[["TCGA-LIHC"]],
    reads_count_per_project_tpm[["TCGA-LUAD"]],
    reads_count_per_project_tpm[["TCGA-LUSC"]],
    reads_count_per_project_tpm[["TCGA-PRAD"]],
    reads_count_per_project_tpm[["TCGA-READ"]],
    reads_count_per_project_tpm[["TCGA-SKCM"]],
    reads_count_per_project_tpm[["TCGA-STAD"]],
    reads_count_per_project_raw[["TCGA-THCA"]])
    
}
# Use only data that has read counts
merged_data_patient_info      <-unique(merged_data_patient_info[which( merged_data_patient_info$sample_id %in%  all_samples),])   
merged_data_patient_info_count<-unique(merged_data_patient_info[which( merged_data_patient_info$sample_id %in%  all_samples),c("project_id","stages", "Sample.ID","sample_id", "Sample.Type")])   

# Cases per stage
table_cases_per_stage<-table(merged_data_patient_info_count$project_id)

# Write to to file the number of cases
cat(paste("Number of cases: ", sep=""),file=results_files,sep="\n", append=TRUE)
cat(paste(names(table_cases_per_stage), collapse=" "),file=results_files,sep="\n", append=TRUE)
cat(paste(table_cases_per_stage, collapse=" "),file=results_files,sep="\n", append=TRUE)
#####################################################################################################################
# Organize how to send to Carles
write_tsv(merged_data_patient_info, "/home/felipe/Documents/Cancer_staging/merged_data_patient_info.tsv")
#####################################################################################################################
# Total number of samples
#TCGA-BRCA TCGA-LIHC TCGA-LUAD TCGA-LUSC TCGA-PRAD TCGA-READ TCGA-SKCM TCGA-STAD 
#     1180       411       579       536       545       167       104       440 
###########################################################################################################################
#write.table(df_reads_count_all_projects_raw,  paste(output_dir,"df_reads_count_all_projects_raw.tsv",sep="/"), na = "NA", append = TRUE, col.names = TRUE, row.names = TRUE, sep = "\t", quote = TRUE)
#write.table(df_reads_count_all_projects_fpkm, paste(output_dir,"df_reads_count_all_projects_fpkm.tsv",sep="/"), na = "NA", append = TRUE, col.names = TRUE, row.names = TRUE, sep = "\t", quote = TRUE)
#write.table(df_reads_count_all_projects_tpm,  paste(output_dir,"df_reads_count_all_projects_tpm.tsv",sep="/"), na = "NA", append = TRUE, col.names = TRUE, row.names = TRUE, sep = "\t", quote = TRUE)
###########################################################################################################################
