# A R script to create metadata from gdc files.
# Inputs:
# gdc_sample_sheet.2024-03-08.tsv
# clinical_"database".txt
# sample_"database".txt
# exposure_"database".txt
# merged_data_patient_info_"database".txt.tsv
##########################################################################################################################################################################################################
# Reading the contents of TSV file using read_tsv() method
gdc_sample_sheet_file<-"/home/felipe/googledrive/Cancer_staging/gdc_sample_sheet.2024-07-18.tsv"

# Read data
gdc_sample_sheet_data<-read.table(file = gdc_sample_sheet_file, sep = '\t', header = TRUE,fill=TRUE)  

# Add collumn sample_id
gdc_sample_sheet_data$sample_submitter_id<-gdc_sample_sheet_data$Sample.ID
#####################################################################################################################
# Set file name variable 
gdc_sample_sheet_data<-gdc_sample_sheet_data[gdc_sample_sheet_data$Data.Category=="Transcriptome Profiling",]

# Check which entries contains the words .rna_seq.augmented_star_gene_counts.tsv
gdc_sample_sheet_data<-gdc_sample_sheet_data[which(grepl(pattern="*.rna_seq.augmented_star_gene_counts.tsv", x=gdc_sample_sheet_data$File.Name)),]

# From the File.ID, only the ID is kept in the variable sample_id
gdc_sample_sheet_data$sample_id<-gsub(".rna_seq.augmented_star_gene_counts.tsv", "", gdc_sample_sheet_data$File.Name)


#####################################################################################################################
# Set path to files                                                                                                 
# Check ajcc_pathologic_stage
# primary_diagnosis
# site_of_resection_or_biopsy
#cat /home/felipe/googledrive/Cancer_staging/clinical.tsv | sed "s/'--/ /g" > /home/felipe/googledrive/Cancer_staging/clinical.txt
#cat /home/felipe/googledrive/Cancer_staging/sample.tsv | sed "s/'--/ /g" > /home/felipe/googledrive/Cancer_staging/sample.txt
#cat /home/felipe/googledrive/Cancer_staging/exposure.tsv | sed "s/'--/ /g" > /home/felipe/googledrive/Cancer_staging/exposure.txt

# Set path to files                                                                                                 
clinical_file="/home/felipe/googledrive/Cancer_staging/clinical.txt" 
sample_file="/home/felipe/googledrive/Cancer_staging/sample.txt"    
exposure_file="/home/felipe/googledrive/Cancer_staging/exposure.txt"                                                  

# Load data
clinical_data<-read.table(file = clinical_file, sep = '\t', header = TRUE,fill=TRUE)    
sample_data<-read.table(file = sample_file, sep = '\t', header = TRUE,fill=TRUE)                                    
exposure_data<-read.table(file = exposure_file, sep = '\t', header = TRUE,fill=TRUE)                                #

#####################################################################################################################
# Create field merge_id
gdc_sample_sheet_data$merge_id<-gdc_sample_sheet_data$Case.ID 
clinical_file$merge_id        <-clinical_file$case_submitter_id
#####################################################################################################################

# Merge data
merged_sample_clinical_data<-merge(sample_data,clinical_data,by="case_id")

# Merge all
merged_sample_clinical_data<-merge(merged_sample_clinical_data,exposure_data,by="case_id")

# Merge tables
merged_data_patient_info<-merge(merged_sample_clinical_data,gdc_sample_sheet_data,by.x="sample_submitter_id",by.y="sample_submitter_id" )
#####################################################################################################################
sum(unique(merged_data_patient_info[,c("sample_id.x","Sample.Type")])[,2]=="Primary Tumor") # Number of Primary Tumor
sum(unique(merged_data_patient_info[,c("sample_id.x","Sample.Type")])[,2]=="Solid Tissue Normal") # Number of Solid Tissue Normal
sum(unique(merged_data_patient_info[,c("sample_id.x","Sample.Type")])[,2]=="Metastatic") # Number of Solid Tissue Normal

# A field to store 
merged_data_patient_info$stages<-merged_data_patient_info$ajcc_pathologic_stage

# Group stages I,II,III and IV
merged_data_patient_info$stages<-gsub("Stage IA", "Stage I", merged_data_patient_info$stages)
merged_data_patient_info$stages<-gsub("Stage IB", "Stage I", merged_data_patient_info$stages)
merged_data_patient_info$stages<-gsub("Stage IIA", "Stage II", merged_data_patient_info$stages)
merged_data_patient_info$stages<-gsub("Stage IIB", "Stage II", merged_data_patient_info$stages)
merged_data_patient_info$stages<-gsub("Stage IIC", "Stage II", merged_data_patient_info$stages)
merged_data_patient_info$stage<-gsub("Stage IIIA", "Stage III", merged_data_patient_info$stages)
merged_data_patient_info$stages<-gsub("Stage IIIB", "Stage III", merged_data_patient_info$stages)
merged_data_patient_info$stages<-gsub("Stage IIIC", "Stage III", merged_data_patient_info$stages)
merged_data_patient_info$stages<-gsub("Stage IVA", "Stage IV", merged_data_patient_info$stages)
merged_data_patient_info$stages<-gsub("Stage IVB", "Stage IV", merged_data_patient_info$stages)

# Cases per stage
table_cases_per_stage<-table(merged_data_patient_info$project_id, merged_data_patient_info$stages)

# Cases per stage
table_cases_per_stage<-table_cases_per_stage[,c("Stage I","Stage II","Stage III","Stage IV")]

# Organize how to send to Carles
write_tsv(merged_data_patient_info, "/home/felipe/googledrive/Cancer_staging/merged_data_patient_info.tsv")
#####################################################################################################################
write_tsv(data.frame(merged_data_patient_info$File.Name), "/home/felipe/googledrive/Cancer_staging/used_file_names.tsv")
#####################################################################################################################
