# A R script to create metadata from gdc files.
# Inputs:
# gdc_sample_sheet.2024-03-08.tsv
# clinical_"database".txt
# sample_"database".txt
# exposure_"database".txt
# merged_data_patient_info_"database".txt.tsv
##########################################################################################################################################################################################################
# Reading the contents of TSV file using read_tsv() method
gdc_sample_sheet_file<-"/home/felipe/Documents/Cancer_staging/gdc_sample_sheet.2024-07-18.tsv"

# Read data
gdc_sample_sheet_data<-read.table(file = gdc_sample_sheet_file, sep = '\t', header = TRUE,fill=TRUE)  

# Add collumn sample_id
gdc_sample_sheet_data$sample_submitter_id<-gdc_sample_sheet_data$Sample.ID
#####################################################################################################################
# Set path to files                                                                                                 
clinical_file="/home/felipe/Documents/Cancer_staging/clinical.txt" 
sample_file="/home/felipe/Documents/Cancer_staging/sample.txt"    
exposure_file="/home/felipe/Documents/Cancer_staging/exposure.txt"                                                  

# Load data
clinical_data<-read.table(file = clinical_file, sep = '\t', header = TRUE,fill=TRUE)    
sample_data<-read.table(file = sample_file, sep = '\t', header = TRUE,fill=TRUE)                                    
exposure_data<-read.table(file = exposure_file, sep = '\t', header = TRUE,fill=TRUE)                                #

# Merge data
merged_sample_clinical_data<-merge(sample_data,clinical_data,by="case_id")

# Merge all
merged_sample_clinical_data<-merge(merged_sample_clinical_data,exposure_data,by="case_id")

# Merge tables
merged_data_patient_info<-merge(merged_sample_clinical_data,gdc_sample_sheet_data,by="sample_submitter_id")
#####################################################################################################################
# Set file name variable 
merged_data_patient_info<-merged_data_patient_info[merged_data_patient_info$Data.Category=="Transcriptome Profiling",]

# Check which entries contains the words .rna_seq.augmented_star_gene_counts.tsv
merged_data_patient_info<-merged_data_patient_info[which(grepl(pattern="*.rna_seq.augmented_star_gene_counts.tsv", x=merged_data_patient_info$File.Name)),]

# From the File.ID, only the ID is kept in the variable sample_id
merged_data_patient_info$sample_id<-gsub(".rna_seq.augmented_star_gene_counts.tsv", "", merged_data_patient_info$File.Name)
#####################################################################################################################
# length(unique(merged_data_patient_info$case_id)) # Number of cases
# length(unique(merged_data_patient_info$sample_id)) # Number of samples
# sum(unique(merged_data_patient_info[,c("sample_id","Sample.Type")])[,2]=="Primary Tumor") # Number of Primary Tumor
# sum(unique(merged_data_patient_info[,c("sample_id","Sample.Type")])[,2]=="Solid Tissue Normal") # Number of Solid Tissue Normal

# Filter tumor and normal samples
primary_tumor<-merged_data_patient_info[merged_data_patient_info[,c("sample_id","Sample.Type")][,2]=="Primary Tumor",]
solid_tissue<-merged_data_patient_info[merged_data_patient_info[,c("sample_id","Sample.Type")][,2]=="Solid Tissue Normal",]
merged_data_patient_info<-rbind(primary_tumor,solid_tissue)

# length(unique(primary_tumor$sample_id))
# length(unique(solid_tissue$sample_id))

# Population demographic
# table(unique(merged_data_patient_info[,c("sample_id","primary_diagnosis")])$primary_diagnosis)
#merged_data_patient_info<-merged_data_patient_info[merged_data_patient_info$primary_diagnosis=="Squamous cell carcinoma, NOS",]

# Filter tumor and normal samples
primary_tumor<-merged_data_patient_info[merged_data_patient_info[,c("sample_id","Sample.Type")][,2]=="Primary Tumor",]
solid_tissue<-merged_data_patient_info[merged_data_patient_info[,c("sample_id","Sample.Type")][,2]=="Solid Tissue Normal",]
merged_data_patient_info<-rbind(primary_tumor,solid_tissue)

# table(unique(merged_data_patient_info[,c("sample_id","ethnicity")])$ethnicity)
# table(unique(merged_data_patient_info[,c("sample_id","gender")])$gender)
# table(unique(merged_data_patient_info[,c("sample_id","vital_status")])$vital_status)
# min(merged_data_patient_info[!is.na(merged_data_patient_info$age_at_index),"age_at_index"])
# max(merged_data_patient_info[!is.na(merged_data_patient_info$age_at_index),"age_at_index"])
# mean(merged_data_patient_info[!is.na(merged_data_patient_info$age_at_index),"age_at_index"])
#####################################################################################################################
# Summary of counts of samples per stage
# Number of Primary Tumor
# Number of Solid Tissue Normal
sum(unique(merged_data_patient_info[,c("sample_id","Sample.Type")])[,2]=="Primary Tumor") # Number of Primary Tumor
sum(unique(merged_data_patient_info[,c("sample_id","Sample.Type")])[,2]=="Solid Tissue Normal") # Number of Solid Tissue Normal

# A field to store 
merged_data_patient_info$stages<-merged_data_patient_info$ajcc_pathologic_stage

# Group stages I,II,III and IV
merged_data_patient_info$stages<-gsub("Stage IA", "Stage I", merged_data_patient_info$stages)
merged_data_patient_info$stages<-gsub("Stage IB", "Stage I", merged_data_patient_info$stages)
merged_data_patient_info$stages<-gsub("Stage IIA", "Stage II", merged_data_patient_info$stages)
merged_data_patient_info$stages<-gsub("Stage IIB", "Stage II", merged_data_patient_info$stages)
merged_data_patient_info$stages<-gsub("Stage IIC", "Stage II", merged_data_patient_info$stages)
merged_data_patient_info$stages<-gsub("Stage IIIA", "Stage III", merged_data_patient_info$stages)
merged_data_patient_info$stages<-gsub("Stage IIIB", "Stage III", merged_data_patient_info$stages)
merged_data_patient_info$stages<-gsub("Stage IIIC", "Stage III", merged_data_patient_info$stages)
merged_data_patient_info$stages<-gsub("Stage IVA", "Stage IV", merged_data_patient_info$stages)
merged_data_patient_info$stages<-gsub("Stage IVB", "Stage IV", merged_data_patient_info$stages)

# Count samples
merged_data_patient_count<-unique(merged_data_patient_info[,c("sample_id","Sample.Type","project_id","stages")])

# Cases per stage
table_cases_per_stage<-table(merged_data_patient_count$project_id)

#####################################################################################################################
# Total of 4075 samples
#TCGA-BRCA TCGA-LIHC TCGA-LUAD TCGA-LUSC TCGA-PRAD TCGA-READ TCGA-SKCM TCGA-STAD 
#     1223       421       598       553       553       175       104       448 
#####################################################################################################################
print("\nCancerStaging_CreateMetadataFromGDCFiles")

# Organize how to send to Carles
write_tsv(merged_data_patient_info, paste(output_dir,"merged_data_patient_info.tsv",sep="/"))
