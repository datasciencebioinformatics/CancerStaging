###########################################################################################################################
# Comments : from the GDC metadata, samples from each case are split according to tissue_type=Normal,tissue_type=Tumor.
# the paired samples are combinations of tissue_type=Normal,tissue_type=Tumor from the same case.
###########################################################################################################################
# A R script to create metadata from gdc files.
# Inputs:
# gdc_sample_sheet.2024-03-08.tsv
# /home/felipe/Documentos/LungPortal/clinical.txt
# /home/felipe/Documentos/LungPortal/sample.txt
# /home/felipe/Documentos/LungPortal/exposure.txt
# Output : merged_data_patient_info.tsv
###########################################################################################################################
merged_data_patient_info_file       <- "/home/felipe/googledrive/Cancer_staging/merged_data_patient_info.tsv"             #
###########################################################################################################################
merged_data_patient_info_data      <-read.table(file = merged_data_patient_info_file, sep = '\t', header = TRUE,fill=TRUE)#
###########################################################################################################################
merged_data_patient_info_data$patient_id  <- merged_data_patient_info_data$File.ID
merged_data_patient_info_data$case_id     <- merged_data_patient_info_data$Case.ID
merged_data_patient_info_data$sample_id   <- merged_data_patient_info_data$Sample.ID
###########################################################################################################################
# Paired samples                                                                                                         
paired_sample_df<-data.frame(normal=c(),tumor=c(),case=c(),project=c())                                                              
                                                                                                                         
# For each case, find the pairs                                                                                          
for (case in unique(merged_data_patient_info_data$Case.ID))                                                                 
{                                                                                                                        
    # All samples for case id = "case"                                                                                   
    case_samples<-merged_data_patient_info_data[merged_data_patient_info_data$Case.ID==case,]                               
                                                                                                                         
    # Take the tumor samples                                                                                           
    tumor_sampĺes <-case_samples[case_samples$tissue_type=="Tumor",]
    normal_sampĺes<-case_samples[case_samples$tissue_type=="Normal",]

    # Take the project id
    stored_project_id<-unique(case_samples$project_id)

    # if vector contains at least one tumor and one normal
    if(length(unique(normal_sampĺes$sample_id))>0 && length(unique(tumor_sampĺes$sample_id))>0)
    {            
            # For each tumor sample
            for (tumor_solid_sample_id in tumor_sampĺes$patient_id)
            {
                # for each normal sample, compile a paired samples
                for (normal_samples_id in normal_sampĺes$patient_id)
                {                  
                    # Contatenate                     
                    paired_sample_df<-rbind(data.frame(normal=c(normal_samples_id),tumor=c(tumor_solid_sample_id),case=case, project=stored_project_id),paired_sample_df)
                }
            }                
    }
}
#######################################################################################################################################
# Verify paired sample
# Verifgy paired samples
# One case, paired tumor control
# One tumor - one control
#######################################################################################################################################
