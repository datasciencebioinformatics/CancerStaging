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
merged_data_patient_info$patient_id  <- merged_data_patient_info$File.ID
merged_data_patient_info$case_id     <- merged_data_patient_info$Case.ID
#merged_data_patient_info$sample_id  <- merged_data_patient_info$Sample.ID
###########################################################################################################################
# Paired samples                                                                                                         
paired_sample_df<-data.frame(normal=c(),tumor=c(),case=c(),project=c())                                                              
                                                                                                                         
# For each case, find the pairs                                                                                          
for (case in unique(merged_data_patient_info$Case.ID))                                                                 
{                                                                                                                        
    # All samples for case id = "case"                                                                                   
    case_samples<-merged_data_patient_info[merged_data_patient_info$Case.ID==case,]             
                                                                                                                         
    # Take the tumor samples                                                                                           
    tumor_sampĺes <-case_samples[case_samples$tissue_type=="Tumor",]
    normal_sampĺes<-case_samples[case_samples$tissue_type=="Normal",]

    # Project id
    project_id<-unique(case_samples[case_samples$Case.ID==case,"project_id"])                               

    # if vector contains at least one tumor and one normal
    if(length(unique(normal_sampĺes$sample_id))>0 && length(unique(tumor_sampĺes$sample_id))>0)
    {            
            # For each tumor sample
            for (tumor_solid_sample_id in tumor_sampĺes$sample_id)
            {
                # for each normal sample, compile a paired samples
                for (normal_samples_id in normal_sampĺes$sample_id)
                {                  
                    # Contatenate                     
                    paired_sample_df<-rbind(data.frame(normal=c(normal_samples_id),tumor=c(tumor_solid_sample_id),case=case, project=project_id),paired_sample_df)
                }
            }                
    }
}
#######################################################################################################################################
paired_sample_df<-unique(paired_sample_df)
#######################################################################################################################################
# table(paired_sample_df$project)
#TCGA-BRCA TCGA-LIHC TCGA-LUAD TCGA-LUSC TCGA-PRAD TCGA-READ TCGA-STAD 
#      119        50        70        51        54         9        33
#######################################################################################################################################
# Verify paired sample
# Verifgy paired samples
# One case, paired tumor control
# One tumor - one control
#######################################################################################################################################
# Save normalized data                                                                                             
merged_data_patient_info_count<-merged_data_patient_info_count[merged_data_patient_info_count$sample_id %in% colnames(df_reads_count_all_projects_fpkm),]
###########################################################################################################################
# All tumor and control samples
unpaired_tumor_samples  <-merged_data_patient_info_count[merged_data_patient_info_count$Sample.Type=="Primary Tumor","sample_id"]
unpaired_control_samples<-merged_data_patient_info_count[merged_data_patient_info_count$Sample.Type=="Solid Tissue Normal","sample_id"]

# Paired samples only
paired_normal_samples    <- paired_sample_df$normal
paired_tumor_samples     <- paired_sample_df$tumor
#######################################################################################################################################
# logchange_tumor_control
list_logchange_tumor_control<-list()

# for each normalization scheme
for (normalized_table_names in names(df_reads_count_all_projects))
{  
  # Store normalized table
  normalized_table<-df_reads_count_all_projects[[normalized_table_names]]

  # folchange=Expr(Stage i)/Expr(Stage ii and II)
  # Paired t-test, RPKM of paired tumor/normal samples
  # Plot with 15208 genes.
  # Log2foldchange
  LOG_CONSTANT=0.001
  log2change       =log( (rowMeans(normalized_table[,unpaired_tumor_samples]+LOG_CONSTANT)/rowMeans(normalized_table[,unpaired_control_samples]+LOG_CONSTANT)),2)	
  log2change_paired=log( (rowMeans(normalized_table[,paired_tumor_samples]+LOG_CONSTANT)/rowMeans(normalized_table[,paired_normal_samples]+LOG_CONSTANT)),2)	
  
  # log2change data
  log2change_tumor_control=na.omit(data.frame(gene=names(log2change),log2change=log2change))
  log2change_tumor_control_paired=na.omit(data.frame(gene=names(log2change_paired),log2change=log2change_paired))
  
  # First, the log2foldchane tumor/normal samples is used
  log2change_tumor_control$Pvalue<-1
  log2change_tumor_control_paired$Pvalue<-1
  
  # For each genes in the tabe
  for (gene in log2change_tumor_control$gene)
  {

	# Take the expression of genes above expression threshold for Stage_i_samples
	expr_unpaired_tumor_samples<-normalized_expression_table[which(normalized_expression_table[gene,unpaired_tumor_samples]>list_threshold_filters[[normalization_scheme]]),]
	
	# Take the expression of genes above expression threshold for Stage_i_samples
	expr_unpaired_control_samples<-normalized_expression_table[which(normalized_expression_table[gene,unpaired_control_samples]>list_threshold_filters[[normalization_scheme]]),]			
	
	# Take the expression of genes above expression threshold for Stage_i_samples
	expr_paired_tumor_samples<-normalized_expression_table[which(normalized_expression_table[gene,paired_tumor_samples]>list_threshold_filters[[normalization_scheme]]),]
	
	# Take the expression of genes above expression threshold for Stage_i_samples
	expr_paired_control_samples<-normalized_expression_table[which(normalized_expression_table[gene,paired_control_samples]>list_threshold_filters[[normalization_scheme]]),]			    
    
	# Take p-value
	log2change_tumor_control[gene,"Pvalue"]<-t.test(x=as.numeric(expr_unpaired_tumor_samples), y=as.numeric(expr_unpaired_control_samples), paired = FALSE, alternative = "two.sided")$p.value	
	log2change_tumor_control_paired[gene,"Pvalue"]<-t.test(x=as.numeric(expr_paired_tumor_samples), y=as.numeric(expr_paired_control_samples), paired = TRUE, alternative = "two.sided")$p.value	
  }
  # FRD 
  log2change_tumor_control$FDR<-p.adjust(log2change_tumor_control$Pvalue, method="fdr")
  log2change_tumor_control_paired$FDR<-p.adjust(log2change_tumor_control_paired$Pvalue, method="fdr")
  #######################################################################################################################################
  colnames(log2change_tumor_control)       <- c("gene","log2change_all_samples","pvalue_all_samples","fdr_all_samples")                 #
  colnames(log2change_tumor_control_paired)<- c("gene","log2change_paired","pvalue_paired","fdr_paired")                                #
  #######################################################################################################################################
  logchange_tumor_control<-merge(log2change_tumor_control,log2change_tumor_control_paired,by="gene")                                    #
  #######################################################################################################################################
  list_logchange_tumor_control[[normalized_table_names]]<-logchange_tumor_control                                                       #
  #######################################################################################################################################  
  # Write TSV
  write_tsv(data.frame(logchange_tumor_control), paste(output_dir,"df_statistics_all_projects_",normalized_table_names,".tsv",sep="")) #
  ####################################################################################################################################### 
}
print("\nCancerStaging_CreatePairedSamplesTumorNormal")

