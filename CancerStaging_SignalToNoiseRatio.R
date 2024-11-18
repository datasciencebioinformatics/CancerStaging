###################################################################################################################################################
source("/home/felipe/Documents/github/CancerStaging/CancerStaging_SetupAllParamters.R")
source("/home/felipe/Documents/github/CancerStaging/CancerStaging_LoadRPackages.R")
source("/home/felipe/Documents/github/CancerStaging/CancerStaging_CreateMetadataFromGDCFiles.R")
source("/home/felipe/Documents/github/CancerStaging/CancerStaging_LoadInteractomeData.R")
source("/home/felipe/Documents/github/CancerStaging/CancerStaging_LoadExpressionData.R")
source("/home/felipe/Documents/github/CancerStaging/CancerStaging_ExpressionDataNormalization_Carels.R")
source("/home/felipe/Documents/github/CancerStaging/CancerStaging_CreatePairedSamplesTumorNormal.R")
source("/home/felipe/Documents/github/CancerStaging/CancerStaging_FindTumorGenes.R")
###################################################################################################################################################
#Moreover, through investigation of biomarkers in up-regulated genes, we found 52 genes with fold changes larger than 50 in tumor, and average TPM ≤ 10 (Table S1). In table 2, we present the biomarkers whose fold change is ≥100 and average TPM ≤ 2, which should warrant significant signal-to-noise ratio (Figure 2). The markers that met these criteria had average T`PM ≥ 100.
normalization_scheme<-"tpm"
###################################################################################################################################################
# All tumor and control samples
colData_tumor  <-unique(merged_data_patient_info_count[merged_data_patient_info_count$Sample.Type=="Primary Tumor",])
colData_normal <-unique(merged_data_patient_info_count[merged_data_patient_info_count$Sample.Type=="Solid Tissue Normal",])
###################################################################################################################################################

# Samples of each stage stored in colData                                                                                             #
sample_stage_I  <-unique(colData_tumor[colData_tumor$stages=="Stage I","sample_id"])                                                          #
sample_stage_II <-unique(colData_tumor[colData_tumor$stages=="Stage II","sample_id"])                                                         #
sample_stage_III<-unique(colData_tumor[colData_tumor$stages=="Stage III","sample_id"])                                                        #
sample_normal   <-unique(colData_normal[,"sample_id"])

normal_samples<-sample_normal
tumor_samples <-c(sample_stage_I,sample_stage_II,sample_stage_III)

# Take the average expression
# Select only tumor genes
normalized_expression_table<-df_reads_count_all_projects[["tpm"]]

# Calculate rowmeans
rowMeans_normalized_expression_table<-rowMeans(normalized_expression_table[,normal_samples])
rowMeans_normalized_expression_table_tumor<-rowMeans(normalized_expression_table[,tumor_samples])


# All tumor and control samples
colData_tumor  <-unique(merged_data_patient_info_count[merged_data_patient_info_count$Sample.Type=="Primary Tumor",])
colData_normal <-unique(merged_data_patient_info_count[merged_data_patient_info_count$Sample.Type=="Solid Tissue Normal",])

tumor_samples<-colData_tumor$sample_id
normal_samples<-colData_normal$sample_id

# Data.frame
df_FC<-data.frame(Gene=c(),Mean_normal=c(),sd_normal=c(),Mean_tumor=c(),sd_tumor=c(),FC=c(),Mean_Stage_I=c(),
 sd_Stage_I=c(),Mean_Stage_II=c(), sd_Stage_II=c(), Mean_Stage_II=c(), sd_Stage_III=c())
 
# For each genes in the tabe
for (gene in names(rowMeans_normalized_expression_table))
{
	# Take expression only for the gene
	mean_Stage_tumor <-mean(as.vector(unlist(normalized_expression_table[gene,tumor_samples])))
	mean_Stage_normal<-mean(as.vector(unlist(normalized_expression_table[gene,normal_samples])))
	mean_Stage_I   <-mean(as.vector(unlist(normalized_expression_table[gene,sample_stage_I])))
	mean_Stage_II  <-mean(as.vector(unlist(normalized_expression_table[gene,sample_stage_II])))
	mean_Stage_III <-mean(as.vector(unlist(normalized_expression_table[gene,sample_stage_III])))
	
	# Take expression only for the gene
	sd_Stage_tumor <-sd(as.vector(unlist(normalized_expression_table[gene,tumor_samples])))
	sd_Stage_normal<-sd(as.vector(unlist(normalized_expression_table[gene,normal_samples])))
	sd_Stage_I     <-sd(as.vector(unlist(normalized_expression_table[gene,sample_stage_I])))
	sd_Stage_II    <-sd(as.vector(unlist(normalized_expression_table[gene,sample_stage_II])))
	sd_Stage_III   <-sd(as.vector(unlist(normalized_expression_table[gene,sample_stage_III])))	
	
	# Filter by threshold_filters
	#mean_Stage_tumor<-mean_Stage_tumor[mean_Stage_tumor > list_threshold_filters[[normalization_scheme]]]

	# Filter by threshold_filters
	#mean_Stage_normal<-mean_Stage_normal[mean_Stage_normal > list_threshold_filters[[normalization_scheme]]]
	
	# Filter by threshold_filters
	#mean_Stage_I<-mean_Stage_I[mean_Stage_I > list_threshold_filters[[normalization_scheme]]]
	#mean_Stage_II<-mean_Stage_II[mean_Stage_II > list_threshold_filters[[normalization_scheme]]]	
	#mean_Stage_III<-mean_Stage_III[mean_Stage_III > list_threshold_filters[[normalization_scheme]]]
	
	FC=mean_Stage_tumor/mean_Stage_normal	
	
	if(length(mean_Stage_tumor) && length(mean_Stage_normal) && length(mean_Stage_I) && length(mean_Stage_II) && length(mean_Stage_III) )
	{		
		df_FC<-rbind(df_FC,data.frame(Gene=gene,Mean_normal=mean_Stage_normal,sd_normal=sd_Stage_normal,Mean_tumor=mean_Stage_I,sd_tumor=sd_Stage_tumor,FC=FC,Mean_Stage_I=mean_Stage_I,
	 sd_Stage_I=sd_Stage_I,Mean_Stage_II=mean_Stage_II, sd_Stage_II=sd_Stage_II, Mean_Stage_III=mean_Stage_III, sd_Stage_III=sd_Stage_III))			
	}	
}
# Set rownames
rownames(df_FC)<-df_FC$Gene

selected_genes<-df_FC[df_FC$FC>=50 & df_FC$Mean_normal<=10,]

# Save TSV file with genes from Stage3
write_tsv(selected_genes, paste(output_dir,"/selected_genes.tsv",sep=""))

# Take selected genes
selected_genes<-df_FC[which(df_FC$FC>100 & df_FC$Mean_normal<=10),]
