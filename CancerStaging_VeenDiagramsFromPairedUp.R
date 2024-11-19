####################################################################################################################################################
# Interactomes_GC3_T2.csv file has 15650 entries. The number of annotated genes with gene length geneLength_ENTREZID_ENSEMBL is 14609. Among these, 14726 are common to Interactomes_GC3_T2 and geneLength_ENTREZID_ENSEMBL and will be used to create the maps. 
# Consitency - check filters meticulously.
# FPKM, TPM  - take these as robust.
# Paramter to set the normalization_scheme
normalization_schemes<-c("tpm","fpkm","tmm","rpkm","tpm_calc")
normalization_schemes<-c("tpm")

# For each normlization normalization_scheme
for (normalization_scheme in normalization_schemes)
{
  #######################################################################################################################################
  # Path to files of selected_genes                                                                                                             # 
  # genes_stages_I
  selected_genes_Stage_I_data    <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_I",".tsv",sep=""), sep = '\t', header = TRUE) #
  selected_genes_Stage_II_data   <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_II",".tsv",sep=""), sep = '\t', header = TRUE) #
  selected_genes_Stage_III_data  <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_III",".tsv",sep=""), sep = '\t', header = TRUE) #

  rownames(selected_genes_Stage_I_data)<-selected_genes_Stage_I_data$gene
  rownames(selected_genes_Stage_II_data)<-selected_genes_Stage_II_data$gene
  rownames(selected_genes_Stage_III_data)<-selected_genes_Stage_III_data$gene

  selected_genes_Stage_I_gene      <- selected_genes_Stage_I_data$gene
  selected_genes_Stage_II_gene     <- selected_genes_Stage_II_data$gene
  selected_genes_Stage_III_gene    <- selected_genes_Stage_III_data$gene
  #######################################################################################################################################                                                                                                                                     #
  unique_stage_I  =intersect(setdiff(selected_genes_Stage_I_gene, c(selected_genes_Stage_II_gene,selected_genes_Stage_III_gene)),selected_genes_Stage_I_gene)
  unique_stage_II =intersect(setdiff(selected_genes_Stage_II_gene, c(selected_genes_Stage_I_gene,selected_genes_Stage_III_gene)),selected_genes_Stage_II_gene)
  unique_stage_III=intersect(setdiff(selected_genes_Stage_III_gene, c(selected_genes_Stage_I_gene,selected_genes_Stage_II_gene)),selected_genes_Stage_III_gene)
  #######################################################################################################################################
  # Select stages 
  stages_I_II_III_unique<-ggVennDiagram(list(Stage_I=unique_stage_I,Stage_II=unique_stage_II,Stage_III=unique_stage_III), label_alpha = 0.9,set_color = c("grey50","grey50","grey50")) +  scale_fill_gradient(low = "white", high = "white") + theme_bw() + ggtitle("Stages I, II and III")+ guides(fill="none")
  stages_I_II_III_unique<-ggVennDiagram(list(Stage_I=selected_genes_Stage_I_gene,Stage_II=selected_genes_Stage_II_gene,Stage_III=selected_genes_Stage_III_gene), label_alpha = 0.9,set_color = c("grey50","grey50","grey50")) +  scale_fill_gradient(low = "white", high = "white") + theme_bw() + ggtitle("Stages I, II and III")+ guides(fill="none")

  # Subset taqble
  merged_data_patient_analysis<-merged_data_patient_info[merged_data_patient_info$stages %in% c("Stage I","Stage II","Stage III"),]

  # Split normal and tumor samples
  tumor_samnples<-unique(data.frame(sample_id=merged_data_patient_analysis[merged_data_patient_analysis$tissue_type=="Tumor","sample_id",""],tissue_type="Tumor"))
  normal_samnples<-unique(data.frame(sample_id=merged_data_patient_analysis[merged_data_patient_analysis$tissue_type=="Normal","sample_id"],tissue_type="Normal"))

  # Merge tumor and normal and samples
  tissue_type<-rbind(tumor_samnples,normal_samnples)

  # Store normalized table
	normalized_table<-df_reads_count_all_projects[[normalized_table_names]][unique(c(selected_genes_Stage_I_gene,selected_genes_Stage_II_gene,selected_genes_Stage_III_gene)),tissue_type$sample_id]

  # Melt data.frame
  melt_normalized_table_pca<-melt(normalized_table_pca)

  # Rename colllumns
  colnames(melt_normalized_table_pca)<-c("ENSEMBL","sample_id","TPM")

  # Merge information
  melt_normalized_table_pca<-merge(melt_normalized_table_pca,tissue_type,by="sample_id")

  # calculate the pca
  pca_res_df_genes_tumor <- prcomp(normalized_table, scale. = TRUE) 
  autoplot(pca_res_df_genes_tumor, data = melt_normalized_table_pca, colour = 'tissue_type')

  
  # FindClusters_resolution
  png(filename=paste(output_dir,"Ven_Diagrams.png",sep=""), width = 28, height = 14, res=600, units = "cm")
    grid.arrange(stages_I_II_III, stages_I_II_III_unique, nrow = 1)
  dev.off()
  #############################################################################################################################
  write_tsv(selected_genes_Stage_I_data[unique_stage_I,],     paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","unique_stage_I",".tsv",sep=""))
  write_tsv(selected_genes_Stage_II_data[unique_stage_II,],   paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","unique_stage_II",".tsv",sep=""))
  write_tsv(selected_genes_Stage_III_data[unique_stage_III,], paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","unique_stage_III",".tsv",sep=""))
  #############################################################################################################################
  cat(print(paste("\n", normalization_scheme, " Number of tumor genes per stage for Stage I : ",paste(length(rownames(selected_genes_Stage_I_data)),"/",length(unique_stage_I),sep=""))),file=paste(output_dir,"outfile.txt",sep="/"),append=TRUE)
  cat(print(paste("\n", normalization_scheme, " Number of tumor genes per stage for Stage II : ",paste(length(rownames(selected_genes_Stage_II_data)),"/",length(unique_stage_II),sep=""))),file=paste(output_dir,"outfile.txt",sep="/"),append=TRUE)
  cat(print(paste("\n", normalization_scheme, " Number of tumor genes per stage for Stage III : ",paste(length(rownames(selected_genes_Stage_III_data)),"/",length(unique_stage_III),sep=""))),file=paste(output_dir,"outfile.txt",sep="/"),append=TRUE)
}
