####################################################################################################################################################
# Interactomes_GC3_T2.csv file has 15650 entries. The number of annotated genes with gene length geneLength_ENTREZID_ENSEMBL is 14609. Among these, 14726 are common to Interactomes_GC3_T2 and geneLength_ENTREZID_ENSEMBL and will be used to create the maps. 
# Consitency - check filters meticulously.
# FPKM, TPM  - take these as robust.
# Paramter to set the normalization_scheme
normalization_schemes<-c("tpm","fpkm","tmm","rpkm","tpm_calc")

# For each normlization normalization_scheme
for (normalization_scheme in normalization_schemes)
{
  #######################################################################################################################################
  # Path to files of selected_genes                                                                                                             # 
  # genes_stages_I
  selected_genes_Stage_I_data    <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_I",".tsv",sep=""), sep = '\t', header = TRUE) #
  selected_genes_Stage_II_data   <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_II",".tsv",sep=""), sep = '\t', header = TRUE) #
  selected_genes_Stage_III_data  <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_III",".tsv",sep=""), sep = '\t', header = TRUE) #

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
  stages_I_II_III<-ggVennDiagram(list(Stage_I=selected_genes_Stage_I_data,Stage_II=selected_genes_Stage_II_data,Stage_III=selected_genes_Stage_III_data), label_alpha = 0.9,set_color = c("grey50","grey50","grey50")) + scale_fill_gradient(low = "white", high = "white") + theme_bw() + ggtitle("Stages I, II and III")+ guides(fill="none")
  
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
