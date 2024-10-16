####################################################################################################################################################
# Interactomes_GC3_T2.csv file has 15650 entries. The number of annotated genes with gene length geneLength_ENTREZID_ENSEMBL is 14609. Among these, 14726 are common to Interactomes_GC3_T2 and geneLength_ENTREZID_ENSEMBL and will be used to create the maps. 
# Consitency - check filters meticulously.
# FPKM, TPM  - take these as robust.
# Paramter to set the normalization_scheme
normalization_schemes<-c("tpm","fpkm","tmm","rpkm","tpm_calc")
normalization_schemes<-c("tpm","tmm")

# For each normlization normalization_scheme
for (normalization_scheme in normalization_schemes)
{
  #######################################################################################################################################
  # Set the field stages
  Interactomes_GC3_T2_merged$Stages<-"overllapping"
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
  Interactomes_GC3_T2_merged[unique_stage_I,"Stages"]<-"Stage I"
  Interactomes_GC3_T2_merged[unique_stage_II,"Stages"]<-"Stage II"
  Interactomes_GC3_T2_merged[unique_stage_III,"Stages"]<-"Stage III"
  Interactomes_GC3_T2_merged<-na.omit(Interactomes_GC3_T2_merged)

  # Visualize: Specify the comparisons you want
  my_comparisons <- list( c("Stage I", "Stage II"), c("Stage I", "Stage III"), c("Stage II", "Stage III"), c("Stage I", "overllapping"),c("Stage II", "overllapping"),c("Stage III", "overllapping"))
  
  # Add global p-value
  # Add pairwise comparisons p-value
  t2_with_outliers<-ggboxplot(Interactomes_GC3_T2_merged, x = "Stages", y = "T2",color = "Stages", palette = "jco")+  stat_compare_means(comparisons = my_comparisons, method = "t.test")+ stat_compare_means(label.y = 50)  + ggtitle(paste(TCGA_project,normalization_scheme,"T2",sep=" : "))
  connections_with_outliers<-ggboxplot(Interactomes_GC3_T2_merged, x = "Stages", y = "Conections", color = "Stages", palette = "jco")+  stat_compare_means(comparisons = my_comparisons, method = "t.test")+ stat_compare_means(label.y = 50)  + ggtitle(paste(TCGA_project,normalization_scheme,"Conections",sep=" : "))
  
  # FindClusters_resolution
  png(filename=paste(output_dir,"T2_Connection_Stat_",normalization_scheme,".png",sep=""), width = 28, height = 14, res=1200, units = "cm")
    ggarrange(t2_with_outliers, connections_with_outliers, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
  dev.off()  
}




  
