list_stage_specific_genes[[paste(normalization_scheme,"_",substring(Stage_i,8,14),sep="")]]<-log2change_Stage_i

# for each  normalization scheme
for (normalization_scheme in normalization_schemes)
{  
  list_stage_specific_genes[[paste(normalization_scheme,"_","stage_I",sep="")]]
  list_stage_specific_genes[[paste(normalization_scheme,"_","stage_II",sep="")]]
  list_stage_specific_genes[[paste(normalization_scheme,"_","stage_III",sep="")]]
  #######################################################################################################################################
  # The stage specific genes in list                                                                                                    #
  # paste(normalization_scheme,"_",substring(Stage_i,8,14),sep="")                                                                      #
  # list_stage_specific_genes                                                                                                           #
  #######################################################################################################################################
  # Path to files of selected_genes                                                                                                     # 
  selected_genes_Stage_I_file       <-paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_Stage_","sample_stage_I",".tsv",sep="")       #
  selected_genes_Stage_II_file      <-paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_Stage_","sample_stage_II",".tsv",sep="")      #
  selected_genes_Stage_III_file     <-paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_Stage_","sample_stage_III",".tsv",sep="")     #
  #######################################################################################################################################
  # Load data                                                                                                                           #
  selected_genes_Stage_I_data       <-read.table(file = selected_genes_Stage_I_file, sep = '\t', header = TRUE,fill=TRUE)               #
  selected_genes_Stage_II_data      <-read.table(file = selected_genes_Stage_II_file, sep = '\t', header = TRUE,fill=TRUE)              #
  selected_genes_Stage_III_data     <-read.table(file = selected_genes_Stage_III_file, sep = '\t', header = TRUE,fill=TRUE)             #
                                                                                                                                        #
  # Set rownames                                                                                                                        #
  rownames(selected_genes_Stage_I_data)<-selected_genes_Stage_I_data$gene                                                               #
  rownames(selected_genes_Stage_II_data)<-selected_genes_Stage_II_data$gene                                                             #
  rownames(selected_genes_Stage_III_data)<-selected_genes_Stage_III_data$gene                                                           #
  #######################################################################################################################################
  unique_stage_I  =intersect(setdiff(selected_genes_Stage_I_data$gene, c(selected_genes_Stage_II_data$gene,selected_genes_Stage_III_data$gene)),selected_genes_Stage_I_data$gene)
  unique_stage_II =intersect(setdiff(selected_genes_Stage_II_data$gene, c(selected_genes_Stage_I_data$gene,selected_genes_Stage_III_data$gene)),selected_genes_Stage_II_data$gene)
  unique_stage_III=intersect(setdiff(selected_genes_Stage_III_data$gene, c(selected_genes_Stage_I_data$gene,selected_genes_Stage_II_data$gene)),selected_genes_Stage_III_data$gene)
  #######################################################################################################################################
  # Select stages 
  stages_I_II_III_unique<-ggVennDiagram(list(Stage_I=unique_stage_I,Stage_II=unique_stage_II,Stage_III=unique_stage_III), label_alpha = 0.9,set_color = c("grey50","grey50","grey50")) +  scale_fill_gradient(low = "white", high = "white") + theme_bw() + ggtitle("Stages I, II and III")+ guides(fill="none")
  stages_I_II_III<-ggVennDiagram(list(Stage_I=selected_genes_Stage_I_data$gene,Stage_II=selected_genes_Stage_II_data$gene,Stage_III=selected_genes_Stage_III_data$gene), label_alpha = 0.9,set_color = c("grey50","grey50","grey50")) + scale_fill_gradient(low = "white", high = "white") + theme_bw() + ggtitle("Stages I, II and III")+ guides(fill="none")
  
  # FindClusters_resolution
  png(filename=paste(output_dir,"Ven_Diagrams.png",sep=""), width = 28, height = 14, res=600, units = "cm")
    grid.arrange(stages_I_II_III, stages_I_II_III_unique, nrow = 1)
  dev.off()
  #############################################################################################################################
  write_tsv(selected_genes_Stage_I_data[unique_stage_I,], paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_unique_stage_I",".tsv",sep=""))
  write_tsv(selected_genes_Stage_II_data[unique_stage_II,], paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_unique_stage_II",".tsv",sep=""))
  write_tsv(selected_genes_Stage_III_data[unique_stage_III,], paste(output_dir,"DE_GenesPerStageMeansFromPairedUp_unique_stage_III",".tsv",sep=""))
  #############################################################################################################################
  cat(print(paste("\nNumber of tumor genes per stage for Stage I : ",paste(length(rownames(selected_genes_Stage_I_data)),"/",length(unique_stage_I),sep=""))),file=paste(output_dir,"outfile.txt",sep="/"),append=TRUE)
  cat(print(paste("\nNumber of tumor genes per stage for Stage II : ",paste(length(rownames(selected_genes_Stage_II_data)),"/",length(unique_stage_II),sep=""))),file=paste(output_dir,"outfile.txt",sep="/"),append=TRUE)
  cat(print(paste("\nNumber of tumor genes per stage for Stage III : ",paste(length(rownames(selected_genes_Stage_III_data)),"/",length(unique_stage_III),sep=""))),file=paste(output_dir,"outfile.txt",sep="/"),append=TRUE)
}
