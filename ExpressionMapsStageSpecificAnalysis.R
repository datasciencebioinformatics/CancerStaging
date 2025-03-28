###########################################################################################################
source("/home/felipe/Documents/github/CancerStaging/CancerStaging_SetupAllParamters.R")                   #
source("/home/felipe/Documents/github/CancerStaging/CancerStaging_LoadRPackages.R")                       #
###########################################################################################################
# Restore the object                                                                                      #
Interactomes_GC3_T2_merged <-readRDS(file = paste(output_dir,"Interactomes_GC3_T2_merged.rds",sep=""))    #
normalization_schemes      <-readRDS(file = paste(output_dir,"normalization_schemes.rds",sep=""))         #
df_reads_count_all_projects<-readRDS(file = paste(output_dir,"df_reads_count_all_projects.rds",sep=""))   #
list_of_comparisson        <-readRDS(file = paste(output_dir,"list_of_comparisson.rds",sep=""))           #
###########################################################################################################
# Data frame to store genes and stages                                                                    #
df_genes_stage<-data.frame(ENSEMBL=c(),Normalization_scheme=c())                                            #

# For each normlization normalization_scheme
for (normalization_scheme in normalization_schemes)
{     
    # genes_stages_I
    genes_stages_I    <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_I",".tsv",sep=""), sep = '\t', header = TRUE)$gene #
    genes_stages_II   <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_II",".tsv",sep=""), sep = '\t', header = TRUE)$gene #
    genes_stages_III  <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_III",".tsv",sep=""), sep = '\t', header = TRUE)$gene #
  
    # Take also expression data from the normalization scheme set by "normalization_scheme"
    expression_table_normalized<-df_reads_count_all_projects[[normalization_scheme]]

    # Stage specific genes from each stage
    expression_table_normalized_stage_I  <-expression_table_normalized[genes_stages_I,list_of_comparisson[["sample_stage_I"]]]
    expression_table_normalized_stage_II <-expression_table_normalized[genes_stages_II,list_of_comparisson[["sample_stage_II"]]]
    expression_table_normalized_stage_III<-expression_table_normalized[genes_stages_III,list_of_comparisson[["sample_stage_III"]]]

    # ENSEMBL_ids
    ENSEMBL_ids_stage_I  <-unique(intersect(rownames(expression_table_normalized_stage_I),  Interactomes_GC3_T2_merged$ENSEMBL))
    ENSEMBL_ids_stage_II <-unique(intersect(rownames(expression_table_normalized_stage_II), Interactomes_GC3_T2_merged$ENSEMBL))
    ENSEMBL_ids_stage_III<-unique(intersect(rownames(expression_table_normalized_stage_III),Interactomes_GC3_T2_merged$ENSEMBL))

    # Set AveExp to zero Interactomes_GC3_T2_merged
    Interactomes_GC3_T2_merged_Stage_I   <-Interactomes_GC3_T2_merged[ENSEMBL_ids_stage_I,]
    Interactomes_GC3_T2_merged_Stage_II  <-Interactomes_GC3_T2_merged[ENSEMBL_ids_stage_II,]
    Interactomes_GC3_T2_merged_Stage_III <-Interactomes_GC3_T2_merged[ENSEMBL_ids_stage_III,]
  
    Interactomes_GC3_T2_merged_Stage_I$AveExp<-0
    Interactomes_GC3_T2_merged_Stage_II$AveExp<-0
    Interactomes_GC3_T2_merged_Stage_III$AveExp<-0
    
    # Calculate the average expression for the epression of each g
    Interactomes_GC3_T2_merged_Stage_I[ENSEMBL_ids_stage_I,"AveExp"]<-rowMeans(expression_table_normalized_stage_I[ENSEMBL_ids_stage_I,])
    Interactomes_GC3_T2_merged_Stage_II[ENSEMBL_ids_stage_II,"AveExp"]<-rowMeans(expression_table_normalized_stage_II[ENSEMBL_ids_stage_II,])
    Interactomes_GC3_T2_merged_Stage_III[ENSEMBL_ids_stage_III,"AveExp"]<-rowMeans(expression_table_normalized_stage_III[ENSEMBL_ids_stage_III,])     

    Interactomes_GC3_T2_merged_Stage_I$Stages<-"Stage I"
    Interactomes_GC3_T2_merged_Stage_II$Stages<-"Stage II"
    Interactomes_GC3_T2_merged_Stage_III$Stages<-"Stage III"
    #########################################################################################################################################
    # Visualize: Specify the comparisons you want
    my_comparisons <- list( c("Stage I", "Stage II"), c("Stage I", "Stage III"), c("Stage II", "Stage III"))

    unique_stage_I    =intersect(setdiff(genes_stages_I, c(genes_stages_II,genes_stages_III)),genes_stages_I)
    unique_stage_II   =intersect(setdiff(genes_stages_II, c(genes_stages_I,genes_stages_III)),genes_stages_II)
    unique_stage_III  =intersect(setdiff(genes_stages_III, c(genes_stages_I,genes_stages_II)),genes_stages_III)

    # Data frame to store genes and stages                                                                    #
    df_genes_stage<-rbind(df_genes_stage,data.frame(ENSEMBL=unique(c(unique_stage_I,unique_stage_II,unique_stage_III)),Normalization_scheme=normalization_scheme))
    ###########################################################################################################

    Interactomes_GC3_T2_merged_Stage_I    <-na.omit(Interactomes_GC3_T2_merged_Stage_I[unique_stage_I,])
    Interactomes_GC3_T2_merged_Stage_II   <-na.omit(Interactomes_GC3_T2_merged_Stage_II[unique_stage_II,])
    Interactomes_GC3_T2_merged_Stage_III  <-na.omit(Interactomes_GC3_T2_merged_Stage_III[unique_stage_III,])

    # Merge Stages
    Interactomes_GC3_T2_merged_all<-na.omit(rbind(Interactomes_GC3_T2_merged_Stage_I,Interactomes_GC3_T2_merged_Stage_II, Interactomes_GC3_T2_merged_Stage_III))  
    #########################################################################################################################################   
    pca_res_T2_Connections <- prcomp(Interactomes_GC3_T2_merged_all[,c("T2","Conections","AveExp")], scale. = TRUE) 
    pca_res_GC3_Connections <- prcomp(Interactomes_GC3_T2_merged_all[,c("GC3","Conections","AveExp")], scale. = TRUE) 
    
    pca_T2  <-autoplot(pca_res_T2_Connections, data = Interactomes_GC3_T2_merged_all, colour = 'Stages') +  theme_bw() + ggtitle("T2, Conections, AveExp")
    pca_GC3 <-autoplot(pca_res_GC3_Connections, data = Interactomes_GC3_T2_merged_all, colour = 'Stages') +  theme_bw() + ggtitle("GC3, Conections, AveExp")

    # Arrange density plot
    pcas_plot<-ggarrange(pca_T2, pca_GC3, nrow = 1,ncol = 2, common.legend = TRUE, legend="bottom")  
       
    # FindClusters_resolution          
    png(filename=paste(output_dir,"countour_T2_Coonections_melt_",normalization_scheme,"_",TCGA_project,"_pcas.png",sep=""), width = 20, height = 10, res=1200, units = "cm")          
        print(annotate_figure(pcas_plot, top = text_grob(paste(normalization_scheme,TCGA_project,sep=" : "), face = "bold", size = 14)))
    dev.off()    
    #########################################################################################################################################   
    
}
# Transform values into z-score
Interactomes_GC3_T2_merged_all$T2_score  <-calculate_z(Interactomes_GC3_T2_merged_all$T2,mean(Interactomes_GC3_T2_merged_all$T2),sd(Interactomes_GC3_T2_merged_all$T2))
Interactomes_GC3_T2_merged_all$GC3_score <-calculate_z(Interactomes_GC3_T2_merged_all$GC3,mean(Interactomes_GC3_T2_merged_all$GC3),sd(Interactomes_GC3_T2_merged_all$GC3))
Interactomes_GC3_T2_merged_all$AveExp_score <-calculate_z(Interactomes_GC3_T2_merged_all$AveExp,mean(Interactomes_GC3_T2_merged_all$AveExp),sd(Interactomes_GC3_T2_merged_all$AveExp))
Interactomes_GC3_T2_merged_all$Conections_score <-calculate_z(Interactomes_GC3_T2_merged_all$Conections,mean(Interactomes_GC3_T2_merged_all$Conections),sd(Interactomes_GC3_T2_merged_all$Conections))

#########################################################################################################################################   
# Merge data.frame to analyze normalization schemes
df_genes_stage<-merge(df_genes_stage,Interactomes_GC3_T2_merged_all,by="ENSEMBL")

# Merge data.frame to analyze normalization schemes
pca_res_df_genes_stage <- prcomp(df_genes_stage[,c("T2","Conections","AveExp")], scale. = TRUE) 

# Merge data.frame to analyze normalization schemes
pca_res_df_genes_stage  <-autoplot(pca_res_df_genes_stage, data = df_genes_stage, colour = 'Normalization_scheme') +  theme_bw() + ggtitle("T2, Conections, AveExp")

# FindClusters_resolution          
png(filename=paste(output_dir,"countour_T2_Coonections_melt_",normalization_scheme,"_",TCGA_project,"_pcas_tmm_tpm.png",sep=""), width = 15, height = 10, res=600, units = "cm")          
    pca_res_df_genes_stage
dev.off()    


#########################################################################################################################################   
# Merge data.frame to analyze normalization schemes
pca_res_df_genes_stage <- prcomp(df_genes_stage[,c("T2_score","Conections_score","AveExp_score")], scale. = TRUE) 

# Merge data.frame to analyze normalization schemes
pca_res_df_genes_stage  <-autoplot(pca_res_df_genes_stage, data = df_genes_stage, colour = 'Normalization_scheme') +  theme_bw() + ggtitle("Z-score T2, Conections, AveExp")

# FindClusters_resolution          
png(filename=paste(output_dir,"countour_T2_Coonections_melt_",normalization_scheme,"_",TCGA_project,"_pcas_tmm_tpm_zscore.png",sep=""), width = 15, height = 10, res=600, units = "cm")          
    pca_res_df_genes_stage
dev.off()    
