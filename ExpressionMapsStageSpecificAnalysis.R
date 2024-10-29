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

    Interactomes_GC3_T2_merged_Stage_I    <-na.omit(Interactomes_GC3_T2_merged_Stage_I[unique_stage_I,])
    Interactomes_GC3_T2_merged_Stage_II   <-na.omit(Interactomes_GC3_T2_merged_Stage_II[unique_stage_II,])
    Interactomes_GC3_T2_merged_Stage_III  <-na.omit(Interactomes_GC3_T2_merged_Stage_III[unique_stage_III,])

    # Merge Stages
    Interactomes_GC3_T2_merged_all<-na.omit(rbind(Interactomes_GC3_T2_merged_Stage_I,Interactomes_GC3_T2_merged_Stage_II, Interactomes_GC3_T2_merged_Stage_III))  
    #########################################################################################################################################
    Interactomes_GC3_T2_merged_all$AveExp
    Interactomes_GC3_T2_merged_all[,c("T2","GC3","Conections","AveExp")]
    
    pca_res_T2_Connections <- prcomp(Interactomes_GC3_T2_merged_all[,c("T2","Conections","AveExp")], scale. = TRUE) 
    pca_res_GC3_Connections <- prcomp(Interactomes_GC3_T2_merged_all[,c("GC3","Conections","AveExp")], scale. = TRUE) 
    
    pca_T2  <-autoplot(pca_res, data = Interactomes_GC3_T2_merged_all, colour = 'Stages') +  theme_bw() + ggtitle("T2, Conections, AveExp")
    pca_GC3 <-autoplot(pca_res, data = Interactomes_GC3_T2_merged_all, colour = 'Stages') +  theme_bw() + ggtitle("GC3, Conections, AveExp")

    pca_res <- prcomp(Interactomes_GC3_T2_merged_all[,c("T2","Conections","AveExp")], scale. = TRUE) 
    autoplot(pca_res, data = Interactomes_GC3_T2_merged_all, colour = 'Stages') +  theme_bw()


    
}
