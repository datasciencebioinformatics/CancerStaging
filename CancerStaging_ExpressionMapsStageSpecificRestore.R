###########################################################################################################
source("/home/felipe/Documents/github/CancerStaging/CancerStaging_SetupAllParamters.R")                   #
source("/home/felipe/Documents/github/CancerStaging/CancerStaging_LoadRPackages.R")                       #
###########################################################################################################
# Restore the object                                                                                      #
Interactomes_GC3_T2_merged <-readRDS(file = paste(output_dir,"Interactomes_GC3_T2_merged.rds",sep=""))    #
#normalization_schemes      <-readRDS(file = paste(output_dir,"normalization_schemes.rds",sep=""))         #
df_reads_count_all_projects<-readRDS(file = paste(output_dir,"df_reads_count_all_projects.rds",sep=""))   #
list_of_comparisson        <-readRDS(file = paste(output_dir,"list_of_comparisson.rds",sep=""))           #
###########################################################################################################
normalization_schemes<-c("tpm")

# For each normlization normalization_scheme
for (normalization_scheme in normalization_schemes)
{     
    # genes_stages_I
    genes_stages_I    <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_I",".tsv",sep=""), sep = '\t', header = TRUE)$gene #
    genes_stages_II   <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_II",".tsv",sep=""), sep = '\t', header = TRUE)$gene #
    genes_stages_III  <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_III",".tsv",sep=""), sep = '\t', header = TRUE)$gene #

    unique_stage_I  =intersect(setdiff(genes_stages_I, c(genes_stages_II,genes_stages_III)),genes_stages_I)
    unique_stage_II =intersect(setdiff(genes_stages_II, c(genes_stages_I,genes_stages_III)),genes_stages_II)
    unique_stage_III=intersect(setdiff(genes_stages_III, c(genes_stages_I,genes_stages_II)),genes_stages_III)
  
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
    
    ENSEMBL_ids_stage_I<-unique_stage_I
    ENSEMBL_ids_stage_II<-unique_stage_II
    ENSEMBL_ids_stage_III<-unique_stage_III

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

    # Calculate the average expression for the epression of each g
    Interactomes_GC3_T2_merged_Stage_I[ENSEMBL_ids_stage_I,"AveExp"]<-rowMeans(expression_table_normalized_stage_I[unique_stage_I,])
    Interactomes_GC3_T2_merged_Stage_II[ENSEMBL_ids_stage_II,"AveExp"]<-rowMeans(expression_table_normalized_stage_II[unique_stage_II,])
    Interactomes_GC3_T2_merged_Stage_III[ENSEMBL_ids_stage_III,"AveExp"]<-rowMeans(expression_table_normalized_stage_III[unique_stage_III,])     
    

    Interactomes_GC3_T2_merged_Stage_I$Stages<-"Stage I"
    Interactomes_GC3_T2_merged_Stage_II$Stages<-"Stage II"
    Interactomes_GC3_T2_merged_Stage_III$Stages<-"Stage III"
    #########################################################################################################################################
    # Visualize: Specify the comparisons you want
    my_comparisons <- list( c("Stage I", "Stage II"), c("Stage I", "Stage III"), c("Stage II", "Stage III"))


    # Data frame to store genes and stages                                                                    #
    df_genes_stage<-data.frame(Genes=c(),Stage=c())      

    Interactomes_GC3_T2_merged_Stage_I<-na.omit(Interactomes_GC3_T2_merged_Stage_I[unique_stage_I,])
    Interactomes_GC3_T2_merged_Stage_II<-na.omit(Interactomes_GC3_T2_merged_Stage_II[unique_stage_II,])
    Interactomes_GC3_T2_merged_Stage_III<-na.omit(Interactomes_GC3_T2_merged_Stage_III[unique_stage_III,])

    #Interactomes_GC3_T2_merged_Stage_I<-na.omit(Interactomes_GC3_T2_merged_Stage_I[ENSEMBL_ids_stage_I,])
    #Interactomes_GC3_T2_merged_Stage_II<-na.omit(Interactomes_GC3_T2_merged_Stage_II[ENSEMBL_ids_stage_II,])
    #Interactomes_GC3_T2_merged_Stage_III<-na.omit(Interactomes_GC3_T2_merged_Stage_III[ENSEMBL_ids_stage_III,])
    
    
    # Merge Stages
    Interactomes_GC3_T2_merged_all<-na.omit(rbind(Interactomes_GC3_T2_merged_Stage_I,Interactomes_GC3_T2_merged_Stage_II, Interactomes_GC3_T2_merged_Stage_III))  
    ###########################################################################################################    
    # Three countour plots will be created
    # One with the average expression
    # Second with the expression per patient
    # Third with the z-score 
    # Conections, T2, AvgExpression
    # Combine AvgExpression, Conections, T2
    # harmonic mean
    # z-core : Composite Scores
    # Implementing the Z score formula in R is quite straightforward. 
    # To reuse code, we will create a function called calculate_z using the mean and sd base functions to calculate Z. 
    # sd calculates the standard deviation in R.
    # weighted average
    # Z-score for AveExp_expression
    # Z-score for AveExp_expression  
    # Z-score for AveExp_expression
    # Z-score for AveExp_expression       
    m1<-ggplot(Interactomes_GC3_T2_merged_Stage_I, aes(Conections, T2, z = AveExp))   + geom_density_2d_filled() + theme_bw() + ggtitle("Stage I")      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))      + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + geom_hline(yintercept=30, colour="yellow") + theme(legend.position = "none") + ggtitle("F")
    m2<-ggplot(Interactomes_GC3_T2_merged_Stage_II, aes(Conections, T2, z = AveExp))   + geom_density_2d_filled() + theme_bw() + ggtitle("Stage II")      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))    + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + geom_hline(yintercept=30, colour="yellow") + theme(legend.position = "none") + ggtitle("G")  
    m3<-ggplot(Interactomes_GC3_T2_merged_Stage_III, aes(Conections, T2, z = AveExp))   + geom_density_2d_filled() + theme_bw() + ggtitle("Stage III")      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))  + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + geom_hline(yintercept=30, colour="yellow") + theme(legend.position = "none") + ggtitle("H")  
    
    # Arrange density plot
    density_plot<-ggarrange(m1, m2, m3, nrow = 1,ncol = 3, common.legend = TRUE, legend="none")
    
    m7<-ggplot(Interactomes_GC3_T2_merged_Stage_I, aes(Conections, T2, z = AveExp))    + geom_point(aes(colour=AveExp),size=3)   +  theme_bw() + ggtitle("Stage I")+  xlim(0, 50)     + ylim(10,40) + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(0, 40, by = 10), limits = c(10, 40))      + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + scale_colour_gradient(low = "red", high = "yellow") + geom_hline(yintercept=30, colour="yellow")   + ggtitle("A")
    m8<-ggplot(Interactomes_GC3_T2_merged_Stage_II, aes(Conections, T2, z = AveExp))   + geom_point(aes(colour=AveExp),size=3)  +  theme_bw() + ggtitle("Stage II")+ xlim(0, 50)       + ylim(10,40) + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(0, 40, by = 10), limits = c(10, 40))      + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  + scale_colour_gradient(low = "red", high = "yellow") + geom_hline(yintercept=30, colour="yellow")   + ggtitle("B")
    m9<-ggplot(Interactomes_GC3_T2_merged_Stage_III, aes(Conections, T2, z = AveExp))  + geom_point(aes(colour=AveExp),size=3) + theme_bw() + ggtitle("Stage III")+ xlim(0, 50)         + ylim(10,40) +  guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE)) + ylim(10,40)  + xlim(0, 50)     + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(0, 40, by = 10), limits = c(10, 40))   + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + scale_colour_gradient(low = "red", high = "yellow") + geom_hline(yintercept=30, colour="yellow") + ggtitle("C")

    # Arrange density plot
    dotplot_plot<-ggarrange(m7, m8, m9, nrow = 1,ncol = 3, common.legend = TRUE, legend="bottom")        
  
    m10<-ggplot(Interactomes_GC3_T2_merged_Stage_I, aes(Conections, T2, z = AveExp))   +  theme_bw() + ggtitle("Stage I")   + geom_density_2d(bin=10)      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE)) + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(0, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))   + scale_colour_gradient(low = "red", high = "yellow")   + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  + theme(legend.position = "none") + geom_hline(yintercept=30, colour="yellow") + ggtitle("D")
    m11<-ggplot(Interactomes_GC3_T2_merged_Stage_II, aes(Conections, T2, z = AveExp))  +  theme_bw() + ggtitle("Stage II")   + geom_density_2d(bin=10)      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(0, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))   + scale_colour_gradient(low = "red", high = "yellow")   + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) + theme(legend.position = "none") + geom_hline(yintercept=30, colour="yellow") + ggtitle("E")
    m12<-ggplot(Interactomes_GC3_T2_merged_Stage_III, aes(Conections, T2, z = AveExp))  +  theme_bw() + ggtitle("Stage III")   + geom_density_2d(bin=10)      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(0, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40)) + scale_colour_gradient(low = "red", high = "yellow") + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + theme(legend.position = "none") + geom_hline(yintercept=30, colour="yellow") + ggtitle("F")

    # Arrange density plot
    countour_plot<-ggarrange(m10, m11, m12, nrow = 1,ncol = 3, common.legend = TRUE, legend="bottom")            
      
    # FindClusters_resolution          
    # Supplmental_Figure
    png(filename=paste(output_dir,"countour_T2_Coonections_melt_",normalization_scheme,"_",TCGA_project,"_Stage_all_T2_perPatient_unique_paper.png",sep=""), width = 30, height = 30, res=600, units = "cm")  
            ggarrange(dotplot_plot, countour_plot, density_plot,  nrow = 3,ncol = 1, common.legend = TRUE, legend="bottom")
            #print(annotate_figure(plot, top = text_grob(TCGA_project, face = "bold", size = 14)))
    dev.off()



    ###########################################################################################################
    # Take ensembl ids
    ENSEMBL_IDs<-intersect(rownames(Interactomes_GC3_T2_merged_all),rownames(expression_table_normalized))

    # Copy the epxression table
    expression_table_normalized_cp<-expression_table_normalized[ENSEMBL_IDs,]

    # Add ENSEMBL to the table
    expression_table_normalized_cp$ENSEMBL<-rownames(expression_table_normalized_cp)

    # Melt table
    expression_table_normalized_melt<-melt(expression_table_normalized_cp,ids="ENSEMBL")

    # Merge tables
    expression_table_normalized_melt<-merge(expression_table_normalized_melt,Interactomes_GC3_T2_merged_all,by="ENSEMBL")

    # Rename collumns
    colnames(expression_table_normalized_melt)[3]<-normalization_scheme

    # Take values from each stage
    Interactomes_GC3_T2_merged_Stage_I_patient<-expression_table_normalized_melt[expression_table_normalized_melt$Stages == "Stage I",]
    Interactomes_GC3_T2_merged_Stage_II_patient<-expression_table_normalized_melt[expression_table_normalized_melt$Stages == "Stage II",]
    Interactomes_GC3_T2_merged_Stage_III_patient<-expression_table_normalized_melt[expression_table_normalized_melt$Stages == "Stage III",]
    #########################################################################################################################################

    ###########################################################################################################    
    # Three countour plots will be created
    # One with the average expression
    # Second with the expression per patient
    # Third with the z-score 
    # Conections, T2, AvgExpression
    # Combine AvgExpression, Conections, T2
    # harmonic mean
    # z-core : Composite Scores
    # Implementing the Z score formula in R is quite straightforward. 
    # To reuse code, we will create a function called calculate_z using the mean and sd base functions to calculate Z. 
    # sd calculates the standard deviation in R.
    # weighted average
    # Z-score for AveExp_expression
    # Z-score for AveExp_expression  
    # Z-score for AveExp_expression
    # Z-score for AveExp_expression  
    
    m1<-ggplot(Interactomes_GC3_T2_merged_Stage_I_patient, aes(Conections, T2, z = AveExp))   + geom_density_2d_filled() + theme_bw() + ggtitle("Stage I")      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))      + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + geom_hline(yintercept=30, colour="yellow") + theme(legend.position = "none") + ggtitle("F") + geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_I_patient$T2), colour="red")
    m2<-ggplot(Interactomes_GC3_T2_merged_Stage_II_patient, aes(Conections, T2, z = AveExp))   + geom_density_2d_filled() + theme_bw() + ggtitle("Stage II")      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))    + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + geom_hline(yintercept=30, colour="yellow") + theme(legend.position = "none") + ggtitle("G") + geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_II_patient$T2), colour="red")
    m3<-ggplot(Interactomes_GC3_T2_merged_Stage_III_patient, aes(Conections, T2, z = AveExp))   + geom_density_2d_filled() + theme_bw() + ggtitle("Stage III")      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))  + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + geom_hline(yintercept=30, colour="yellow") + theme(legend.position = "none") + ggtitle("H") + geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_III_patient$T2), colour="red")
    
    # Arrange density plot
    density_plot<-ggarrange(m1, m2, m3, nrow = 1,ncol = 3, common.legend = TRUE, legend="none")
    
    m7<-ggplot(Interactomes_GC3_T2_merged_Stage_I_patient, aes(Conections, T2, z = AveExp))    + geom_point(aes(colour=AveExp),size=3)   +  theme_bw() + ggtitle("Stage I")+  xlim(0, 50)     + ylim(10,40) + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(0, 40, by = 10), limits = c(10, 40))      + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + scale_colour_gradient(low = "red", high = "yellow") + geom_hline(yintercept=30, colour="yellow")   + ggtitle("A") + geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_I_patient$T2), colour="red")
    m8<-ggplot(Interactomes_GC3_T2_merged_Stage_II_patient, aes(Conections, T2, z = AveExp))   + geom_point(aes(colour=AveExp),size=3)  +  theme_bw() + ggtitle("Stage II")+ xlim(0, 50)       + ylim(10,40) + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(0, 40, by = 10), limits = c(10, 40))      + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  + scale_colour_gradient(low = "red", high = "yellow") + geom_hline(yintercept=30, colour="yellow")   + ggtitle("B") + geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_II_patient$T2), colour="red")
    m9<-ggplot(Interactomes_GC3_T2_merged_Stage_III_patient, aes(Conections, T2, z = AveExp))  + geom_point(aes(colour=AveExp),size=3) + theme_bw() + ggtitle("Stage III")+ xlim(0, 50)         + ylim(10,40) +  guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE)) + ylim(10,40)  + xlim(0, 50)     + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(0, 40, by = 10), limits = c(10, 40))   + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + scale_colour_gradient(low = "red", high = "yellow") + geom_hline(yintercept=30, colour="yellow") + ggtitle("C") + geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_III_patient$T2), colour="red")

    # Arrange density plot
    dotplot_plot<-ggarrange(m7, m8, m9, nrow = 1,ncol = 3, common.legend = TRUE, legend="bottom")        
  
    m10<-ggplot(Interactomes_GC3_T2_merged_Stage_I_patient, aes(Conections, T2, z = AveExp))   +  theme_bw() + ggtitle("Stage I")   + geom_density_2d(bin=10)      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE)) + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(0, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))   + scale_colour_gradient(low = "red", high = "yellow")   + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  + theme(legend.position = "none") + geom_hline(yintercept=30, colour="yellow") + ggtitle("D") + geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_I_patient$T2), colour="red")
    m11<-ggplot(Interactomes_GC3_T2_merged_Stage_II_patient, aes(Conections, T2, z = AveExp))  +  theme_bw() + ggtitle("Stage II")   + geom_density_2d(bin=10)      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(0, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))   + scale_colour_gradient(low = "red", high = "yellow")   + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) + theme(legend.position = "none") + geom_hline(yintercept=30, colour="yellow") + ggtitle("E") + geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_II_patient$T2), colour="red")
    m12<-ggplot(Interactomes_GC3_T2_merged_Stage_III_patient, aes(Conections, T2, z = AveExp))  +  theme_bw() + ggtitle("Stage III")   + geom_density_2d(bin=10)      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(0, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40)) + scale_colour_gradient(low = "red", high = "yellow") + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + theme(legend.position = "none") + geom_hline(yintercept=30, colour="yellow") + ggtitle("F") + geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_III_patient$T2), colour="red")

    # Arrange density plot
    countour_plot<-ggarrange(m10, m11, m12, nrow = 1,ncol = 3, common.legend = TRUE, legend="bottom")            
      
    # FindClusters_resolution
    # Supplmental_Figure
    png(filename=paste(output_dir,"countour_T2_Coonections_melt_",normalization_scheme,"_",TCGA_project,"_Stage_all_T2_perPatient_patient_paper.png",sep=""), width = 30, height = 30, res=600, units = "cm")  
            ggarrange(dotplot_plot, countour_plot, density_plot,  nrow = 3,ncol = 1, common.legend = TRUE, legend="bottom")
            #print(annotate_figure(plot, top = text_grob(TCGA_project, face = "bold", size = 14)))
    dev.off()

    m1<-ggplot(Interactomes_GC3_T2_merged_Stage_I_patient, aes(Conections, T2, z = AveExp))   + geom_density_2d_filled() + theme_bw() + ggtitle("Stage I")      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))      + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + geom_hline(yintercept=30, colour="yellow") + theme(legend.position = "none") + ggtitle("A") + geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_I_patient$T2), colour="red")
    m2<-ggplot(Interactomes_GC3_T2_merged_Stage_II_patient, aes(Conections, T2, z = AveExp))   + geom_density_2d_filled() + theme_bw() + ggtitle("Stage II")      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))    + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + geom_hline(yintercept=30, colour="yellow") + theme(legend.position = "none") + ggtitle("B") + geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_II_patient$T2), colour="red")
    m3<-ggplot(Interactomes_GC3_T2_merged_Stage_III_patient, aes(Conections, T2, z = AveExp))   + geom_density_2d_filled() + theme_bw() + ggtitle("Stage III")      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))  + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + geom_hline(yintercept=30, colour="yellow") + theme(legend.position = "none") + ggtitle("C") + geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_III_patient$T2), colour="red")
    
    # Arrange density plot
    density_plot<-ggarrange(m1, m2, m3, nrow = 1,ncol = 3, common.legend = TRUE, legend="none")
    

    # FindClusters_resolution          
    png(filename=paste(output_dir,"countour_T2_Coonections_melt_",normalization_scheme,"_",TCGA_project,"_Stage_all_T2_perPatient_patient_paper.png",sep=""), width = 30, height = 10, res=600, units = "cm")  
            ggarrange(density_plot,  nrow = 1,ncol = 1, common.legend = TRUE, legend="bottom")
            #print(annotate_figure(plot, top = text_grob(TCGA_project, face = "bold", size = 14)))
    dev.off()

}


































































# For each normlization normalization_scheme
for (normalization_scheme in normalization_schemes)
{     
    # genes_stages_I
    genes_stages_I    <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_I",".tsv",sep=""), sep = '\t', header = TRUE)$gene #
    genes_stages_II   <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_II",".tsv",sep=""), sep = '\t', header = TRUE)$gene #
    genes_stages_III  <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_III",".tsv",sep=""), sep = '\t', header = TRUE)$gene #

    unique_stage_I  =intersect(setdiff(genes_stages_I, c(genes_stages_II,genes_stages_III)),genes_stages_I)
    unique_stage_II =intersect(setdiff(genes_stages_II, c(genes_stages_I,genes_stages_III)),genes_stages_II)
    unique_stage_III=intersect(setdiff(genes_stages_III, c(genes_stages_I,genes_stages_II)),genes_stages_III)
  
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
    
    ENSEMBL_ids_stage_I<-unique_stage_I
    ENSEMBL_ids_stage_II<-unique_stage_II
    ENSEMBL_ids_stage_III<-unique_stage_III

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

    # Calculate the average expression for the epression of each g
    Interactomes_GC3_T2_merged_Stage_I[ENSEMBL_ids_stage_I,"AveExp"]<-rowMeans(expression_table_normalized_stage_I[unique_stage_I,])
    Interactomes_GC3_T2_merged_Stage_II[ENSEMBL_ids_stage_II,"AveExp"]<-rowMeans(expression_table_normalized_stage_II[unique_stage_II,])
    Interactomes_GC3_T2_merged_Stage_III[ENSEMBL_ids_stage_III,"AveExp"]<-rowMeans(expression_table_normalized_stage_III[unique_stage_III,])     
    

    Interactomes_GC3_T2_merged_Stage_I$Stages<-"Stage I"
    Interactomes_GC3_T2_merged_Stage_II$Stages<-"Stage II"
    Interactomes_GC3_T2_merged_Stage_III$Stages<-"Stage III"
    #########################################################################################################################################
    # Visualize: Specify the comparisons you want
    my_comparisons <- list( c("Stage I", "Stage II"), c("Stage I", "Stage III"), c("Stage II", "Stage III"))


    # Data frame to store genes and stages                                                                    #
    df_genes_stage<-data.frame(Genes=c(),Stage=c())      

    Interactomes_GC3_T2_merged_Stage_I<-na.omit(Interactomes_GC3_T2_merged_Stage_I[unique_stage_I,])
    Interactomes_GC3_T2_merged_Stage_II<-na.omit(Interactomes_GC3_T2_merged_Stage_II[unique_stage_II,])
    Interactomes_GC3_T2_merged_Stage_III<-na.omit(Interactomes_GC3_T2_merged_Stage_III[unique_stage_III,])

    #Interactomes_GC3_T2_merged_Stage_I<-na.omit(Interactomes_GC3_T2_merged_Stage_I[ENSEMBL_ids_stage_I,])
    #Interactomes_GC3_T2_merged_Stage_II<-na.omit(Interactomes_GC3_T2_merged_Stage_II[ENSEMBL_ids_stage_II,])
    #Interactomes_GC3_T2_merged_Stage_III<-na.omit(Interactomes_GC3_T2_merged_Stage_III[ENSEMBL_ids_stage_III,])
    
    
    # Merge Stages
    Interactomes_GC3_T2_merged_all<-na.omit(rbind(Interactomes_GC3_T2_merged_Stage_I,Interactomes_GC3_T2_merged_Stage_II, Interactomes_GC3_T2_merged_Stage_III))  
    ###########################################################################################################
    # Take ensembl ids
    ENSEMBL_IDs<-intersect(rownames(Interactomes_GC3_T2_merged_all),rownames(expression_table_normalized))

    # Copy the epxression table
    expression_table_normalized_cp<-expression_table_normalized[ENSEMBL_IDs,]

    # Add ENSEMBL to the table
    expression_table_normalized_cp$ENSEMBL<-rownames(expression_table_normalized_cp)

    # Melt table
    expression_table_normalized_melt<-melt(expression_table_normalized_cp,ids="ENSEMBL")

    # Merge tables
    expression_table_normalized_melt<-merge(expression_table_normalized_melt,Interactomes_GC3_T2_merged_all,by="ENSEMBL")

    # Rename collumns
    colnames(expression_table_normalized_melt)[3]<-normalization_scheme

    # Take values from each stage
    Interactomes_GC3_T2_merged_Stage_I_patient<-expression_table_normalized_melt[expression_table_normalized_melt$Stages == "Stage I",]
    Interactomes_GC3_T2_merged_Stage_II_patient<-expression_table_normalized_melt[expression_table_normalized_melt$Stages == "Stage II",]
    Interactomes_GC3_T2_merged_Stage_III_patient<-expression_table_normalized_melt[expression_table_normalized_melt$Stages == "Stage III",]
    #########################################################################################################################################
    
    # Three countour plots will be created
    # One with the average expression
    # Second with the expression per patient
    # Third with the z-score 
    # Conections, T2, AvgExpression
    # Combine AvgExpression, Conections, T2
    # harmonic mean
    # z-core : Composite Scores
    # Implementing the Z score formula in R is quite straightforward. 
    # To reuse code, we will create a function called calculate_z using the mean and sd base functions to calculate Z. 
    # sd calculates the standard deviation in R.
    # weighted average
    # Z-score for AveExp_expression
    # Z-score for AveExp_expression  
    # Z-score for AveExp_expression
    # Z-score for AveExp_expression       
    m1<-ggplot(Interactomes_GC3_T2_merged_Stage_I_patient, aes(Conections, T2, z = AveExp))   + geom_density_2d_filled() + theme_bw() + ggtitle("Stage I")      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))      + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + geom_hline(yintercept=30, colour="yellow") + theme(legend.position = "none") + ggtitle("F")
    m2<-ggplot(Interactomes_GC3_T2_merged_Stage_II_patient, aes(Conections, T2, z = AveExp))   + geom_density_2d_filled() + theme_bw() + ggtitle("Stage II")      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))    + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + geom_hline(yintercept=30, colour="yellow") + theme(legend.position = "none") + ggtitle("G")  
    m3<-ggplot(Interactomes_GC3_T2_merged_Stage_III_patient, aes(Conections, T2, z = AveExp))   + geom_density_2d_filled() + theme_bw() + ggtitle("Stage III")      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))  + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + geom_hline(yintercept=30, colour="yellow") + theme(legend.position = "none") + ggtitle("H")  
    
    # Arrange density plot
    density_plot<-ggarrange(m1, m2, m3, nrow = 1,ncol = 3, common.legend = TRUE, legend="none")
    
    m7<-ggplot(Interactomes_GC3_T2_merged_Stage_I_patient, aes(Conections, T2, z = AveExp))    + geom_point(aes(colour=AveExp),size=3)   +  theme_bw() + ggtitle("Stage I")+  xlim(0, 50)     + ylim(10,40) + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(0, 40, by = 10), limits = c(10, 40))      + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + scale_colour_gradient(low = "red", high = "yellow") + geom_hline(yintercept=30, colour="yellow")   + ggtitle("A")
    m8<-ggplot(Interactomes_GC3_T2_merged_Stage_II_patient, aes(Conections, T2, z = AveExp))   + geom_point(aes(colour=AveExp),size=3)  +  theme_bw() + ggtitle("Stage II")+ xlim(0, 50)       + ylim(10,40) + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(0, 40, by = 10), limits = c(10, 40))      + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  + scale_colour_gradient(low = "red", high = "yellow") + geom_hline(yintercept=30, colour="yellow")   + ggtitle("B")
    m9<-ggplot(Interactomes_GC3_T2_merged_Stage_III_patient, aes(Conections, T2, z = AveExp))  + geom_point(aes(colour=AveExp),size=3) + theme_bw() + ggtitle("Stage III")+ xlim(0, 50)         + ylim(10,40) +  guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE)) + ylim(10,40)  + xlim(0, 50)     + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(0, 40, by = 10), limits = c(10, 40))   + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + scale_colour_gradient(low = "red", high = "yellow") + geom_hline(yintercept=30, colour="yellow") + ggtitle("C")

    # Arrange density plot
    dotplot_plot<-ggarrange(m7, m8, m9, nrow = 1,ncol = 3, common.legend = TRUE, legend="bottom")        
  
    m10<-ggplot(Interactomes_GC3_T2_merged_Stage_I_patient, aes(Conections, T2, z = AveExp))   +  theme_bw() + ggtitle("Stage I")   + geom_density_2d(bin=10)      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE)) + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(0, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))   + scale_colour_gradient(low = "red", high = "yellow")   + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  + theme(legend.position = "none") + geom_hline(yintercept=30, colour="yellow") + ggtitle("D")
    m11<-ggplot(Interactomes_GC3_T2_merged_Stage_II_patient, aes(Conections, T2, z = AveExp))  +  theme_bw() + ggtitle("Stage II")   + geom_density_2d(bin=10)      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(0, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))   + scale_colour_gradient(low = "red", high = "yellow")   + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) + theme(legend.position = "none") + geom_hline(yintercept=30, colour="yellow") + ggtitle("E")
    m12<-ggplot(Interactomes_GC3_T2_merged_Stage_III_patient, aes(Conections, T2, z = AveExp))  +  theme_bw() + ggtitle("Stage III")   + geom_density_2d(bin=10)      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(0, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40)) + scale_colour_gradient(low = "red", high = "yellow") + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + theme(legend.position = "none") + geom_hline(yintercept=30, colour="yellow") + ggtitle("F")

    # Arrange density plot
    countour_plot<-ggarrange(m10, m11, m12, nrow = 1,ncol = 3, common.legend = TRUE, legend="bottom")            
      
    # FindClusters_resolution          
    png(filename=paste(output_dir,"countour_T2_Coonections_melt_",normalization_scheme,"_",TCGA_project,"_Stage_all_T2_perPatient_unique_paper.png",sep=""), width = 30, height = 30, res=600, units = "cm")  
            ggarrange(dotplot_plot, countour_plot, density_plot,  nrow = 3,ncol = 1, common.legend = TRUE, legend="bottom")
            #print(annotate_figure(plot, top = text_grob(TCGA_project, face = "bold", size = 14)))
    dev.off()

}











    m4 <- ggplot(Interactomes_GC3_T2_merged_all, aes(x=Stages, y=T2)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()            + ggtitle(paste("T2", normalization_scheme,sep=" : "))                   +  stat_compare_means(comparisons = my_comparisons, method = "t.test")  + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))
    m5 <- ggplot(Interactomes_GC3_T2_merged_all, aes(x=Stages, y=Conections)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()    + ggtitle(paste("Connectivity", normalization_scheme,sep=" : ")) +  stat_compare_means(comparisons = my_comparisons, method = "t.test")          + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) 
    m6 <- ggplot(Interactomes_GC3_T2_merged_all, aes(x=Stages, y=AveExp)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()        + ggtitle(paste("Expr.", normalization_scheme,sep=" : "))          +  stat_compare_means(comparisons = my_comparisons, method = "t.test")        + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))
    m7 <- ggplot(Interactomes_GC3_T2_merged_all, aes(x=Stages, y=GC3)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()           + ggtitle(paste("GC3", normalization_scheme,sep=" : "))          +  stat_compare_means(comparisons = my_comparisons, method = "t.test")             + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))

    # Arrange density plot
    boxplots_plot<-ggarrange(m4, m5, m6, m7, nrow = 1,ncol = 4, common.legend = TRUE, legend="bottom")    


    # FindClusters_resolution          
    png(filename=paste(output_dir,"countour_T2_Coonections_melt_",normalization_scheme,"_",TCGA_project,"_Stage_all_boxplots.png",sep=""), width = 30, height = 15, res=600, units = "cm")      
        boxplots_plot
    dev.off()
    
    
    #########################################################################################################################################
    # Three countour plots will be created
    # One with the average expression
    # Second with the expression per patient
    # Third with the z-score 
    # Conections, GC3, AvgExpression
    # Combine AvgExpression, Conections, T2
    # harmonic mean
    # z-core : Composite Scores
    # Implementing the Z score formula in R is quite straightforward. 
    # To reuse code, we will create a function called calculate_z using the mean and sd base functions to calculate Z. 
    # sd calculates the standard deviation in R.
    # weighted average
    # Z-score for AveExp_expression
    # Z-score for AveExp_expression  
    m1<-ggplot(Interactomes_GC3_T2_merged_Stage_I, aes(Conections, GC3, z = AveExp))   + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85))      + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  
    m2<-ggplot(Interactomes_GC3_T2_merged_Stage_II, aes(Conections, GC3, z = AveExp))   + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage II Expr. ",sep=""))      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85))    + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))
    m3<-ggplot(Interactomes_GC3_T2_merged_Stage_III, aes(Conections, GC3, z = AveExp))   + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage III Expr. ",sep=""))      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85))  + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))

    # Arrange density plot
    density_plot<-ggarrange(m1, m2, m3, nrow = 1,ncol = 3, common.legend = TRUE, legend="bottom")
    
    m7<-ggplot(Interactomes_GC3_T2_merged_Stage_I, aes(Conections, GC3, z = AveExp))    + geom_point(aes(colour=AveExp),size=3)   +  theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep="")) + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85))  + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))                       
    m8<-ggplot(Interactomes_GC3_T2_merged_Stage_II, aes(Conections, GC3, z = AveExp))   + geom_point(aes(colour=AveExp),size=3)  +  theme_bw() + ggtitle(paste(normalization_scheme,   ": All points Stage II Expr. ",sep=""))  + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85))  + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))
    m9<-ggplot(Interactomes_GC3_T2_merged_Stage_III, aes(Conections, GC3, z = AveExp))  + geom_point(aes(colour=AveExp),size=3) + theme_bw() + ggtitle(paste(normalization_scheme,  ": All points Stage III Expr. ",sep=""))    + xlim(0, 50)  + ylim(25,85)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85))  + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))

    # Arrange density plot
    dotplot_plot<-ggarrange(m7, m8, m9, nrow = 1,ncol = 3, common.legend = TRUE, legend="bottom")        
  
    m10<-ggplot(Interactomes_GC3_T2_merged_Stage_I, aes(Conections, GC3, z = AveExp))   +  theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))   + geom_density_2d(bin=10)           +   guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE)) + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85)) + geom_point(aes(colour=AveExp), size=3)  + scale_colour_gradient(low = "red", high = "yellow")    + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))
    m11<-ggplot(Interactomes_GC3_T2_merged_Stage_II, aes(Conections, GC3, z = AveExp))   +  theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage II Expr. ",sep=""))   + geom_density_2d(bin=10)           +   guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE)) + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85)) + geom_point(aes(colour=AveExp), size=3) + scale_colour_gradient(low = "red", high = "yellow")   + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))
    m12<-ggplot(Interactomes_GC3_T2_merged_Stage_III, aes(Conections, GC3, z = AveExp))   +  theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage III Expr. ",sep=""))   + geom_density_2d(bin=10)           +   guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE)) + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85)) + geom_point(aes(colour=AveExp), size=3) + scale_colour_gradient(low = "red", high = "yellow") + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))

    # Arrange density plot
    countour_plot<-ggarrange(m10, m11, m12, nrow = 1,ncol = 3, common.legend = TRUE, legend="bottom")            
      
    # FindClusters_resolution               
    png(filename=paste(output_dir,"countour_T2_Coonections_melt_",normalization_scheme,"_",TCGA_project,"_Stage_all_GC3.png",sep=""), width = 30, height = 30, res=600, units = "cm")  
            plot<-ggarrange(countour_plot,density_plot,  nrow = 2,ncol = 1, common.legend = TRUE, legend="bottom")
            print(annotate_figure(plot, top = text_grob(TCGA_project, face = "bold", size = 14)))
    dev.off()
    
    h1<-ggplot(unique(Interactomes_GC3_T2_merged_all[,c("ENSEMBL","T2","GC3","Conections","Stages")]), aes(x=T2, color=Stages, after_stat(count))) + geom_density() +   theme_bw()     + xlim(0, 50) + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + ggtitle("A")
    h2<-ggplot(unique(Interactomes_GC3_T2_merged_all[,c("ENSEMBL","T2","GC3","Conections","Stages")]), aes(x=GC3, color=Stages, after_stat(count))) +geom_density() +  theme_bw()    + xlim(0, 100)  + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  + ggtitle("B")
    h3<-ggplot(unique(Interactomes_GC3_T2_merged_all[,c("ENSEMBL","T2","AveExp","Conections","Stages")]), aes(x=AveExp, color=Stages, after_stat(count))) + geom_density() +  theme_bw()  + xlim(0, 10000) + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  
    h4<-ggplot(unique(Interactomes_GC3_T2_merged_all[,c("ENSEMBL","T2","AveExp","Conections","Stages")]), aes(x=Conections, color=Stages, after_stat(count))) + geom_density() +  theme_bw()  + xlim(0, 10000) + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  

    i1<-ggplot(unique(Interactomes_GC3_T2_merged_all[,c("ENSEMBL","T2","GC3","Conections","Stages")]), aes(x=T2, color=Stages)) + geom_histogram(fill="white", alpha=0.5, position="identity", bins=20)  +   theme_bw()     + xlim(0, 50) + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  
    i2<-ggplot(unique(Interactomes_GC3_T2_merged_all[,c("ENSEMBL","T2","GC3","Conections","Stages")]), aes(x=GC3, color=Stages)) + geom_histogram(fill="white", alpha=0.5, position="identity", bins=20)  +  theme_bw()    + xlim(0, 100) + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  
    i3<-ggplot(unique(Interactomes_GC3_T2_merged_all[,c("ENSEMBL","T2","AveExp","Conections","Stages")]), aes(x=AveExp, color=Stages)) + geom_histogram(fill="white", alpha=0.5, position="identity", bins=20)  +  theme_bw()  + xlim(0, 10000) + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  
    i4<-ggplot(unique(Interactomes_GC3_T2_merged_all[,c("ENSEMBL","T2","AveExp","Conections","Stages")]), aes(x=Conections, color=Stages)) + geom_histogram(fill="white", alpha=0.5, position="identity", bins=20)  +  theme_bw()  + xlim(0, 10000)  + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  

    p1<-ggplot(unique(Interactomes_GC3_T2_merged_all[,c("ENSEMBL","T2","GC3","Conections","Stages")]), aes(x=T2, color=Stages)) + geom_histogram(fill="white", alpha=0.5, position="identity", bins=20, aes(y = after_stat(count / sum(count))))  +   theme_bw()     + xlim(0, 50) + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  +   scale_y_continuous(labels = scales::percent)  + ylab("relative frequency (%)") + ggtitle("A")
    p2<-ggplot(unique(Interactomes_GC3_T2_merged_all[,c("ENSEMBL","T2","GC3","Conections","Stages")]), aes(x=GC3, color=Stages)) + geom_histogram(fill="white", alpha=0.5, position="identity", bins=20, aes(y = after_stat(count / sum(count))))  +  theme_bw()    + xlim(0, 100) + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   +   scale_y_continuous(labels = scales::percent) + ylab("")                       + ggtitle("B")
    
    # FindClusters_resolution               
    png(filename=paste(output_dir,"countour_T2_Coonections_melt_",normalization_scheme,"_",TCGA_project,"_paper.png",sep=""), width = 20, height = 12, res=1200, units = "cm")            
          ggarrange(p1, p2, nrow = 1, ncol=2, common.legend = TRUE, legend="bottom") 
    dev.off()        

    
    # FindClusters_resolution               
    png(filename=paste(output_dir,"countour_T2_Coonections_melt_",normalization_scheme,"_",TCGA_project,"_unique_histogram.png",sep=""), width = 30, height = 15, res=600, units = "cm")            
          plot<-ggarrange(h1, h2, h3,h4,i1, i2, i3, i4, nrow = 2, ncol=4, common.legend = TRUE, legend="bottom") 
          annotate_figure(plot, top = text_grob(paste(TCGA_project,normalization_scheme,sep=" "), face = "bold", size = 14))  
    dev.off()        

    #########################################################################################################################################
    # Implementing the Z score formula in R is quite straightforward. 
    # To reuse code, we will create a function called calculate_z using the mean and sd base functions to calculate Z. 
    # sd calculates the standard deviation in R.
    # weighted average    # FindClusters_resolution               
    png(filename=paste(output_dir,"countour_T2_Coonections_melt_",normalization_scheme,"_",TCGA_project,"_unique_histogram.png",sep=""), width = 10, height = 20, res=1200, units = "cm")            
          plot<-ggarrange(h1, h2, nrow = 1, ncol=2, common.legend = TRUE, legend="bottom") 
    dev.off()        

    
    # FindClusters_resolution               
    png(filename=paste(output_dir,"countour_T2_Coonections_melt_",normalization_scheme,"_",TCGA_project,"_unique_histogram.png",sep=""), width = 30, height = 15, res=600, units = "cm")            
          plot<-ggarrange(h1, h2, h3,h4,i1, i2, i3, i4, nrow = 2, ncol=4, common.legend = TRUE, legend="bottom") 
          annotate_figure(plot, top = text_grob(paste(TCGA_project,normalization_scheme,sep=" "), face = "bold", size = 14))  
    dev.off()        

    #########################################################################################################################################
    # Implementing the Z score formula in R is quite straightforward. 
    # To reuse code, we will create a function called calculate_z using the mean and sd base functions to calculate Z. 
    # sd calculates the standard deviation in R.
    # weighted average
    # Z-score for AveExp_expression
    # Z-score for AveExp_expression       
    m1<-ggplot(Interactomes_GC3_T2_merged_Stage_I_patient, aes(Conections, T2, z = AveExp))   + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))      + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + theme(legend.position = "none")
    m2<-ggplot(Interactomes_GC3_T2_merged_Stage_II_patient, aes(Conections, T2, z = AveExp))   + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage II Expr. ",sep=""))      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))    + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + theme(legend.position = "none")
    m3<-ggplot(Interactomes_GC3_T2_merged_Stage_III_patient, aes(Conections, T2, z = AveExp))   + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage III Expr. ",sep=""))      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))  + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + theme(legend.position = "none")

    # Arrange density plot
    density_plot<-ggarrange(m1, m2, m3, nrow = 1,ncol = 3, common.legend = TRUE, legend="none")
  
    m4 <- ggplot(expression_table_normalized_melt, aes(x=Stages, y=T2)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()            + ggtitle(paste("T2 All points: ", normalization_scheme,sep=""))                   +  stat_compare_means(comparisons = my_comparisons, method = "t.test")  + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  
    m5 <- ggplot(expression_table_normalized_melt, aes(x=Stages, y=Conections)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()    + ggtitle(paste("Connectivity All points: ", normalization_scheme,sep="")) +  stat_compare_means(comparisons = my_comparisons, method = "t.test")          + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  
    m6 <- ggplot(expression_table_normalized_melt, aes(x=Stages, y=AveExp)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()        + ggtitle(paste("Expr. All points: ", normalization_scheme,sep=""))          +  stat_compare_means(comparisons = my_comparisons, method = "t.test")        + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  

    # Arrange density plot
    boxplots_plot<-ggarrange(m4, m5, m6, nrow = 1,ncol = 3, common.legend = TRUE, legend="bottom")    
  
    m7<-ggplot(Interactomes_GC3_T2_merged_Stage_I_patient, aes(Conections, T2, z = AveExp))    + geom_point(aes(colour=AveExp),size=3)   +  theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))+  xlim(0, 50)     + ylim(10,40) + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(0, 40, by = 10), limits = c(10, 40))      + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  
    m8<-ggplot(Interactomes_GC3_T2_merged_Stage_II_patient, aes(Conections, T2, z = AveExp))   + geom_point(aes(colour=AveExp),size=3)  +  theme_bw() + ggtitle(paste(normalization_scheme,   ": All points Stage II Expr. ",sep=""))+ xlim(0, 50)       + ylim(10,40) + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(0, 40, by = 10), limits = c(10, 40))      + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  
    m9<-ggplot(Interactomes_GC3_T2_merged_Stage_III_patient, aes(Conections, T2, z = AveExp))  + geom_point(aes(colour=AveExp),size=3) + theme_bw() + ggtitle(paste(normalization_scheme,  ": All points Stage III Expr. ",sep=""))+ xlim(0, 50)         + ylim(10,40) +  guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE)) + ylim(10,40)  + xlim(0, 50)     + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(0, 40, by = 10), limits = c(10, 40))   + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  

 
    # Z-score for AveExp_expression
    # Z-score for AveExp_expression       
    m1<-ggplot(Interactomes_GC3_T2_merged_Stage_I_patient, aes(Conections, T2, z = AveExp))   + geom_density_2d_filled() + theme_bw() + ggtitle("Stage I")      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))      + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + geom_hline(yintercept=30, colour="yellow") + theme(legend.position = "none") + ggtitle("F")
    m2<-ggplot(Interactomes_GC3_T2_merged_Stage_II_patient, aes(Conections, T2, z = AveExp))   + geom_density_2d_filled() + theme_bw() + ggtitle("Stage II")      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))    + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + geom_hline(yintercept=30, colour="yellow") + theme(legend.position = "none") + ggtitle("G")  
    m3<-ggplot(Interactomes_GC3_T2_merged_Stage_III_patient, aes(Conections, T2, z = AveExp))   + geom_density_2d_filled() + theme_bw() + ggtitle("Stage III")      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))  + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + geom_hline(yintercept=30, colour="yellow") + theme(legend.position = "none") + ggtitle("H")  
    
    # Arrange density plot
    density_plot<-ggarrange(m1, m2, m3, nrow = 1,ncol = 3, common.legend = TRUE, legend="none")
    
    m7<-ggplot(Interactomes_GC3_T2_merged_Stage_I_patient, aes(Conections, T2, z = AveExp))    + geom_point(aes(colour=AveExp),size=3)   +  theme_bw() + ggtitle("Stage I")+  xlim(0, 50)     + ylim(10,40) + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(0, 40, by = 10), limits = c(10, 40))      + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + scale_colour_gradient(low = "red", high = "yellow") + geom_hline(yintercept=30, colour="yellow")   + ggtitle("A")
    m8<-ggplot(Interactomes_GC3_T2_merged_Stage_II_patient, aes(Conections, T2, z = AveExp))   + geom_point(aes(colour=AveExp),size=3)  +  theme_bw() + ggtitle("Stage II")+ xlim(0, 50)       + ylim(10,40) + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(0, 40, by = 10), limits = c(10, 40))      + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  + scale_colour_gradient(low = "red", high = "yellow") + geom_hline(yintercept=30, colour="yellow")   + ggtitle("B")
    m9<-ggplot(Interactomes_GC3_T2_merged_Stage_III_patient, aes(Conections, T2, z = AveExp))  + geom_point(aes(colour=AveExp),size=3) + theme_bw() + ggtitle("Stage III")+ xlim(0, 50)         + ylim(10,40) +  guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE)) + ylim(10,40)  + xlim(0, 50)     + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(0, 40, by = 10), limits = c(10, 40))   + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + scale_colour_gradient(low = "red", high = "yellow") + geom_hline(yintercept=30, colour="yellow") + ggtitle("C")

    # Arrange density plot
    dotplot_plot<-ggarrange(m7, m8, m9, nrow = 1,ncol = 3, common.legend = TRUE, legend="bottom")        
  
    m10<-ggplot(Interactomes_GC3_T2_merged_Stage_I_patient, aes(Conections, T2, z = AveExp))   +  theme_bw() + ggtitle("Stage I")   + geom_density_2d(bin=10)      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE)) + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(0, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))   + scale_colour_gradient(low = "red", high = "yellow")   + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  + theme(legend.position = "none") + geom_hline(yintercept=30, colour="yellow") + ggtitle("D")
    m11<-ggplot(Interactomes_GC3_T2_merged_Stage_II_patient, aes(Conections, T2, z = AveExp))  +  theme_bw() + ggtitle("Stage II")   + geom_density_2d(bin=10)      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(0, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))   + scale_colour_gradient(low = "red", high = "yellow")   + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) + theme(legend.position = "none") + geom_hline(yintercept=30, colour="yellow") + ggtitle("E")
    m12<-ggplot(Interactomes_GC3_T2_merged_Stage_III_patient, aes(Conections, T2, z = AveExp))  +  theme_bw() + ggtitle("Stage III")   + geom_density_2d(bin=10)      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(0, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40)) + scale_colour_gradient(low = "red", high = "yellow") + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + theme(legend.position = "none") + geom_hline(yintercept=30, colour="yellow") + ggtitle("F")

    # Arrange density plot
    countour_plot<-ggarrange(m10, m11, m12, nrow = 1,ncol = 3, common.legend = TRUE, legend="bottom")            
      
    # FindClusters_resolution          
    png(filename=paste(output_dir,"countour_T2_Coonections_melt_",normalization_scheme,"_",TCGA_project,"_Stage_all_T2_perPatient_paper.png",sep=""), width = 30, height = 30, res=600, units = "cm")  
            ggarrange(dotplot_plot, countour_plot, density_plot,  nrow = 3,ncol = 1, common.legend = TRUE, legend="bottom")
            #print(annotate_figure(plot, top = text_grob(TCGA_project, face = "bold", size = 14)))
    dev.off() 
    ############################################################################################################################
    write_tsv(Interactomes_GC3_T2_merged_Stage_I_patient, paste(output_dir,"/Interactomes_GC3_T2_merged_Stage_I_patient_",normalization_scheme,"_Stage_I",".tsv",sep=""))					
    write_tsv(Interactomes_GC3_T2_merged_Stage_II_patient, paste(output_dir,"/Interactomes_GC3_T2_merged_Stage_II_patient_",normalization_scheme,"_Stage_II",".tsv",sep=""))					
    write_tsv(Interactomes_GC3_T2_merged_Stage_III_patient, paste(output_dir,"/Interactomes_GC3_T2_merged_Stage_III_patient_",normalization_scheme,"_Stage_III",".tsv",sep=""))					
    ############################################################################################################################



    ############################################################################################################################

    # Z-score for AveExp_expression
    # Z-score for AveExp_expression       
    m1<-ggplot(Interactomes_GC3_T2_merged_Stage_I_patient, aes(Conections, T2, z = AveExp))   + geom_density_2d_filled() + theme_bw() + ggtitle("Stage I")      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))      + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + geom_hline(yintercept=30, colour="yellow")
    m2<-ggplot(Interactomes_GC3_T2_merged_Stage_II_patient, aes(Conections, T2, z = AveExp))   + geom_density_2d_filled() + theme_bw() + ggtitle("Stage II")      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))    + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + geom_hline(yintercept=30, colour="yellow")
    m3<-ggplot(Interactomes_GC3_T2_merged_Stage_III_patient, aes(Conections, T2, z = AveExp))   + geom_density_2d_filled() + theme_bw() + ggtitle("Stage III")      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))  + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   + geom_hline(yintercept=30, colour="yellow")

    # Arrange density plot
    density_plot<-ggarrange(m1, m2, m3, nrow = 1,ncol = 3, common.legend = TRUE, legend="bottom")
    

    # FindClusters_resolution          
    png(filename=paste(output_dir,"countour_T2_Coonections_melt_",normalization_scheme,"_",TCGA_project,"_Stage_all_T2_perPatient_paper.png",sep=""), width = 24, height = 12, res=600, units = "cm")  
            density_plot
    dev.off()        
    
    
    # Implementing the Z score formula in R is quite straightforward. 
    # To reuse code, we will create a function called calculate_z using the mean and sd base functions to calculate Z. 
    # sd calculates the standard deviation in R.
    # weighted average
    # Z-score for AveExp_expression
    # Z-score for AveExp_expression       
    m1<-ggplot(Interactomes_GC3_T2_merged_Stage_I_patient, aes(Conections, GC3, z = AveExp))   + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85)) + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  
    m2<-ggplot(Interactomes_GC3_T2_merged_Stage_II_patient, aes(Conections, GC3, z = AveExp))   + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage II Expr. ",sep=""))      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85)) + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  
    m3<-ggplot(Interactomes_GC3_T2_merged_Stage_III_patient, aes(Conections, GC3, z = AveExp))   + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage III Expr. ",sep=""))      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85)) + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  

    # Arrange density plot
    density_plot<-ggarrange(m1, m2, m3, nrow = 1,ncol = 3, common.legend = TRUE, legend="bottom")
  
    m4 <- ggplot(expression_table_normalized_melt, aes(x=Stages, y=GC3)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()            + ggtitle(paste("GC3 All points: ", normalization_scheme,sep=""))                   +  stat_compare_means(comparisons = my_comparisons, method = "t.test")   + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  
    m5 <- ggplot(expression_table_normalized_melt, aes(x=Stages, y=Conections)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()    + ggtitle(paste("Connectivity All points: ", normalization_scheme,sep="")) +  stat_compare_means(comparisons = my_comparisons, method = "t.test")  + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  
    m6 <- ggplot(expression_table_normalized_melt, aes(x=Stages, y=AveExp)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()        + ggtitle(paste("Expr. All points: ", normalization_scheme,sep=""))          +  stat_compare_means(comparisons = my_comparisons, method = "t.test")  + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))        

    # Arrange density plot
    boxplots_plot<-ggarrange(m4, m5, m6, nrow = 1,ncol = 3, common.legend = TRUE, legend="bottom")    
  
    m7<-ggplot(Interactomes_GC3_T2_merged_Stage_I_patient, aes(Conections, GC3, z = AveExp))    + geom_point(aes(colour=AveExp),size=3)   +  theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))+  xlim(0, 50)     + ylim(25,85) + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(25,85)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(0, 40, by = 10), limits = c(25, 85)) + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  
    m8<-ggplot(Interactomes_GC3_T2_merged_Stage_II_patient, aes(Conections, GC3, z = AveExp))   + geom_point(aes(colour=AveExp),size=3)  +  theme_bw() + ggtitle(paste(normalization_scheme,   ": All points Stage II Expr. ",sep=""))+ xlim(0, 50)       + ylim(25,85) + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(25,85)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(0, 40, by = 10), limits = c(25, 85)) + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  
    m9<-ggplot(Interactomes_GC3_T2_merged_Stage_III_patient, aes(Conections, GC3, z = AveExp))  + geom_point(aes(colour=AveExp),size=3) + theme_bw() + ggtitle(paste(normalization_scheme,  ": All points Stage III Expr. ",sep=""))+ xlim(0, 50)         + ylim(25,85) +  guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE)) + ylim(25,85)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(0, 40, by = 10), limits = c(25, 85)) + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  

    # Arrange density plot
    dotplot_plot<-ggarrange(m7, m8, m9, nrow = 1,ncol = 3, common.legend = TRUE, legend="bottom")        
  
    m10<-ggplot(Interactomes_GC3_T2_merged_Stage_I_patient, aes(Conections,  GC3, z = AveExp))   +  theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))   + geom_density_2d(bin=10)      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE)) + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(0, 40, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85))  + geom_point(aes(colour=AveExp), size=3)  + scale_colour_gradient(low = "red", high = "yellow") + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  
    m11<-ggplot(Interactomes_GC3_T2_merged_Stage_II_patient, aes(Conections, GC3, z = AveExp))  +  theme_bw()  + ggtitle(paste(normalization_scheme,   ": All points Stage II Expr. ",sep=""))   + geom_density_2d(bin=10)      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(25,85)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(0, 40, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85))  + geom_point(aes(colour=AveExp), size=3)  + scale_colour_gradient(low = "red", high = "yellow") + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  
    m12<-ggplot(Interactomes_GC3_T2_merged_Stage_III_patient, aes(Conections,GC3, z = AveExp))  +  theme_bw()  + ggtitle(paste(normalization_scheme,   ": All points Stage III Expr. ",sep=""))   + geom_density_2d(bin=10)      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(25,85)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(0, 40, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85))  + geom_point(aes(colour=AveExp), size=3)  + scale_colour_gradient(low = "red", high = "yellow") + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  

    # Arrange density plot
    countour_plot<-ggarrange(m10, m11, m12, nrow = 1,ncol = 3, common.legend = TRUE, legend="bottom")            
      

    # FindClusters_resolution          
    png(filename=paste(output_dir,"countour_T2_Coonections_melt_",normalization_scheme,"_",TCGA_project,"_Stage_all_GC3_perPatient.png",sep=""), width = 30, height = 30, res=600, units = "cm")  
            plot<-ggarrange(countour_plot,density_plot,  nrow = 2,ncol = 1, common.legend = TRUE, legend="bottom")
            print(annotate_figure(plot, top = text_grob(TCGA_project, face = "bold", size = 14)))
    dev.off()       



    m4 <- ggplot(expression_table_normalized_melt, aes(x=Stages, y=T2)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()            + ggtitle(paste("T2", normalization_scheme,sep=" : "))                   +  stat_compare_means(comparisons = my_comparisons, method = "t.test")  + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))
    m5 <- ggplot(expression_table_normalized_melt, aes(x=Stages, y=Conections)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()    + ggtitle(paste("Connectivity", normalization_scheme,sep=" : ")) +  stat_compare_means(comparisons = my_comparisons, method = "t.test")          + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) 
    m6 <- ggplot(expression_table_normalized_melt, aes(x=Stages, y=AveExp)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()        + ggtitle(paste("Expr.", normalization_scheme,sep=" : "))          +  stat_compare_means(comparisons = my_comparisons, method = "t.test")        + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))
    m7 <- ggplot(expression_table_normalized_melt, aes(x=Stages, y=GC3)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()           + ggtitle(paste("GC3", normalization_scheme,sep=" : "))          +  stat_compare_means(comparisons = my_comparisons, method = "t.test")             + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))

    # Arrange density plot
    boxplots_plot<-ggarrange(m4, m5, m6, m7, nrow = 1,ncol = 4, common.legend = TRUE, legend="bottom")    



    # FindClusters_resolution          
    png(filename=paste(output_dir,"countour_T2_Coonections_melt_",normalization_scheme,"_",TCGA_project,"_Stage_all_boxplots_perpatient.png",sep=""), width = 30, height = 15, res=600, units = "cm")  
        print(boxplots_plot)
    dev.off()
    
    h1<-ggplot(unique(expression_table_normalized_melt[,c("ENSEMBL","T2","GC3","Conections","Stages")]), aes(x=T2, color=Stages,after_stat(count))) + geom_density() +   theme_bw()     + xlim(0, 50)  + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  
    h2<-ggplot(unique(expression_table_normalized_melt[,c("ENSEMBL","T2","GC3","Conections","Stages")]), aes(x=GC3, color=Stages,after_stat(count))) +geom_density() +  theme_bw()    + xlim(0, 100)   + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  
    h3<-ggplot(unique(expression_table_normalized_melt[,c("ENSEMBL","T2","AveExp","Conections","Stages")]), aes(x=AveExp, color=Stages,after_stat(count))) + geom_density() +  theme_bw()  + xlim(0, 10000) + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  
    h4<-ggplot(unique(expression_table_normalized_melt[,c("ENSEMBL","T2","AveExp","Conections","Stages")]), aes(x=Conections, color=Stages,after_stat(count))) + geom_density() +  theme_bw()  + xlim(0, 10000) + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  

    i1<-ggplot(unique(expression_table_normalized_melt[,c("ENSEMBL","T2","GC3","Conections","Stages")]), aes(x=T2, color=Stages)) + geom_histogram(fill="white", alpha=0.5, position="identity", bins=20)  +   theme_bw()     + xlim(0, 50) + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  
    i2<-ggplot(unique(expression_table_normalized_melt[,c("ENSEMBL","T2","GC3","Conections","Stages")]), aes(x=GC3, color=Stages)) + geom_histogram(fill="white", alpha=0.5, position="identity", bins=20)  +  theme_bw()    + xlim(0, 100) + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))   
    i3<-ggplot(unique(expression_table_normalized_melt[,c("ENSEMBL","T2","AveExp","Conections","Stages")]), aes(x=AveExp, color=Stages)) + geom_histogram(fill="white", alpha=0.5, position="identity", bins=20)  +  theme_bw()  + xlim(0, 10000) + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))  
    i4<-ggplot(unique(expression_table_normalized_melt[,c("ENSEMBL","T2","AveExp","Conections","Stages")]), aes(x=Conections, color=Stages)) + geom_histogram(fill="white", alpha=0.5, position="identity", bins=20)  +  theme_bw()  + xlim(0, 10000) + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))     
    
    # FindClusters_resolution               
    png(filename=paste(output_dir,"countour_T2_Coonections_melt_",normalization_scheme,"_",TCGA_project,"_unique_histogram_perpatient.png",sep=""), width = 25, height = 25, res=600, units = "cm")            
          plot<-ggarrange(h1, h2, h3,h4,i1, i2, i3, i4, nrow = 2, ncol=4, common.legend = TRUE, legend="bottom") 
          annotate_figure(plot, top = text_grob(paste(TCGA_project,normalization_scheme,sep=" "), face = "bold", size = 14))  
    dev.off()            
}
