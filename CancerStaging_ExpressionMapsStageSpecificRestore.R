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

    unique_stage_I  =intersect(setdiff(genes_stages_I, c(genes_stages_II,genes_stages_III)),genes_stages_I)
    unique_stage_II =intersect(setdiff(genes_stages_II, c(genes_stages_I,genes_stages_III)),genes_stages_II)
    unique_stage_III=intersect(setdiff(genes_stages_III, c(genes_stages_I,genes_stages_II)),genes_stages_III)

    Interactomes_GC3_T2_merged_Stage_I<-na.omit(Interactomes_GC3_T2_merged_Stage_I[unique_stage_I,])
    Interactomes_GC3_T2_merged_Stage_II<-na.omit(Interactomes_GC3_T2_merged_Stage_II[unique_stage_II,])
    Interactomes_GC3_T2_merged_Stage_III<-na.omit(Interactomes_GC3_T2_merged_Stage_III[unique_stage_III,])

    # Merge Stages
    Interactomes_GC3_T2_merged_all<-na.omit(rbind(Interactomes_GC3_T2_merged_Stage_I,Interactomes_GC3_T2_merged_Stage_II, Interactomes_GC3_T2_merged_Stage_III))  
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
    m1<-ggplot(Interactomes_GC3_T2_merged_Stage_I, aes(Conections, T2, z = AveExp))   + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))
    m2<-ggplot(Interactomes_GC3_T2_merged_Stage_II, aes(Conections, T2, z = AveExp))   + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage II Expr. ",sep=""))      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))
    m3<-ggplot(Interactomes_GC3_T2_merged_Stage_III, aes(Conections, T2, z = AveExp))   + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage III Expr. ",sep=""))      + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))

    # Arrange density plot
    density_plot<-ggarrange(m1, m2, m3, nrow = 1,ncol = 3, common.legend = TRUE, legend="bottom")
  
    m4 <- ggplot(Interactomes_GC3_T2_merged_all, aes(x=Stages, y=T2)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()           +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_all$T2), linetype="dashed", color = "red")    + ggtitle(paste("T2 All points: ", normalization_scheme,sep=""))                   +  stat_compare_means(comparisons = my_comparisons, method = "t.test") 
    m5 <- ggplot(Interactomes_GC3_T2_merged_all, aes(x=Stages, y=Conections)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()   +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_all$Conections), linetype="dashed", color = "red")    + ggtitle(paste("Connectivity All points: ", normalization_scheme,sep="")) +  stat_compare_means(comparisons = my_comparisons, method = "t.test") 
    m6 <- ggplot(Interactomes_GC3_T2_merged_all, aes(x=Stages, y=AveExp)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()         +      geom_hline(yintercept=median(Interactomes_GC3_T2_merged_all$Expr), linetype="dashed", color = "red")    + ggtitle(paste("Expr. All points: ", normalization_scheme,sep=""))          +  stat_compare_means(comparisons = my_comparisons, method = "t.test")       

    # Arrange density plot
    boxplots_plot<-ggarrange(m4, m5, m6, nrow = 1,ncol = 3, common.legend = TRUE, legend="bottom")    
  
    m7<-ggplot(Interactomes_GC3_T2_merged_Stage_I, aes(Conections, T2, z = AveExp))    + geom_point(aes(colour=AveExp))   +  theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))+  xlim(0, 50)     + ylim(10,40) + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(0, 40, by = 10), limits = c(10, 40))
    m8<-ggplot(Interactomes_GC3_T2_merged_Stage_II, aes(Conections, T2, z = AveExp))   + geom_point(aes(colour=AveExp))  +  theme_bw() + ggtitle(paste(normalization_scheme,   ": All points Stage II Expr. ",sep=""))+ xlim(0, 50)       + ylim(10,40) + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(0, 40, by = 10), limits = c(10, 40))
    m9<-ggplot(Interactomes_GC3_T2_merged_Stage_III, aes(Conections, T2, z = AveExp))  + geom_point(aes(colour=AveExp)) + theme_bw() + ggtitle(paste(normalization_scheme,  ": All points Stage III Expr. ",sep=""))+ xlim(0, 50)         + ylim(10,40) +  guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE)) + ylim(10,40)  + xlim(0, 50)     + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(0, 40, by = 10), limits = c(10, 40))

    # Arrange density plot
    dotplot_plot<-ggarrange(m7, m8, m9, nrow = 1,ncol = 3, common.legend = TRUE, legend="bottom")        
  
    m10<-ggplot(Interactomes_GC3_T2_merged_Stage_I, aes(Conections, T2, z = AveExp))   +  theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))   + geom_density_2d(bin=10)         + geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_I$Conections), linetype="dashed", color = "red")   +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_I$T2), linetype="dashed", color = "red")     + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_I$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_I$T2), linetype="dashed", color =   "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE)) + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(0, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))
    m11<-ggplot(Interactomes_GC3_T2_merged_Stage_II, aes(Conections, T2, z = AveExp))  +  theme_bw() + ggtitle(paste(normalization_scheme,   ": All points Stage II Expr. ",sep=""))   + geom_density_2d(bin=10)        + geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_II$Conections), linetype="dashed", color = "red")  +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_II$T2), linetype="dashed", color = "red")   + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_II$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_II$T2), linetype="dashed", color =  "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(0, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))
    m12<-ggplot(Interactomes_GC3_T2_merged_Stage_III, aes(Conections, T2, z = AveExp)) +  theme_bw() + ggtitle(paste(normalization_scheme,  ": All points Stage III Expr. ",sep=""))   + geom_density_2d(bin=10)       + geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_III$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_III$T2), linetype="dashed", color = "red") + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_III$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_III$T2), linetype="dashed", color = "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(0, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))  

    # Arrange density plot
    countour_plot<-ggarrange(m10, m11, m12, nrow = 1,ncol = 3, common.legend = TRUE, legend="bottom")            
      
    # FindClusters_resolution               
    png(filename=paste(output_dir,"countour_T2_Coonections_melt_",normalization_scheme,"_",TCGA_project,"_Stage_all_T2.png",sep=""), width = 30, height = 50, res=600, units = "cm")  
            plot<-ggarrange(m4, m5, m6, m7, m8, m9, m10, m11,m12, m1,m2,m3, nrow = 4,ncol = 3, common.legend = TRUE, legend="bottom")
            print(annotate_figure(plot, top = text_grob(TCGA_project, face = "bold", size = 14)))
    dev.off()
    #########################################################################################################################################
    m1<-ggplot(Interactomes_GC3_T2_merged_Stage_I, aes(Conections, GC3, z = AveExp))   + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))+  geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_I$Conections), linetype="dashed", color = "red")   +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_I$GC3), linetype="dashed", color = "red")     + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_I$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_I$GC3), linetype="dashed", color =   "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85)) 
    m2<-ggplot(Interactomes_GC3_T2_merged_Stage_II, aes(Conections, GC3, z = AveExp))  + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,   ": All points Stage II Expr. ",sep=""))+  geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_II$Conections), linetype="dashed", color = "red")  +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_II$GC3), linetype="dashed", color = "red")   + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_II$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_II$GC3), linetype="dashed", color =  "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(25,85)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85))
    m3<-ggplot(Interactomes_GC3_T2_merged_Stage_III, aes(Conections, GC3, z = AveExp)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,  ": All points Stage III Expr. ",sep=""))+  geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_III$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_III$GC3), linetype="dashed", color = "red") + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_III$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_III$GC3), linetype="dashed", color = "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(25,85)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85))
  
    m4 <- ggplot(Interactomes_GC3_T2_merged_all, aes(x=Stages, y=GC3)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()           +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_all$GC3), linetype="dashed", color = "red")    + ggtitle(paste("T2 All points: ", normalization_scheme,sep=""))                   +  stat_compare_means(comparisons = my_comparisons, method = "t.test") 
    m5 <- ggplot(Interactomes_GC3_T2_merged_all, aes(x=Stages, y=Conections)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()   +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_all$Conections), linetype="dashed", color = "red")    + ggtitle(paste("Connectivity All points: ", normalization_scheme,sep="")) +  stat_compare_means(comparisons = my_comparisons, method = "t.test") 
    m6 <- ggplot(Interactomes_GC3_T2_merged_all, aes(x=Stages, y=AveExp)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()         +      geom_hline(yintercept=median(Interactomes_GC3_T2_merged_all$Expr), linetype="dashed", color = "red")    + ggtitle(paste("Expr. All points: ", normalization_scheme,sep=""))          +  stat_compare_means(comparisons = my_comparisons, method = "t.test")       
  
    m7<-ggplot(Interactomes_GC3_T2_merged_Stage_I, aes(Conections, GC3, z = AveExp))    + geom_point(aes(colour=AveExp))   +  theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))+  xlim(0, 50)       + geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_I$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_I$GC3), linetype="dashed", color = "red")     + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_I$Conections), linetype="dashed", color = "yellow") +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_I$GC3), linetype="dashed", color =   "yellow")   + ylim(25,85) + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85)) 
    m8<-ggplot(Interactomes_GC3_T2_merged_Stage_II, aes(Conections, GC3, z = AveExp))   + geom_point(aes(colour=AveExp))  +  theme_bw() + ggtitle(paste(normalization_scheme,   ": All points Stage II Expr. ",sep=""))+ xlim(0, 50) + geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_II$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_II$GC3), linetype="dashed", color = "red")   + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_II$Conections), linetype="dashed", color = "yellow") +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_II$GC3), linetype="dashed", color =  "yellow")   + ylim(25,85) + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))         + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85)) 
    m9<-ggplot(Interactomes_GC3_T2_merged_Stage_III, aes(Conections, GC3, z = AveExp))  + geom_point(aes(colour=AveExp)) + theme_bw() + ggtitle(paste(normalization_scheme,  ": All points Stage III Expr. ",sep=""))+ xlim(0, 50) +  geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_III$Conections), linetype="dashed", color = "red")  +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_III$GC3), linetype="dashed", color = "red") + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_III$Conections), linetype="dashed", color = "yellow") +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_III$GC3), linetype="dashed", color = "yellow")   + ylim(25,85) +  guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))       + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85)) 
  
    m10<-ggplot(Interactomes_GC3_T2_merged_Stage_I, aes(Conections, GC3, z = AveExp))   +  theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))   + geom_density_2d(bin=10)         + geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_I$Conections), linetype="dashed", color = "red")   +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_I$GC3), linetype="dashed", color = "red")     + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_I$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_I$GC3), linetype="dashed", color =   "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE)) + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85)) 
    m11<-ggplot(Interactomes_GC3_T2_merged_Stage_II, aes(Conections, GC3, z = AveExp))  +  theme_bw() + ggtitle(paste(normalization_scheme,   ": All points Stage II Expr. ",sep=""))   + geom_density_2d(bin=10)        + geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_II$Conections), linetype="dashed", color = "red")  +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_II$GC3), linetype="dashed", color = "red")   + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_II$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_II$GC3), linetype="dashed", color =  "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85)) 
    m12<-ggplot(Interactomes_GC3_T2_merged_Stage_III, aes(Conections, GC3, z = AveExp)) +  theme_bw() + ggtitle(paste(normalization_scheme,  ": All points Stage III Expr. ",sep=""))   + geom_density_2d(bin=10)       + geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_III$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_III$GC3), linetype="dashed", color = "red") + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_III$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_III$GC3), linetype="dashed", color = "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))   + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85)) 
      
    # FindClusters_resolution               
    png(filename=paste(output_dir,"countour_T2_Coonections_melt_",normalization_scheme,"_",TCGA_project,"_Stage_all_GC3.png",sep=""), width = 30, height = 50, res=600, units = "cm")  
            plot<-ggarrange(m4, m5, m6, m7, m8, m9, m10, m11,m12, m1,m2,m3, nrow = 4,ncol = 3, common.legend = TRUE, legend="bottom")
            print(annotate_figure(plot, top = text_grob(TCGA_project, face = "bold", size = 14)))
    dev.off()  

    h1<-ggplot(unique(Interactomes_GC3_T2_merged_all[,c("ENSEMBL","T2","GC3","Conections","Stages")]), aes(x=T2, color=Stages)) + geom_density() +   theme_bw()     + xlim(0, 50)
    h2<-ggplot(unique(Interactomes_GC3_T2_merged_all[,c("ENSEMBL","T2","GC3","Conections","Stages")]), aes(x=GC3, color=Stages)) +geom_density() +  theme_bw()    + xlim(0, 100)
    h3<-ggplot(unique(Interactomes_GC3_T2_merged_all[,c("ENSEMBL","T2","AveExp","Conections","Stages")]), aes(x=AveExp, color=Stages)) + geom_density() +  theme_bw()  + xlim(0, 10000)
    h4<-ggplot(unique(Interactomes_GC3_T2_merged_all[,c("ENSEMBL","T2","AveExp","Conections","Stages")]), aes(x=Conections, color=Stages)) + geom_density() +  theme_bw()  + xlim(0, 10000)

    i1<-ggplot(unique(Interactomes_GC3_T2_merged_all[,c("ENSEMBL","T2","GC3","Conections","Stages")]), aes(x=T2, color=Stages)) + geom_histogram(fill="white", alpha=0.5, position="identity", bins=20)  +   theme_bw()     + xlim(0, 50)
    i2<-ggplot(unique(Interactomes_GC3_T2_merged_all[,c("ENSEMBL","T2","GC3","Conections","Stages")]), aes(x=GC3, color=Stages)) + geom_histogram(fill="white", alpha=0.5, position="identity", bins=20)  +  theme_bw()    + xlim(0, 100)
    i3<-ggplot(unique(Interactomes_GC3_T2_merged_all[,c("ENSEMBL","T2","AveExp","Conections","Stages")]), aes(x=AveExp, color=Stages)) + geom_histogram(fill="white", alpha=0.5, position="identity", bins=20)  +  theme_bw()  + xlim(0, 10000)
    i4<-ggplot(unique(Interactomes_GC3_T2_merged_all[,c("ENSEMBL","T2","AveExp","Conections","Stages")]), aes(x=Conections, color=Stages)) + geom_histogram(fill="white", alpha=0.5, position="identity", bins=20)  +  theme_bw()  + xlim(0, 10000)    
    
    # FindClusters_resolution               
    png(filename=paste(output_dir,"countour_T2_Coonections_melt_",normalization_scheme,"_",TCGA_project,"_histogram.png",sep=""), width = 30, height = 30, res=600, units = "cm")            
          plot<-ggarrange(h1, h2, h3,h4,i1, i2, i3, i4, nrow = 2, ncol=4, common.legend = TRUE, legend="bottom") 
          annotate_figure(plot, top = text_grob(paste(TCGA_project,normalization_scheme,sep=" "), face = "bold", size = 14))  
    dev.off()    
}







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

    unique_stage_I  =intersect(setdiff(genes_stages_I, c(genes_stages_II,genes_stages_III)),genes_stages_I)
    unique_stage_II =intersect(setdiff(genes_stages_II, c(genes_stages_I,genes_stages_III)),genes_stages_II)
    unique_stage_III=intersect(setdiff(genes_stages_III, c(genes_stages_I,genes_stages_II)),genes_stages_III)

    Interactomes_GC3_T2_merged_Stage_I<-na.omit(Interactomes_GC3_T2_merged_Stage_I[unique_stage_I,])
    Interactomes_GC3_T2_merged_Stage_II<-na.omit(Interactomes_GC3_T2_merged_Stage_II[unique_stage_II,])
    Interactomes_GC3_T2_merged_Stage_III<-na.omit(Interactomes_GC3_T2_merged_Stage_III[unique_stage_III,])
    
    
    # Merge Stages
    Interactomes_GC3_T2_merged_all<-na.omit(rbind(Interactomes_GC3_T2_merged_Stage_I,Interactomes_GC3_T2_merged_Stage_II, Interactomes_GC3_T2_merged_Stage_III))  
    #########################################################################################################################################    
    # Merge dataset GC3_T2 with Expression Stages I, II and III
    Interactomes_GC3_T2_expr_stageI_unique<-cbind(Interactomes_GC3_T2_merged[unique_stage_I,],expression_table_normalized[unique_stage_I,list_of_comparisson[["sample_stage_I"]]])
    Interactomes_GC3_T2_expr_stageII_unique<-cbind(Interactomes_GC3_T2_merged[unique_stage_II,],expression_table_normalized[unique_stage_II,list_of_comparisson[["sample_stage_II"]]])
    Interactomes_GC3_T2_expr_stageIII_unique<-cbind(Interactomes_GC3_T2_merged[unique_stage_III,],expression_table_normalized[unique_stage_III,list_of_comparisson[["sample_stage_III"]]])
    
    # Melt data.frames
    Interactomes_GC3_T2_expr_stageI_unique   <-na.omit(melt(Interactomes_GC3_T2_expr_stageI_unique, id = c("T2", "GC3", "Conections", "ENSEMBL")))
    Interactomes_GC3_T2_expr_stageII_unique  <-na.omit(melt(Interactomes_GC3_T2_expr_stageII_unique, id = c("T2", "GC3", "Conections", "ENSEMBL")))
    Interactomes_GC3_T2_expr_stageIII_unique <-na.omit(melt(Interactomes_GC3_T2_expr_stageIII_unique, id = c("T2", "GC3", "Conections", "ENSEMBL")))

    Interactomes_GC3_T2_expr_stageI_unique$Stages<-"Stage I"
    Interactomes_GC3_T2_expr_stageII_unique$Stages<-"Stage II"
    Interactomes_GC3_T2_expr_stageIII_unique$Stages<-"Stage III"    
    
    # ggplot
    t1<-ggplot(Interactomes_GC3_T2_expr_stageI_unique, aes(Conections, T2, z = value))   + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,    ": All points only Stage I genes Expr. ",sep=""))+  geom_vline(xintercept=mean(Interactomes_GC3_T2_expr_stageI_unique$T2), linetype="dashed", color = "red")   +  geom_hline(yintercept=mean(Interactomes_GC3_T2_expr_stageI_unique$T2), linetype="dashed", color = "red")     + geom_vline(xintercept=median(Interactomes_GC3_T2_expr_stageI_unique$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_expr_stageI_unique$T2), linetype="dashed", color =   "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))             + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85)) 
    t2<-ggplot(Interactomes_GC3_T2_expr_stageII_unique, aes(Conections, T2, z = value))   + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,    ": All points only Stage II genes Expr. ",sep=""))+  geom_vline(xintercept=mean(Interactomes_GC3_T2_expr_stageII_unique$T2), linetype="dashed", color = "red")   +  geom_hline(yintercept=mean(Interactomes_GC3_T2_expr_stageII_unique$T2), linetype="dashed", color = "red")     + geom_vline(xintercept=median(Interactomes_GC3_T2_expr_stageII_unique$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_expr_stageII_unique$T2), linetype="dashed", color =   "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))       + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85)) 
    t3<-ggplot(Interactomes_GC3_T2_expr_stageIII_unique, aes(Conections, T2, z = value))   + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,    ": All points only Stage III genes Expr. ",sep=""))+  geom_vline(xintercept=mean(Interactomes_GC3_T2_expr_stageIII_unique$T2), linetype="dashed", color = "red")   +  geom_hline(yintercept=mean(Interactomes_GC3_T2_expr_stageIII_unique$T2), linetype="dashed", color = "red")     + geom_vline(xintercept=median(Interactomes_GC3_T2_expr_stageIII_unique$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_expr_stageIII_unique$T2), linetype="dashed", color =   "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE)) + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85)) 

    # Unique values
    Interactomes_GC3_T2_unique_all<-rbind(Interactomes_GC3_T2_expr_stageI_unique,Interactomes_GC3_T2_expr_stageII_unique, Interactomes_GC3_T2_expr_stageIII_unique)
    #########################################################################################################################################
    m1<-ggplot(Interactomes_GC3_T2_expr_stageI_unique, aes(Conections, GC3, z = value))   + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))+  geom_vline(xintercept=mean(Interactomes_GC3_T2_expr_stageI_unique$Conections), linetype="dashed", color = "red")   +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_I$GC3), linetype="dashed", color = "red")     + geom_vline(xintercept=median(Interactomes_GC3_T2_expr_stageI_unique$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_expr_stageI_unique$GC3), linetype="dashed", color =   "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85)) 
    m2<-ggplot(Interactomes_GC3_T2_expr_stageII_unique, aes(Conections, GC3, z = value))  + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,   ": All points Stage II Expr. ",sep=""))+  geom_vline(xintercept=mean(Interactomes_GC3_T2_expr_stageII_unique$Conections), linetype="dashed", color = "red")  +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_II$GC3), linetype="dashed", color = "red")   + geom_vline(xintercept=median(Interactomes_GC3_T2_expr_stageII_unique$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_expr_stageII_unique$GC3), linetype="dashed", color =  "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85)) 
    m3<-ggplot(Interactomes_GC3_T2_expr_stageIII_unique, aes(Conections, GC3, z = value)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,  ": All points Stage III Expr. ",sep=""))+  geom_vline(xintercept=mean(Interactomes_GC3_T2_expr_stageIII_unique$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_III$GC3), linetype="dashed", color = "red") + geom_vline(xintercept=median(Interactomes_GC3_T2_expr_stageIII_unique$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_expr_stageIII_unique$GC3), linetype="dashed", color = "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85)) 
  
    m4 <- ggplot(Interactomes_GC3_T2_unique_all, aes(x=Stages, y=GC3)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()           +  geom_hline(yintercept=median(Interactomes_GC3_T2_unique_all$GC3), linetype="dashed", color = "red")    + ggtitle(paste("T2 All points: ", normalization_scheme,sep=""))                   +  stat_compare_means(comparisons = my_comparisons, method = "t.test") 
    m5 <- ggplot(Interactomes_GC3_T2_unique_all, aes(x=Stages, y=Conections)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()   +  geom_hline(yintercept=median(Interactomes_GC3_T2_unique_all$Conections), linetype="dashed", color = "red")    + ggtitle(paste("Connectivity All points: ", normalization_scheme,sep="")) +  stat_compare_means(comparisons = my_comparisons, method = "t.test") 
    m6 <- ggplot(Interactomes_GC3_T2_unique_all, aes(x=Stages, y=value)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()         +      geom_hline(yintercept=median(Interactomes_GC3_T2_unique_all$Expr), linetype="dashed", color = "red")    + ggtitle(paste("Expr. All points: ", normalization_scheme,sep=""))          +  stat_compare_means(comparisons = my_comparisons, method = "t.test")       
  
    m7<-ggplot(Interactomes_GC3_T2_expr_stageI_unique, aes(Conections, GC3, z = value))    + geom_point(aes(colour=value))   +  theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))+  xlim(0, 50)       + geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_I$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_I$GC3), linetype="dashed", color = "red")     + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_I$Conections), linetype="dashed", color = "yellow") +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_I$GC3), linetype="dashed", color =   "yellow")   + ylim(25,85) + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85)) 
    m8<-ggplot(Interactomes_GC3_T2_expr_stageII_unique, aes(Conections, GC3, z = value))   + geom_point(aes(colour=value))  +  theme_bw() + ggtitle(paste(normalization_scheme,   ": All points Stage II Expr. ",sep=""))+ xlim(0, 50) + geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_II$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_II$GC3), linetype="dashed", color = "red")   + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_II$Conections), linetype="dashed", color = "yellow") +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_II$GC3), linetype="dashed", color =  "yellow")   + ylim(25,85) + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))         + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85)) 
    m9<-ggplot(Interactomes_GC3_T2_expr_stageIII_unique, aes(Conections, GC3, z = value))  + geom_point(aes(colour=value)) + theme_bw() + ggtitle(paste(normalization_scheme,  ": All points Stage III Expr. ",sep=""))+ xlim(0, 50) +  geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_III$Conections), linetype="dashed", color = "red")  +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_III$GC3), linetype="dashed", color = "red") + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_III$Conections), linetype="dashed", color = "yellow") +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_III$GC3), linetype="dashed", color = "yellow")   + ylim(25,85) +  guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))       + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85)) 
  
    m10<-ggplot(Interactomes_GC3_T2_expr_stageI_unique, aes(Conections, GC3, z = value))   +  theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))   + geom_density_2d(bin=10)         + geom_vline(xintercept=mean(Interactomes_GC3_T2_expr_stageI_unique$Conections), linetype="dashed", color = "red")   +  geom_hline(yintercept=mean(Interactomes_GC3_T2_expr_stageI_unique$GC3), linetype="dashed", color = "red")     + geom_vline(xintercept=median(Interactomes_GC3_T2_expr_stageI_unique$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_expr_stageI_unique$GC3), linetype="dashed", color =   "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE)) + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85)) 
    m11<-ggplot(Interactomes_GC3_T2_expr_stageI_unique, aes(Conections, GC3, z = value))   +  theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))   + geom_density_2d(bin=10)         + geom_vline(xintercept=mean(Interactomes_GC3_T2_expr_stageI_unique$Conections), linetype="dashed", color = "red")   +  geom_hline(yintercept=mean(Interactomes_GC3_T2_expr_stageI_unique$GC3), linetype="dashed", color = "red")     + geom_vline(xintercept=median(Interactomes_GC3_T2_expr_stageI_unique$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_expr_stageI_unique$GC3), linetype="dashed", color =   "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE)) + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85)) 
    m12<-ggplot(Interactomes_GC3_T2_expr_stageIII_unique, aes(Conections, GC3, z = value)) +  theme_bw() + ggtitle(paste(normalization_scheme,  ": All points Stage III Expr. ",sep=""))   + geom_density_2d(bin=10)       + geom_vline(xintercept=mean(Interactomes_GC3_T2_expr_stageIII_unique$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(Interactomes_GC3_T2_expr_stageIII_unique$GC3), linetype="dashed", color = "red") + geom_vline(xintercept=median(Interactomes_GC3_T2_expr_stageIII_unique$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_III$GC3), linetype="dashed", color = "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))       + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85)) 
      
    # FindClusters_resolution               
    png(filename=depravação moral, width = 30, height = 50, res=600, units = "cm")  
            plot<-ggarrange(m4, m5, m6, m7, m8, m9, m10, m11,m12, m1,m2,m3, nrow = 4,ncol = 3, common.legend = TRUE, legend="bottom")
            print(annotate_figure(plot, top = text_grob(TCGA_project, face = "bold", size = 14)))
    dev.off()  

    #########################################################################################################################################
    m1<-ggplot(Interactomes_GC3_T2_expr_stageI_unique, aes(Conections, T2, z = value))   + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))+  geom_vline(xintercept=mean(Interactomes_GC3_T2_expr_stageI_unique$Conections), linetype="dashed", color = "red")   +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_I$T2), linetype="dashed", color = "red")     + geom_vline(xintercept=median(Interactomes_GC3_T2_expr_stageI_unique$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_expr_stageI_unique$T2), linetype="dashed", color =   "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))
    m2<-ggplot(Interactomes_GC3_T2_expr_stageII_unique, aes(Conections, T2, z = value))  + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,   ": All points Stage II Expr. ",sep=""))+  geom_vline(xintercept=mean(Interactomes_GC3_T2_expr_stageII_unique$Conections), linetype="dashed", color = "red")  +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_II$T2), linetype="dashed", color = "red")   + geom_vline(xintercept=median(Interactomes_GC3_T2_expr_stageII_unique$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_expr_stageII_unique$T2), linetype="dashed", color =  "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))
    m3<-ggplot(Interactomes_GC3_T2_expr_stageIII_unique, aes(Conections, T2, z = value)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,  ": All points Stage III Expr. ",sep=""))+  geom_vline(xintercept=mean(Interactomes_GC3_T2_expr_stageIII_unique$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_III$T2), linetype="dashed", color = "red") + geom_vline(xintercept=median(Interactomes_GC3_T2_expr_stageIII_unique$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_expr_stageIII_unique$T2), linetype="dashed", color = "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))
  
    m4 <- ggplot(Interactomes_GC3_T2_unique_all, aes(x=Stages, y=GC3)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()           +  geom_hline(yintercept=median(Interactomes_GC3_T2_unique_all$T2), linetype="dashed", color = "red")    + ggtitle(paste("T2 All points: ", normalization_scheme,sep=""))                   +  stat_compare_means(comparisons = my_comparisons, method = "t.test") 
    m5 <- ggplot(Interactomes_GC3_T2_unique_all, aes(x=Stages, y=Conections)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()   +  geom_hline(yintercept=median(Interactomes_GC3_T2_unique_all$Conections), linetype="dashed", color = "red")    + ggtitle(paste("Connectivity All points: ", normalization_scheme,sep="")) +  stat_compare_means(comparisons = my_comparisons, method = "t.test") 
    m6 <- ggplot(Interactomes_GC3_T2_unique_all, aes(x=Stages, y=value)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()         +      geom_hline(yintercept=median(Interactomes_GC3_T2_unique_all$Expr), linetype="dashed", color = "red")    + ggtitle(paste("Expr. All points: ", normalization_scheme,sep=""))          +  stat_compare_means(comparisons = my_comparisons, method = "t.test")       
  
    m7<-ggplot(Interactomes_GC3_T2_expr_stageI_unique, aes(Conections, T2, z = value))    + geom_point(aes(colour=value))   +  theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))+  xlim(0, 50)       + geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_I$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_I$T2), linetype="dashed", color = "red")     + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_I$Conections), linetype="dashed", color = "yellow") +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_I$T2), linetype="dashed", color =   "yellow")   + ylim(10,40) + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))
    m8<-ggplot(Interactomes_GC3_T2_expr_stageII_unique, aes(Conections, T2, z = value))   + geom_point(aes(colour=value))  +  theme_bw() + ggtitle(paste(normalization_scheme,   ": All points Stage II Expr. ",sep=""))+ xlim(0, 50) + geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_II$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_II$T2), linetype="dashed", color = "red")   + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_II$Conections), linetype="dashed", color = "yellow") +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_II$T2), linetype="dashed", color =  "yellow")   + ylim(10,40) + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))         + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))
    m9<-ggplot(Interactomes_GC3_T2_expr_stageIII_unique, aes(Conections, T2, z = value))  + geom_point(aes(colour=value)) + theme_bw() + ggtitle(paste(normalization_scheme,  ": All points Stage III Expr. ",sep=""))+ xlim(0, 50) +  geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_III$Conections), linetype="dashed", color = "red")  +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_III$T2), linetype="dashed", color = "red") + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_III$Conections), linetype="dashed", color = "yellow") +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_III$T2), linetype="dashed", color = "yellow")   + ylim(10,40) +  guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))       + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))
  
    m10<-ggplot(Interactomes_GC3_T2_expr_stageI_unique, aes(Conections, T2, z = value))   +  theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))   + geom_density_2d(bin=10)         + geom_vline(xintercept=mean(Interactomes_GC3_T2_expr_stageI_unique$Conections), linetype="dashed", color = "red")   +  geom_hline(yintercept=mean(Interactomes_GC3_T2_expr_stageI_unique$T2), linetype="dashed", color = "red")     + geom_vline(xintercept=median(Interactomes_GC3_T2_expr_stageI_unique$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_expr_stageI_unique$T2), linetype="dashed", color =   "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE)) + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))
    m11<-ggplot(Interactomes_GC3_T2_expr_stageI_unique, aes(Conections, T2, z = value))   +  theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))   + geom_density_2d(bin=10)         + geom_vline(xintercept=mean(Interactomes_GC3_T2_expr_stageI_unique$Conections), linetype="dashed", color = "red")   +  geom_hline(yintercept=mean(Interactomes_GC3_T2_expr_stageI_unique$T2), linetype="dashed", color = "red")     + geom_vline(xintercept=median(Interactomes_GC3_T2_expr_stageI_unique$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_expr_stageI_unique$T2), linetype="dashed", color =   "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE)) + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))
    m12<-ggplot(Interactomes_GC3_T2_expr_stageIII_unique, aes(Conections, T2, z = value)) +  theme_bw() + ggtitle(paste(normalization_scheme,  ": All points Stage III Expr. ",sep=""))   + geom_density_2d(bin=10)       + geom_vline(xintercept=mean(Interactomes_GC3_T2_expr_stageIII_unique$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(Interactomes_GC3_T2_expr_stageIII_unique$T2), linetype="dashed", color = "red") + geom_vline(xintercept=median(Interactomes_GC3_T2_expr_stageIII_unique$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_III$T2), linetype="dashed", color = "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))       + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))
      
    # FindClusters_resolution               
    png(filename=paste(output_dir,"countour_T2_Coonections_melt_",normalization_scheme,"_",TCGA_project,"_Stage_unique_T2.png",sep=""), width = 30, height = 50, res=600, units = "cm")  
            plot<-ggarrange(m4, m5, m6, m7, m8, m9, m10, m11,m12, m1,m2,m3, nrow = 4,ncol = 3, common.legend = TRUE, legend="bottom")
            print(annotate_figure(plot, top = text_grob(TCGA_project, face = "bold", size = 14)))
    dev.off()    
    #########################################################################################################################################
    h1<-ggplot(unique(Interactomes_GC3_T2_unique_all[,c("ENSEMBL","T2","GC3","Conections","Stages")]), aes(x=T2, color=Stages)) + geom_density() +   theme_bw()     + xlim(0, 50)
    h2<-ggplot(unique(Interactomes_GC3_T2_unique_all[,c("ENSEMBL","T2","GC3","Conections","Stages")]), aes(x=GC3, color=Stages)) +geom_density() +  theme_bw()    + xlim(0, 100)
    h3<-ggplot(unique(Interactomes_GC3_T2_unique_all[,c("ENSEMBL","T2","value","Conections","Stages")]), aes(x=value, color=Stages)) + geom_density() +  theme_bw()  + xlim(0, 10000)
    h4<-ggplot(unique(Interactomes_GC3_T2_unique_all[,c("ENSEMBL","T2","value","Conections","Stages")]), aes(x=Conections, color=Stages)) + geom_density() +  theme_bw()  + xlim(0, 10000)

    i1<-ggplot(unique(Interactomes_GC3_T2_unique_all[,c("ENSEMBL","T2","GC3","Conections","Stages")]), aes(x=T2, color=Stages)) + geom_histogram(fill="white", alpha=0.5, position="identity", bins=20)  +   theme_bw()     + xlim(0, 50)
    i2<-ggplot(unique(Interactomes_GC3_T2_unique_all[,c("ENSEMBL","T2","GC3","Conections","Stages")]), aes(x=GC3, color=Stages)) + geom_histogram(fill="white", alpha=0.5, position="identity", bins=20)  +  theme_bw()    + xlim(0, 100)
    i3<-ggplot(unique(Interactomes_GC3_T2_unique_all[,c("ENSEMBL","T2","value","Conections","Stages")]), aes(x=value, color=Stages)) + geom_histogram(fill="white", alpha=0.5, position="identity", bins=20)  +  theme_bw()  + xlim(0, 10000)
    i4<-ggplot(unique(Interactomes_GC3_T2_unique_all[,c("ENSEMBL","T2","value","Conections","Stages")]), aes(x=Conections, color=Stages)) + geom_histogram(fill="white", alpha=0.5, position="identity", bins=20)  +  theme_bw()  + xlim(0, 10000)    
    
    # FindClusters_resolution               
    png(filename=paste(output_dir,"countour_T2_Coonections_melt_",normalization_scheme,"_",TCGA_project,"_unique_histogram.png",sep=""), width = 30, height = 30, res=600, units = "cm")            
          plot<-ggarrange(h1, h2, h3,h4,i1, i2, i3, i4, nrow = 2, ncol=4, common.legend = TRUE, legend="bottom") 
          annotate_figure(plot, top = text_grob(paste(TCGA_project,normalization_scheme,sep=" "), face = "bold", size = 14))  
    dev.off()    
}
