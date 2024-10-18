#######################################################################################################################
#source("/home/felipe/Documents/github/CancerStaging/CancerStaging_SetupAllParamters.R")                               #
#source("/home/felipe/Documents/github/CancerStaging/CancerStaging_LoadRPackages.R")                                   #
#######################################################################################################################
# saveRDS                                                                                                             #
#saveRDS(object = df_reads_count_all_projects, file = paste(output_dir,"df_reads_count_all_projects.rds",sep=""))      #
#saveRDS(object = geneLength_ENTREZID_ENSEMBL, file = paste(output_dir,"geneLength_ENTREZID_ENSEMBL.rds",sep=""))      #
#saveRDS(object = Interactomes_GC3_T2_merged, file = paste(output_dir,"Interactomes_GC3_T2_merged.rds",sep=""))        #
                                                                                                                      #
# Restore the object                                                                                                  #
#geneLength_ENTREZID_ENSEMBL<-readRDS(file = paste(output_dir,"geneLength_ENTREZID_ENSEMBL.rds",sep=""))               #
#Interactomes_GC3_T2_merged <-readRDS(file = paste(output_dir,"Interactomes_GC3_T2_merged.rds",sep=""))                #
#df_reads_count_all_projects<-readRDS(file = paste(output_dir,"df_reads_count_all_projects.rds",sep=""))               ###########################
##merged_data_patient_info   <-read.table(file = "/home/felipe/Documents/Cancer_staging/merged_data_patient_info.tsv", sep = '\t', header = TRUE) #
#################################################################################################################################################
# Repeat for evey normalization scheme                                                                                                          #
################################################################################################################################################# 
# Construnction of 3d coordinates.
# X=T2           : Thymine composition in second codon position (T2)
# Y=connectivity : interactome_data                 # Carels checked    .
# Z=expression   : reads_count_all_projects         # Use the fpkm
# Reading the contents of TSV file using read_tsv() method
Interactomes_GC3_T2_file <-"/home/felipe/Documents/github/CancerStaging/Interactomes_GC3_T2.csv"

# Read the Interactomes_GC3_T2_data
# First  collumn : UniprotKB
# Second collumn : SYMBOL
# Third  collumn : GC3
# Fourth collumn : T2
# Fifth  collumn : Conections
Interactomes_GC3_T2_data <-read.table(file = Interactomes_GC3_T2_file, sep = '\t', header = TRUE,fill=TRUE) 

# Set the second collum from GeneSymbol to SYMBOL
colnames(Interactomes_GC3_T2_data)[2]<-"SYMBOL"

# Merge Interactomes_GC3_T2 and geneLength_ENTREZID_ENSEMBL by gene SYMBOL.
# geneLength_ENTREZID_ENSEMBL contains the following collumn:
# ENTREZID   entrezid identifier
# geneLength genelenght calculated by Carels
# SYMBOL     gene symbol
# ENSEMBL    ENSEMBL symbol
Interactomes_GC3_T2_merged<-merge(geneLength_ENTREZID_ENSEMBL,Interactomes_GC3_T2_data,by="SYMBOL")

# Transform the Interactomes_GC3_T2_merged to as.data.table
Interactomes_GC3_T2_merged<-as.data.table(Interactomes_GC3_T2_merged)

# Set ensembl ids
ENSEMBL_ids<-unique(Interactomes_GC3_T2_merged$ENSEMBL)

# Selecte collumns to be extracted, Interactomes_GC3_T2_merged
# "T2","GC3","Conections" and "ENSEMBL"
Interactomes_GC3_T2_merged<-data.frame(Interactomes_GC3_T2_merged[which(Interactomes_GC3_T2_merged$ENSEMBL %in% ENSEMBL_ids),c("T2","GC3","Conections","ENSEMBL")])

# In the file "Interactomes_GC3_T2" each genes can have multiple entries with multiple values for the variable connections.
# Only the occurance with greatest number of Conections will be used.
Interactomes_GC3_T2_df<-data.frame(T2=c(),GC3=c(),Conections=c(),ENSEMBL=c())

# Interactomes_GC3_T2_merged
for (gene_ENSEMBL in unique(Interactomes_GC3_T2_merged$ENSEMBL))
{
  # Take all entries of this 
  Interactomes_GC3_T2_all<-Interactomes_GC3_T2_merged[which(Interactomes_GC3_T2_merged$ENSEMBL==gene_ENSEMBL),]

  # Interactomes_GC3_T2_df
  Interactomes_GC3_T2_df<-rbind(Interactomes_GC3_T2_df,Interactomes_GC3_T2_all[which.max(Interactomes_GC3_T2_all$Conections),])
}
# Replace the data.frames
Interactomes_GC3_T2_merged<-Interactomes_GC3_T2_df

# Set rownames
rownames(Interactomes_GC3_T2_merged)<-Interactomes_GC3_T2_merged$ENSEMBL

# Filter by T2
Interactomes_GC3_T2_merged<-Interactomes_GC3_T2_merged[Interactomes_GC3_T2_merged$T2 <=85,]
Interactomes_GC3_T2_merged<-Interactomes_GC3_T2_merged[Interactomes_GC3_T2_merged$GC3 <=100,]
#limit_expr=10000
#limit_expr=100

# Filter by T2
Interactomes_GC3_T2_merged<-na.omit(Interactomes_GC3_T2_merged)
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
    
    # merged_expression_interactomes
    merged_expression_table_normalized_stage_I <-cbind(expression_table_normalized_stage_I[ENSEMBL_ids_stage_I,],Interactomes_GC3_T2_merged_Stage_I[Interactomes_GC3_T2_merged_Stage_I$ENSEMBL %in% ENSEMBL_ids_stage_I,c("T2","GC3","Conections","ENSEMBL")])
    merged_expression_table_normalized_stage_II<-cbind(expression_table_normalized_stage_II[ENSEMBL_ids_stage_II,],Interactomes_GC3_T2_merged_Stage_II[Interactomes_GC3_T2_merged_Stage_II$ENSEMBL %in% ENSEMBL_ids_stage_II,c("T2","GC3","Conections","ENSEMBL")])  
    merged_expression_table_normalized_stage_III<-cbind(expression_table_normalized_stage_III[ENSEMBL_ids_stage_III,],Interactomes_GC3_T2_merged_Stage_III[Interactomes_GC3_T2_merged_Stage_III$ENSEMBL %in% ENSEMBL_ids_stage_III,c("T2","GC3","Conections","ENSEMBL")])  

    # Rempove NA lines
    merged_expression_table_normalized_stage_I<-merged_expression_table_normalized_stage_I[complete.cases(merged_expression_table_normalized_stage_I), ]    
    merged_expression_table_normalized_stage_II<-merged_expression_table_normalized_stage_II[complete.cases(merged_expression_table_normalized_stage_II), ]    
    merged_expression_table_normalized_stage_III<-merged_expression_table_normalized_stage_III[complete.cases(merged_expression_table_normalized_stage_III), ]    
    
    # Melt data.frame 
    merged_expression_table_normalized_stage_I <- melt(data.frame(merged_expression_table_normalized_stage_I), id=c("T2","GC3","Conections","ENSEMBL"))
    merged_expression_table_normalized_stage_II <- melt(data.frame(merged_expression_table_normalized_stage_II), id=c("T2","GC3","Conections","ENSEMBL"))
    merged_expression_table_normalized_stage_III <- melt(data.frame(merged_expression_table_normalized_stage_III), id=c("T2","GC3","Conections","ENSEMBL"))

    # Change variables
    merged_expression_table_normalized_stage_I$variable<-gsub('\\.', '-', merged_expression_table_normalized_stage_I$variable)  
    merged_expression_table_normalized_stage_II$variable<-gsub('\\.', '-', merged_expression_table_normalized_stage_II$variable)  
    merged_expression_table_normalized_stage_III$variable<-gsub('\\.', '-', merged_expression_table_normalized_stage_III$variable)  
    
    # Set colnames
    colnames(merged_expression_table_normalized_stage_I)[6]<-normalization_scheme
    colnames(merged_expression_table_normalized_stage_II)[6]<-normalization_scheme
    colnames(merged_expression_table_normalized_stage_III)[6]<-normalization_scheme

    #########################################################################################################################################        
    selected_genes_Stage_I_gene      <- genes_stages_I
    selected_genes_Stage_II_gene     <- genes_stages_II
    selected_genes_Stage_III_gene    <- genes_stages_III
    #######################################################################################################################################                                                                                                                                     #    
    unique_stage_I  =intersect(setdiff(selected_genes_Stage_I_gene, c(selected_genes_Stage_II_gene,selected_genes_Stage_III_gene)),selected_genes_Stage_I_gene)
    unique_stage_II =intersect(setdiff(selected_genes_Stage_II_gene, c(selected_genes_Stage_I_gene,selected_genes_Stage_III_gene)),selected_genes_Stage_II_gene)
    unique_stage_III=intersect(setdiff(selected_genes_Stage_III_gene, c(selected_genes_Stage_I_gene,selected_genes_Stage_II_gene)),selected_genes_Stage_III_gene)
    #######################################################################################################################################                                                                                                                                     #  
  
    # Melt data.frame 
    merged_expression_table_normalized_all_stages <- rbind(merged_expression_table_normalized_stage_I,merged_expression_table_normalized_stage_II,merged_expression_table_normalized_stage_III)
  
    # Set the field stages
    Interactomes_GC3_T2_merged$Stages                   <-"overlapping"
    merged_expression_table_normalized_all_stages$Stages<-"overlapping"
  
    # Set stages 
    Interactomes_GC3_T2_merged[unique_stage_I,"Stages"]<-"Stage I"
    Interactomes_GC3_T2_merged[unique_stage_II,"Stages"]<-"Stage II"
    Interactomes_GC3_T2_merged[unique_stage_III,"Stages"]<-"Stage III"

    merged_expression_table_normalized_all_stages[merged_expression_table_normalized_all_stages$ENSEMBL %in% unique_stage_I,"Stages"]<-"Stage I"
    merged_expression_table_normalized_all_stages[merged_expression_table_normalized_all_stages$ENSEMBL %in% unique_stage_II,"Stages"]<-"Stage II"
    merged_expression_table_normalized_all_stages[merged_expression_table_normalized_all_stages$ENSEMBL %in% unique_stage_III,"Stages"]<-"Stage III"    
    #########################################################################################################################################    
  
    ###################################################################################################################################################################
    colnames(merged_expression_table_normalized_stage_I)[6]<-"Expr"
    colnames(merged_expression_table_normalized_stage_II)[6]<-"Expr"
    colnames(merged_expression_table_normalized_stage_III)[6]<-"Expr"
    colnames(merged_expression_table_normalized_all_stages)[6]<-"Expr"    
    #########################################################################################################################################
    # Visualize: Specify the comparisons you want
    my_comparisons <- list( c("Stage I", "Stage II"), c("Stage I", "Stage III"), c("Stage II", "Stage III"), c("Stage I", "overlapping"),c("Stage II", "overlapping"),c("Stage III", "overlapping"))
    #########################################################################################################################################
    # Three countour plots will be created
    # One with the average expression
    # Second with the expression per patient
    # Third with the z-score 
    #########################################################################################################################################
    # Filter up Average expression greater than zero
    merged_expression_table_normalized_stage_I<-merged_expression_table_normalized_all_stages[merged_expression_table_normalized_all_stages$Stages=="Stage I",]
    merged_expression_table_normalized_stage_II<-merged_expression_table_normalized_all_stages[merged_expression_table_normalized_all_stages$Stages=="Stage II",]
    merged_expression_table_normalized_stage_III<-merged_expression_table_normalized_all_stages[merged_expression_table_normalized_all_stages$Stages=="Stage III",]
    merged_expression_table_normalized_stage_overlapping<-merged_expression_table_normalized_all_stages[merged_expression_table_normalized_all_stages$Stages=="overlapping",]
 
    #merged_expression_table_normalized_stage_I  <-merged_expression_table_normalized_stage_I[merged_expression_table_normalized_stage_I$Expr>0 & merged_expression_table_normalized_stage_I$Expr<limit_expr,]
    #merged_expression_table_normalized_stage_II <-merged_expression_table_normalized_stage_II[merged_expression_table_normalized_stage_II$Expr>0 & merged_expression_table_normalized_stage_II$Expr<limit_expr,]
    #merged_expression_table_normalized_stage_III <-merged_expression_table_normalized_stage_III[merged_expression_table_normalized_stage_III$Expr>0  & merged_expression_table_normalized_stage_III$Expr<limit_expr,]
    #merged_expression_table_normalized_all_stages <-merged_expression_table_normalized_all_stages[merged_expression_table_normalized_all_stages$Expr>0  & merged_expression_table_normalized_all_stages$Expr<limit_expr,]    

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
    merged_expression_table_normalized_all_stages<-rbind(merged_expression_table_normalized_stage_I,merged_expression_table_normalized_stage_II,merged_expression_table_normalized_stage_III, merged_expression_table_normalized_stage_overlapping)  
  
    m1<-ggplot(merged_expression_table_normalized_stage_I, aes(Conections, T2, z = Expr))   + geom_density_2d_filled(alpha = 0.5) + theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))+ geom_contour()         + geom_vline(xintercept=mean(merged_expression_table_normalized_stage_I$Conections), linetype="dashed", color = "red")   +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_I$T2), linetype="dashed", color = "red")     + geom_vline(xintercept=median(merged_expression_table_normalized_stage_I$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(merged_expression_table_normalized_stage_I$T2), linetype="dashed", color =   "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE)) + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(0, 40, by = 1), breaks = seq(0, 40, by = 10), limits = c(10, 40))
    m2<-ggplot(merged_expression_table_normalized_stage_II, aes(Conections, T2, z = Expr))  + geom_density_2d_filled(alpha = 0.5) + theme_bw() + ggtitle(paste(normalization_scheme,   ": All points Stage II Expr. ",sep=""))+ geom_contour()        + geom_vline(xintercept=mean(merged_expression_table_normalized_stage_II$Conections), linetype="dashed", color = "red")  +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_II$T2), linetype="dashed", color = "red")   + geom_vline(xintercept=median(merged_expression_table_normalized_stage_II$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(merged_expression_table_normalized_stage_II$T2), linetype="dashed", color =  "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(0, 40, by = 1), breaks = seq(0, 40, by = 10), limits = c(10, 40))
    m3<-ggplot(merged_expression_table_normalized_stage_III, aes(Conections, T2, z = Expr)) + geom_density_2d_filled(alpha = 0.5) + theme_bw() + ggtitle(paste(normalization_scheme,  ": All points Stage III Expr. ",sep=""))+ geom_contour()       + geom_vline(xintercept=mean(merged_expression_table_normalized_stage_III$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_III$T2), linetype="dashed", color = "red") + geom_vline(xintercept=median(merged_expression_table_normalized_stage_III$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(merged_expression_table_normalized_stage_III$T2), linetype="dashed", color = "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(0, 40, by = 1), breaks = seq(0, 40, by = 10), limits = c(10, 40))
    
    m4 <- ggplot(merged_expression_table_normalized_all_stages, aes(x=Stages, y=T2)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()           +  geom_hline(yintercept=median(merged_expression_table_normalized_all_stages$T2), linetype="dashed", color = "red")    + ggtitle(paste("T2 All points: ", normalization_scheme,sep=""))                   +  stat_compare_means(comparisons = my_comparisons, method = "t.test") 
    m5 <- ggplot(merged_expression_table_normalized_all_stages, aes(x=Stages, y=Conections)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()   +  geom_hline(yintercept=median(merged_expression_table_normalized_all_stages$Conections), linetype="dashed", color = "red")    + ggtitle(paste("Connectivity All points: ", normalization_scheme,sep="")) +  stat_compare_means(comparisons = my_comparisons, method = "t.test") 
    m6 <- ggplot(merged_expression_table_normalized_all_stages, aes(x=Stages, y=Expr)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()         +      geom_hline(yintercept=median(merged_expression_table_normalized_all_stages$Expr), linetype="dashed", color = "red")    + ggtitle(paste("Expr. All points: ", normalization_scheme,sep=""))          +  stat_compare_means(comparisons = my_comparisons, method = "t.test") 

    m7<-ggplot(merged_expression_table_normalized_stage_I, aes(Conections, T2, z = Expr))  + geom_point(aes(colour=Expr))   +  theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))+ geom_contour()  + xlim(0, 50)       + geom_vline(xintercept=mean(merged_expression_table_normalized_stage_I$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_I$T2), linetype="dashed", color = "red")     + geom_vline(xintercept=median(merged_expression_table_normalized_stage_I$Conections), linetype="dashed", color = "yellow") +  geom_hline(yintercept=median(merged_expression_table_normalized_stage_I$T2), linetype="dashed", color =   "yellow")   + ylim(10,40) + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(0, 40, by = 1), breaks = seq(0, 40, by = 10), limits = c(10, 40))
    m8<-ggplot(merged_expression_table_normalized_stage_II, aes(Conections, T2, z = Expr))  + geom_point(aes(colour=Expr))  +  theme_bw() + ggtitle(paste(normalization_scheme,   ": All points Stage II Expr. ",sep=""))+ geom_contour()  + xlim(0, 50) + geom_vline(xintercept=mean(merged_expression_table_normalized_stage_II$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_II$T2), linetype="dashed", color = "red")   + geom_vline(xintercept=median(merged_expression_table_normalized_stage_II$Conections), linetype="dashed", color = "yellow") +  geom_hline(yintercept=median(merged_expression_table_normalized_stage_II$T2), linetype="dashed", color =  "yellow")   + ylim(10,40) + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)       + scale_y_continuous(minor_breaks = seq(0, 40, by = 1), breaks = seq(0, 40, by = 10), limits = c(10, 40))
    m9<-ggplot(merged_expression_table_normalized_stage_III, aes(Conections, T2, z = Expr))  + geom_point(aes(colour=Expr)) + theme_bw() + ggtitle(paste(normalization_scheme,  ": All points Stage III Expr. ",sep=""))+ geom_contour()  + xlim(0, 50) +  geom_vline(xintercept=mean(merged_expression_table_normalized_stage_III$Conections), linetype="dashed", color = "red")  +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_III$T2), linetype="dashed", color = "red") + geom_vline(xintercept=median(merged_expression_table_normalized_stage_III$Conections), linetype="dashed", color = "yellow") +  geom_hline(yintercept=median(merged_expression_table_normalized_stage_III$T2), linetype="dashed", color = "yellow")   + ylim(10,40) +  guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE)) + ylim(10,40)  + xlim(0, 50)     + scale_y_continuous(minor_breaks = seq(0, 40, by = 1), breaks = seq(0, 40, by = 10), limits = c(10, 40))
  
    # FindClusters_resolution               
    png(filename=paste(output_dir,"countour_T2_Coonections_melt_",normalization_scheme,"_",TCGA_project,"_Stage_all_T2.png",sep=""), width = 30, height = 50, res=600, units = "cm")  
            plot<-ggarrange(m7, m8, m9, m1,m2,m3,m4,m5,m6, nrow = 3,ncol = 3, common.legend = TRUE, legend="bottom")
            annotate_figure(plot, top = text_grob(TCGA_project, face = "bold", size = 14))
    dev.off()  
  

    m1<-ggplot(merged_expression_table_normalized_stage_I, aes(Conections, GC3, z = Expr))   + geom_density_2d_filled(alpha = 0.5) + theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))+ geom_contour()         + geom_vline(xintercept=mean(merged_expression_table_normalized_stage_I$Conections), linetype="dashed", color = "red")   +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_I$GC3), linetype="dashed", color = "red")     + geom_vline(xintercept=median(merged_expression_table_normalized_stage_I$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(merged_expression_table_normalized_stage_I$GC3), linetype="dashed", color =   "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))   + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 5), limits = c(25, 85))
    m2<-ggplot(merged_expression_table_normalized_stage_II, aes(Conections, GC3, z = Expr))  + geom_density_2d_filled(alpha = 0.5) + theme_bw() + ggtitle(paste(normalization_scheme,   ": All points Stage II Expr. ",sep=""))+ geom_contour()        + geom_vline(xintercept=mean(merged_expression_table_normalized_stage_II$Conections), linetype="dashed", color = "red")  +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_II$GC3), linetype="dashed", color = "red")   + geom_vline(xintercept=median(merged_expression_table_normalized_stage_II$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(merged_expression_table_normalized_stage_II$GC3), linetype="dashed", color =  "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))   + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 5), limits = c(25, 85))
    m3<-ggplot(merged_expression_table_normalized_stage_III, aes(Conections, GC3, z = Expr))  + geom_density_2d_filled(alpha = 0.5) + theme_bw() + ggtitle(paste(normalization_scheme,   ": All points Stage II Expr. ",sep=""))+ geom_contour()        + geom_vline(xintercept=mean(merged_expression_table_normalized_stage_III$Conections), linetype="dashed", color = "red")  +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_III$GC3), linetype="dashed", color = "red")   + geom_vline(xintercept=median(merged_expression_table_normalized_stage_III$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(merged_expression_table_normalized_stage_II$GC3), linetype="dashed", color =  "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))   + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 5), limits = c(25, 85))
    
    m4 <- ggplot(merged_expression_table_normalized_all_stages, aes(x=Stages, y=GC3))        +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()           +  geom_hline(yintercept=median(merged_expression_table_normalized_all_stages$T2), linetype="dashed", color = "red")    + ggtitle(paste("T2 All points: ", normalization_scheme,sep=""))                   +  stat_compare_means(comparisons = my_comparisons, method = "t.test") 
    m5 <- ggplot(merged_expression_table_normalized_all_stages, aes(x=Stages, y=Conections)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()   +  geom_hline(yintercept=median(merged_expression_table_normalized_all_stages$Conections), linetype="dashed", color = "red")    + ggtitle(paste("Connectivity All points: ", normalization_scheme,sep="")) +  stat_compare_means(comparisons = my_comparisons, method = "t.test") 
    m6 <- ggplot(merged_expression_table_normalized_all_stages, aes(x=Stages, y=Expr))       +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()         +      geom_hline(yintercept=median(merged_expression_table_normalized_all_stages$Expr), linetype="dashed", color = "red")    + ggtitle(paste("Expr. All points: ", normalization_scheme,sep=""))          +  stat_compare_means(comparisons = my_comparisons, method = "t.test") 

    m7<-ggplot(merged_expression_table_normalized_stage_I, aes(Conections, GC3, z = Expr))  + geom_point(aes(colour=Expr))   +  theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))+ geom_contour()  + xlim(0, 50)       + geom_vline(xintercept=mean(merged_expression_table_normalized_stage_I$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_I$GC3), linetype="dashed", color = "red")     + geom_vline(xintercept=median(merged_expression_table_normalized_stage_I$Conections), linetype="dashed", color = "yellow") +  geom_hline(yintercept=median(merged_expression_table_normalized_stage_I$T2), linetype="dashed", color =   "yellow")   + ylim(10,40) + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 5), limits = c(25, 85))
    m8<-ggplot(merged_expression_table_normalized_stage_II, aes(Conections, GC3, z = Expr))  + geom_point(aes(colour=Expr))  +  theme_bw() + ggtitle(paste(normalization_scheme,   ": All points Stage II Expr. ",sep=""))+ geom_contour()  + xlim(0, 50) + geom_vline(xintercept=mean(merged_expression_table_normalized_stage_II$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_II$GC3), linetype="dashed", color = "red")   + geom_vline(xintercept=median(merged_expression_table_normalized_stage_II$Conections), linetype="dashed", color = "yellow") +  geom_hline(yintercept=median(merged_expression_table_normalized_stage_II$T2), linetype="dashed", color =  "yellow")   + ylim(10,40) + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))       + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 5), limits = c(25, 85))
    m9<-ggplot(merged_expression_table_normalized_stage_III, aes(Conections, GC3, z = Expr))  + geom_point(aes(colour=Expr)) + theme_bw() + ggtitle(paste(normalization_scheme,  ": All points Stage III Expr. ",sep=""))+ geom_contour()  + xlim(0, 50) +  geom_vline(xintercept=mean(merged_expression_table_normalized_stage_III$Conections), linetype="dashed", color = "red")  +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_III$GC3), linetype="dashed", color = "red") + geom_vline(xintercept=median(merged_expression_table_normalized_stage_III$Conections), linetype="dashed", color = "yellow") +  geom_hline(yintercept=median(merged_expression_table_normalized_stage_III$T2), linetype="dashed", color = "yellow")   + ylim(10,40) +  guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))    + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 5), limits = c(25, 85))
  
    # FindClusters_resolution               
    png(filename=paste(output_dir,"countour_T2_Coonections_melt_",normalization_scheme,"_",TCGA_project,"_Stage_all_GC3.png",sep=""), width = 30, height = 50, res=600, units = "cm")  
            plot<-ggarrange(m7, m8, m9, m1,m2,m3,m4,m5,m6, nrow = 3,ncol = 3, common.legend = TRUE, legend="bottom")
            annotate_figure(plot, top = text_grob(TCGA_project, face = "bold", size = 14))
    dev.off()    
  
  # Histograms for T2
  h1<-ggplot(unique(merged_expression_table_normalized_all_sel[,c("ENSEMBL","T2","GC3","Stages")]), aes(x=T2, color=Stages)) +  geom_histogram(fill="white", alpha=0.5, position="identity") + theme_bw()     + xlim(0, 50)
  h2<-ggplot(unique(merged_expression_table_normalized_all_sel[,c("ENSEMBL","T2","GC3","Stages")]), aes(x=GC3, color=Stages)) +  geom_histogram(fill="white", alpha=0.5, position="identity") + theme_bw()    + xlim(0, 100)
  h3<-ggplot(unique(merged_expression_table_normalized_all_sel[,c("ENSEMBL","T2","Expr","Stages")]), aes(x=Expr, color=Stages)) +  geom_histogram(fill="white", alpha=0.5, position="identity") + theme_bw()  + xlim(0, 10000)

  # FindClusters_resolution               
  png(filename=paste(output_dir,"countour_T2_Coonections_melt_",normalization_scheme,"_",TCGA_project,"_histogram.png",sep=""), width = 25, height = 15, res=600, units = "cm")            
            plot<-ggarrange(h1, h2, h3, nrow = 1, common.legend = TRUE, legend="bottom") 
            annotate_figure(plot, top = text_grob(TCGA_project, face = "bold", size = 14))  
  dev.off()  
  
}


