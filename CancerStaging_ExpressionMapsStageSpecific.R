f#######################################################################################################################
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
Interactomes_GC3_T2_merged<-Interactomes_GC3_T2_merged[Interactomes_GC3_T2_merged$Conections <=100,]

# Filter by T2
Interactomes_GC3_T2_merged<-na.omit(Interactomes_GC3_T2_merged)
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
    # genes_stages_I
    genes_stages_I    <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_I",".tsv",sep=""), sep = '\t', header = TRUE)$gene #
    genes_stages_II   <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_II",".tsv",sep=""), sep = '\t', header = TRUE)$gene #
    genes_stages_III  <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_III",".tsv",sep=""), sep = '\t', header = TRUE)$gene #
  
    # Take also expression data from the normalization scheme set by "normalization_scheme"
    expression_table_normalized<-df_reads_count_all_projects[[normalization_scheme]]

    # Stage specific genes from each stage
    expression_table_normalized_stage_I  <-expression_table_normalized[genes_stages_I,]
    expression_table_normalized_stage_II <-expression_table_normalized[genes_stages_II,]
    expression_table_normalized_stage_III<-expression_table_normalized[genes_stages_III,]

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
    #######################################################################################################################################
    # Melt data.frame 
    merged_expression_table_normalized_all_stages <- rbind(merged_expression_table_normalized_stage_I,merged_expression_table_normalized_stage_II,merged_expression_table_normalized_stage_III)
  
    # Set the field stages
    Interactomes_GC3_T2_merged$Stages                   <-"overlapping"
    merged_expression_table_normalized_all_stages$Stages<-"overlapping"

    merged_expression_table_normalized_all_stages$ENSEMBL

    # Set stages 
    Interactomes_GC3_T2_merged[unique_stage_I,"Stages"]<-"Stage I"
    Interactomes_GC3_T2_merged[unique_stage_II,"Stages"]<-"Stage II"
    Interactomes_GC3_T2_merged[unique_stage_III,"Stages"]<-"Stage III"

    merged_expression_table_normalized_all_stages[merged_expression_table_normalized_all_stages$ENSEMBL %in% unique_stage_I,"Stages"]<-"Stage I"
    merged_expression_table_normalized_all_stages[merged_expression_table_normalized_all_stages$ENSEMBL %in% unique_stage_II,"Stages"]<-"Stage II"
    merged_expression_table_normalized_all_stages[merged_expression_table_normalized_all_stages$ENSEMBL %in% unique_stage_III,"Stages"]<-"Stage III"
    #########################################################################################################################################    
  
    ####################################################################################################################################################
    # I have two tables to be use.                                                                                                                     #
    # First , a table with with "T2", "GC3", "Conections", "ENSEMBL" "AveExp" per gene                                                                 #
    # Interactomes_GC3_T2_merged                                                                                                                       #
    # Second, a table with with "T2", "GC3", "Conections", "ENSEMBL" "Exp"    per patient                                                              #
    # melt_expression_interactomes                                                                                                                     #
    #############################################################################################################################################################
    Interactomes_GC3_T2_selected_Stage_I                       <-Interactomes_GC3_T2_merged_Stage_I[,c("T2","AveExp","Conections","GC3")]   #
    Interactomes_GC3_T2_selected_Stage_II                      <-Interactomes_GC3_T2_merged_Stage_II[,c("T2","AveExp","Conections","GC3")]  #
    Interactomes_GC3_T2_selected_Stage_III                     <-Interactomes_GC3_T2_merged_Stage_III[,c("T2","AveExp","Conections","GC3")] #
                                                                                                                                                                #
    # FindClusters_resolution          #                                                                                                                        #
    png(filename=paste(output_dir,"scatterplot3d_avg_",normalization_scheme,"_stage_I.png",sep=""), width = 24, height = 24, res=600, units = "cm")             #
            scatterplot3d(Interactomes_GC3_T2_selected_Stage_I[,c("T2","Conections","AveExp")], pch = 16,main="Stage I - Average Expression")       #
    dev.off()                                                                                                                                                   #
    png(filename=paste(output_dir,"scatterplot3d_avg_",normalization_scheme,"_stage_II.png",sep=""), width = 24, height = 24, res=600, units = "cm")            #
            scatterplot3d(Interactomes_GC3_T2_selected_Stage_II[,c("T2","Conections","AveExp")], pch = 16,main="Stage II - Average Expression")     #
    dev.off()                                                                                                                                                   #
    png(filename=paste(output_dir,"scatterplot3d_avg_",normalization_scheme,"_stage_III.png",sep=""), width = 24, height = 24, res=600, units = "cm")           #
            scatterplot3d(Interactomes_GC3_T2_selected_Stage_III[,c("T2","Conections","AveExp")], pch = 16,main="Stage III - Average Expression")   #
    dev.off()                                                                                                                                                   #
   ####################################################################################################################################################################
    # FindClusters_resolution          #                                                                                                                              #
    png(filename=paste(output_dir,"scatterplot3d_melt_",normalization_scheme,"_stage_I.png",sep=""), width = 24, height = 24, res=600, units = "cm")                  #
            scatterplot3d(merged_expression_table_normalized_stage_I[,c("T2","Conections",normalization_scheme)], pch = 16,main="Stage I- Expression per patient")    #
    dev.off()                                                                                                                                                         #
    # FindClusters_resolution          #                                                                                                                              #
    png(filename=paste(output_dir,"scatterplot3d_melt_",normalization_scheme,"_stage_II.png",sep=""), width = 24, height = 24, res=600, units = "cm")                 #
            scatterplot3d(merged_expression_table_normalized_stage_II[,c("T2","Conections",normalization_scheme)], pch = 16,main="Stage II- Expression per patient")  #
    dev.off()                                                                                                                                                         #
      # FindClusters_resolution          #                                                                                                                            #
    png(filename=paste(output_dir,"scatterplot3d_melt_",normalization_scheme,"_stage_III.png",sep=""), width = 24, height = 24, res=600, units = "cm")                #
            scatterplot3d(merged_expression_table_normalized_stage_III[,c("T2","Conections",normalization_scheme)], pch = 16,main="Stage III- Expression per patient")#
    dev.off()                                                                                                                                                         #
    ###################################################################################################################################################################
  
    # Only Variable Labels on the outside (no axis labels)
    Interactomes_GC3_T2_melt_Stage_all   <- ggpairs(merged_expression_table_normalized_stage_all[,c("T2","GC3",normalization_scheme,"Conections")], axisLabels = "none")
    
    # FindClusters_resolution
    png(filename=paste(output_dir,"correaltion_matrix_melt_",normalization_scheme,"_stages_all.png",sep=""), width = 20, height = 20, res=600, units = "cm")  
            Interactomes_GC3_T2_melt_Stage_all
    dev.off()
    ###################################################################################################################################################################
    colnames(merged_expression_table_normalized_stage_I)[6]<-"Expr"
    colnames(merged_expression_table_normalized_stage_II)[6]<-"Expr"
    colnames(merged_expression_table_normalized_stage_III)[6]<-"Expr"
    colnames(merged_expression_table_normalized_all_stages)[6]<-"Expr"

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
    merged_expression_table_normalized_stage_I$Exp_z_score        <-  calculate_z(merged_expression_table_normalized_stage_I$Expr,        mean(merged_expression_table_normalized_stage_I$Expr, na.rm = TRUE),sd(merged_expression_table_normalized_stage_I$Expr, na.rm = TRUE))  
    merged_expression_table_normalized_stage_I$Conections_z_score <-  calculate_z(merged_expression_table_normalized_stage_I$Conections,  mean(merged_expression_table_normalized_stage_I$Conections, na.rm = TRUE),sd(merged_expression_table_normalized_stage_I$Conections, na.rm = TRUE))    
    merged_expression_table_normalized_stage_I$T2_z_score         <-  calculate_z(merged_expression_table_normalized_stage_I$T2,          mean(merged_expression_table_normalized_stage_I$T2, na.rm = TRUE),sd(merged_expression_table_normalized_stage_I$T2, na.rm = TRUE))      

    merged_expression_table_normalized_stage_II$Exp_z_score        <-  calculate_z(merged_expression_table_normalized_stage_II$Expr,        mean(merged_expression_table_normalized_stage_II$Expr, na.rm = TRUE),sd(merged_expression_table_normalized_stage_II$Expr, na.rm = TRUE))  
    merged_expression_table_normalized_stage_II$Conections_z_score <-  calculate_z(merged_expression_table_normalized_stage_II$Conections,  mean(merged_expression_table_normalized_stage_II$Conections, na.rm = TRUE),sd(merged_expression_table_normalized_stage_II$Conections, na.rm = TRUE))    
    merged_expression_table_normalized_stage_II$T2_z_score         <-  calculate_z(merged_expression_table_normalized_stage_II$T2,          mean(merged_expression_table_normalized_stage_II$T2, na.rm = TRUE),sd(merged_expression_table_normalized_stage_II$T2, na.rm = TRUE))      

    merged_expression_table_normalized_stage_III$Exp_z_score        <-  calculate_z(merged_expression_table_normalized_stage_III$Expr,        mean(merged_expression_table_normalized_stage_III$Expr, na.rm = TRUE),sd(merged_expression_table_normalized_stage_III$Expr, na.rm = TRUE))  
    merged_expression_table_normalized_stage_III$Conections_z_score <-  calculate_z(merged_expression_table_normalized_stage_III$Conections,  mean(merged_expression_table_normalized_stage_III$Conections, na.rm = TRUE),sd(merged_expression_table_normalized_stage_III$Conections, na.rm = TRUE))    
    merged_expression_table_normalized_stage_III$T2_z_score         <-  calculate_z(merged_expression_table_normalized_stage_III$T2,          mean(merged_expression_table_normalized_stage_III$T2, na.rm = TRUE),sd(merged_expression_table_normalized_stage_III$T2, na.rm = TRUE))        

    #########################################################################################################################################
    # Visualize: Specify the comparisons you want
    my_comparisons <- list( c("Stage I", "Stage II"), c("Stage I", "Stage III"), c("Stage II", "Stage III"), c("Stage I", "overllapping"),c("Stage II", "overllapping"),c("Stage III", "overllapping"))
    #########################################################################################################################################
    # Three countour plots will be created
    # One with the average expression
    # Second with the expression per patient
    # Third with the z-score 
    #########################################################################################################################################
    # Filter up Average expression greater than zero
    merged_expression_table_normalized_stage_I  <-merged_expression_table_normalized_stage_I[merged_expression_table_normalized_stage_I$Expr>0 & merged_expression_table_normalized_stage_I$Expr<100,]
    merged_expression_table_normalized_stage_II <-merged_expression_table_normalized_stage_II[merged_expression_table_normalized_stage_II$Expr>0 & merged_expression_table_normalized_stage_II$Expr<100,]
    merged_expression_table_normalized_stage_III <-merged_expression_table_normalized_stage_III[merged_expression_table_normalized_stage_III$Expr>0  & merged_expression_table_normalized_stage_III$Expr<100,]
        
    m1<-ggplot(merged_expression_table_normalized_stage_I, aes(Conections, T2, z = Expr))  + geom_point(aes(colour=Expr)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))+ geom_contour()  + xlim(0, 50) + ylim(0, 60)        + geom_vline(xintercept=mean(merged_expression_table_normalized_stage_I$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_I$T2), linetype="dashed", color = "red")     + geom_vline(xintercept=median(merged_expression_table_normalized_stage_I$T2), linetype="dashed", color = "yellow")   +  geom_hline(yintercept=median(merged_expression_table_normalized_stage_I$Conections), linetype="dashed", color =   "yellow")   
    m2<-ggplot(merged_expression_table_normalized_stage_II, aes(Conections, T2, z = Expr))  + geom_point(aes(colour=Expr)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,   ": All points Stage II Expr. ",sep=""))+ geom_contour()  + xlim(0, 50) + ylim(0, 60)       + geom_vline(xintercept=mean(merged_expression_table_normalized_stage_II$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_II$T2), linetype="dashed", color = "red")   + geom_vline(xintercept=median(merged_expression_table_normalized_stage_II$T2), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(merged_expression_table_normalized_stage_II$Conections), linetype="dashed", color =  "yellow")
    m3<-ggplot(merged_expression_table_normalized_stage_III, aes(Conections, T2, z = Expr))  + geom_point(aes(colour=Expr)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,  ": All points Stage III Expr. ",sep=""))+ geom_contour()  + xlim(0, 50) + ylim(0, 60)      + geom_vline(xintercept=mean(merged_expression_table_normalized_stage_III$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_III$T2), linetype="dashed", color = "red") + geom_vline(xintercept=median(merged_expression_table_normalized_stage_III$T2), linetype="dashed", color = "yellow") +  geom_hline(yintercept=median(merged_expression_table_normalized_stage_III$Conections), linetype="dashed", color = "yellow")    
    
    m4 <- ggplot(merged_expression_table_normalized_all_stages, aes(x=Stages, y=T2)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()           +  geom_hline(yintercept=median(merged_expression_table_normalized_all_stages$T2), linetype="dashed", color = "red")    + ggtitle(paste("T2 All points: ", normalization_scheme,sep="")) +  stat_compare_means(comparisons = my_comparisons, method = "t.test")
    m5 <- ggplot(merged_expression_table_normalized_all_stages, aes(x=Stages, y=Conections)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()   +  geom_hline(yintercept=median(merged_expression_table_normalized_all_stages$Conections), linetype="dashed", color = "red")    + ggtitle(paste("Connectivity All points: ", normalization_scheme,sep="")) +  stat_compare_means(comparisons = my_comparisons, method = "t.test")
    m6 <- ggplot(merged_expression_table_normalized_all_stages, aes(x=Stages, y=Expr)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()   +  geom_hline(yintercept=median(merged_expression_table_normalized_all_stages$Expr), linetype="dashed", color = "red")    + ggtitle(paste("Expr. All points: ", normalization_scheme,sep="")) +  stat_compare_means(comparisons = my_comparisons, method = "t.test")

    # FindClusters_resolution               
    png(filename=paste(output_dir,"countour_T2_Coonections_melt_",normalization_scheme,"_Stage_all.png",sep=""), width = 25, height = 20, res=600, units = "cm")  
            ggarrange(m1,m2,m3,m4,m5,m6, nrow = 2,ncol = 3, common.legend = TRUE, legend="bottom")
    dev.off()  
      
    # Filter up Average expression greater than zero        
    m1<-ggplot(merged_expression_table_normalized_stage_I, aes(Conections_z_score, T2_z_score, z = Exp_z_score))  + geom_point(aes(colour=Exp_z_score)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,    ": Stage I Z-score",sep=""))+ geom_contour()        + geom_vline(xintercept=mean(merged_expression_table_normalized_stage_I$Conections_z_score), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_I$T2_z_score), linetype="dashed", color = "red")        + geom_vline(xintercept=median(merged_expression_table_normalized_stage_I$T2_z_score), linetype="dashed", color = "yellow") +  geom_hline(yintercept=median(merged_expression_table_normalized_stage_I$Conections_z_score), linetype="dashed", color = "yellow") 
    m2<-ggplot(merged_expression_table_normalized_stage_II, aes(Conections_z_score, T2_z_score, z = Exp_z_score))  + geom_point(aes(colour=Exp_z_score)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,    ": Stage II Z-score",sep=""))+ geom_contour()      + geom_vline(xintercept=mean(merged_expression_table_normalized_stage_II$Conections_z_score), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_II$T2_z_score), linetype="dashed", color = "red")      + geom_vline(xintercept=median(merged_expression_table_normalized_stage_II$T2_z_score), linetype="dashed", color = "yellow") +  geom_hline(yintercept=median(merged_expression_table_normalized_stage_II$Conections_z_score), linetype="dashed", color = "yellow")      
    m3<-ggplot(merged_expression_table_normalized_stage_III, aes(Conections_z_score, T2_z_score, z = Exp_z_score))  + geom_point(aes(colour=Exp_z_score)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,    ": Stage III Z-score",sep=""))+ geom_contour()    + geom_vline(xintercept=mean(merged_expression_table_normalized_stage_III$Conections_z_score), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_III$T2_z_score), linetype="dashed", color = "red")     + geom_vline(xintercept=median(merged_expression_table_normalized_stage_III$T2_z_score), linetype="dashed", color = "yellow") +  geom_hline(yintercept=median(merged_expression_table_normalized_stage_III$Conections_z_score), linetype="dashed", color = "yellow")  

    m4 <- ggplot(merged_expression_table_normalized_all, aes(x=Stage, y=T2_z_score)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()           +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged$T2_z_score), linetype="dashed", color = "red")    + ggtitle(paste("T2 z-score: ", normalization_scheme,sep="")) +  stat_compare_means(comparisons = my_comparisons, method = "t.test")
    m5 <- ggplot(merged_expression_table_normalized_all, aes(x=Stage, y=Conections_z_score)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()   +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged$Conections_z_score), linetype="dashed", color = "red")    + ggtitle(paste("Connectivity z-score: ", normalization_scheme,sep="")) +  stat_compare_means(comparisons = my_comparisons, method = "t.test")
    m6 <- ggplot(merged_expression_table_normalized_all, aes(x=Stage, y=Exp_z_score)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()   +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged$Exp_z_score), linetype="dashed", color = "red")    + ggtitle(paste("Expr. Z-score: ", normalization_scheme,sep="")) +  stat_compare_means(comparisons = my_comparisons, method = "t.test")
      
    # FindClusters_resolution               
    png(filename=paste(output_dir,"countour_T2_Coonections_zscore_",normalization_scheme,"_Stage_all.png",sep=""), width = 25, height = 20, res=600, units = "cm")  
            ggarrange(m1,m2,m3,m4,m5,m6, nrow = 2,ncol = 3, common.legend = TRUE, legend="bottom")
    dev.off()  
    #########################################################################################################################################    
}

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
    # genes_stages_I                                                
    unique_genes_stages_I    <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","unique_stage_I",".tsv",sep=""), sep = '\t', header = TRUE)$gene
    unique_genes_stages_II   <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","unique_stage_II",".tsv",sep=""), sep = '\t', header = TRUE)$gene
    unique_genes_stages_III  <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","unique_stage_III",".tsv",sep=""), sep = '\t', header = TRUE)$gene

    # Stage specific genes from each stage
    expression_table_normalized_stage_I  <-expression_table_normalized[unique_genes_stages_I,]
    expression_table_normalized_stage_II <-expression_table_normalized[unique_genes_stages_II,]
    expression_table_normalized_stage_III<-expression_table_normalized[unique_genes_stages_III,]

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
    ####################################################################################################################################################
    # I have two tables to be use.                                                                                                                     #
    # First , a table with with "T2", "GC3", "Conections", "ENSEMBL" "AveExp" per gene                                                                 #
    # Interactomes_GC3_T2_merged                                                                                                                       #
    # Second, a table with with "T2", "GC3", "Conections", "ENSEMBL" "Exp"    per patient                                                              #
    # melt_expression_interactomes                                                                                                                     #
    #############################################################################################################################################################
    Interactomes_GC3_T2_selected_Stage_I                       <-Interactomes_GC3_T2_merged_Stage_I[,c("T2","AveExp","Conections","GC3")]   #
    Interactomes_GC3_T2_selected_Stage_II                      <-Interactomes_GC3_T2_merged_Stage_II[,c("T2","AveExp","Conections","GC3")]  #
    Interactomes_GC3_T2_selected_Stage_III                     <-Interactomes_GC3_T2_merged_Stage_III[,c("T2","AveExp","Conections","GC3")] #
                                                                                                                                                                #
    # FindClusters_resolution          #                                                                                                                        #
    png(filename=paste(output_dir,"scatterplot3d_unique_avg_",normalization_scheme,"_stage_I.png",sep=""), width = 24, height = 24, res=600, units = "cm")             #
            scatterplot3d(Interactomes_GC3_T2_selected_Stage_I[,c("T2","Conections","AveExp")], pch = 16,main="Stage I - Average Expression")       #
    dev.off()                                                                                                                                                   #
    png(filename=paste(output_dir,"scatterplot3d_unique_avg_",normalization_scheme,"_stage_II.png",sep=""), width = 24, height = 24, res=600, units = "cm")            #
            scatterplot3d(Interactomes_GC3_T2_selected_Stage_II[,c("T2","Conections","AveExp")], pch = 16,main="Stage II - Average Expression")     #
    dev.off()                                                                                                                                                   #
    png(filename=paste(output_dir,"scatterplot3d_unique_avg_",normalization_scheme,"_stage_III.png",sep=""), width = 24, height = 24, res=600, units = "cm")           #
            scatterplot3d(Interactomes_GC3_T2_selected_Stage_III[,c("T2","Conections","AveExp")], pch = 16,main="Stage III - Average Expression")   #
    dev.off()                                                                                                                                                   #
   ####################################################################################################################################################################
    # FindClusters_resolution          #                                                                                                                              #
    png(filename=paste(output_dir,"scatterplot3d_unique_melt_",normalization_scheme,"_stage_I.png",sep=""), width = 24, height = 24, res=600, units = "cm")                  #
            scatterplot3d(merged_expression_table_normalized_stage_I[,c("T2","Conections",normalization_scheme)], pch = 16,main="Stage I- Expression per patient")    #
    dev.off()                                                                                                                                                         #
    # FindClusters_resolution          #                                                                                                                              #
    png(filename=paste(output_dir,"scatterplot3d_unique_melt_",normalization_scheme,"_stage_II.png",sep=""), width = 24, height = 24, res=600, units = "cm")                 #
            scatterplot3d(merged_expression_table_normalized_stage_II[,c("T2","Conections",normalization_scheme)], pch = 16,main="Stage II- Expression per patient")  #
    dev.off()                                                                                                                                                         #
      # FindClusters_resolution          #                                                                                                                            #
    png(filename=paste(output_dir,"scatterplot3d_unique_melt_",normalization_scheme,"_stage_III.png",sep=""), width = 24, height = 24, res=600, units = "cm")                #
            scatterplot3d(merged_expression_table_normalized_stage_III[,c("T2","Conections",normalization_scheme)], pch = 16,main="Stage III- Expression per patient")#
    dev.off()                                                                                                                                                         #
    ###################################################################################################################################################################
    # Melt data.frame 
    merged_expression_table_normalized_stage_all<-rbind(merged_expression_table_normalized_stage_I,merged_expression_table_normalized_stage_II,merged_expression_table_normalized_stage_III)
  
    # Only Variable Labels on the outside (no axis labels)
    Interactomes_GC3_T2_melt_Stage_all   <- ggpairs(merged_expression_table_normalized_stage_all[,c("T2","GC3",normalization_scheme,"Conections")], axisLabels = "none")
    
    # FindClusters_resolution
    png(filename=paste(output_dir,"correaltion_matrix_unique_melt_",normalization_scheme,"_stages_all.png",sep=""), width = 20, height = 20, res=600, units = "cm")  
            Interactomes_GC3_T2_melt_Stage_all
    dev.off()
    ###################################################################################################################################################################
    colnames(merged_expression_table_normalized_stage_I)[6]<-"Expr"
    colnames(merged_expression_table_normalized_stage_II)[6]<-"Expr"
    colnames(merged_expression_table_normalized_stage_III)[6]<-"Expr"

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
    merged_expression_table_normalized_stage_I$Exp_z_score        <-  calculate_z(merged_expression_table_normalized_stage_I$Expr,        mean(merged_expression_table_normalized_stage_I$Expr, na.rm = TRUE),sd(merged_expression_table_normalized_stage_I$Expr, na.rm = TRUE))  
    merged_expression_table_normalized_stage_I$Conections_z_score <-  calculate_z(merged_expression_table_normalized_stage_I$Conections,  mean(merged_expression_table_normalized_stage_I$Conections, na.rm = TRUE),sd(merged_expression_table_normalized_stage_I$Conections, na.rm = TRUE))    
    merged_expression_table_normalized_stage_I$T2_z_score         <-  calculate_z(merged_expression_table_normalized_stage_I$T2,          mean(merged_expression_table_normalized_stage_I$T2, na.rm = TRUE),sd(merged_expression_table_normalized_stage_I$T2, na.rm = TRUE))      

    merged_expression_table_normalized_stage_II$Exp_z_score        <-  calculate_z(merged_expression_table_normalized_stage_II$Expr,        mean(merged_expression_table_normalized_stage_II$Expr, na.rm = TRUE),sd(merged_expression_table_normalized_stage_II$Expr, na.rm = TRUE))  
    merged_expression_table_normalized_stage_II$Conections_z_score <-  calculate_z(merged_expression_table_normalized_stage_II$Conections,  mean(merged_expression_table_normalized_stage_II$Conections, na.rm = TRUE),sd(merged_expression_table_normalized_stage_II$Conections, na.rm = TRUE))    
    merged_expression_table_normalized_stage_II$T2_z_score         <-  calculate_z(merged_expression_table_normalized_stage_II$T2,          mean(merged_expression_table_normalized_stage_II$T2, na.rm = TRUE),sd(merged_expression_table_normalized_stage_II$T2, na.rm = TRUE))      

    merged_expression_table_normalized_stage_III$Exp_z_score        <-  calculate_z(merged_expression_table_normalized_stage_III$Expr,        mean(merged_expression_table_normalized_stage_III$Expr, na.rm = TRUE),sd(merged_expression_table_normalized_stage_III$Expr, na.rm = TRUE))  
    merged_expression_table_normalized_stage_III$Conections_z_score <-  calculate_z(merged_expression_table_normalized_stage_III$Conections,  mean(merged_expression_table_normalized_stage_III$Conections, na.rm = TRUE),sd(merged_expression_table_normalized_stage_III$Conections, na.rm = TRUE))    
    merged_expression_table_normalized_stage_III$T2_z_score         <-  calculate_z(merged_expression_table_normalized_stage_III$T2,          mean(merged_expression_table_normalized_stage_III$T2, na.rm = TRUE),sd(merged_expression_table_normalized_stage_III$T2, na.rm = TRUE))        
            
    #########################################################################################################################################
    # Three countour plots will be created
    # One with the average expression
    # Second with the expression per patient
    # Third with the z-score 
    #########################################################################################################################################
    # Filter up Average expression greater than zero
    merged_expression_table_normalized_stage_I  <-merged_expression_table_normalized_stage_I[merged_expression_table_normalized_stage_I$Expr>0,]
    merged_expression_table_normalized_stage_II <-merged_expression_table_normalized_stage_II[merged_expression_table_normalized_stage_II$Expr>0,]
    merged_expression_table_normalized_stage_III <-merged_expression_table_normalized_stage_III[merged_expression_table_normalized_stage_III$Expr>0,]
        
    m1<-ggplot(merged_expression_table_normalized_stage_I, aes(Conections, T2, z = Expr))  + geom_point(aes(colour=Expr)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,    ": Stage I Expr. ",sep=""))+ geom_contour()  + xlim(0, 50) + ylim(0, 50)      #  + theme(legend.position="none")        a
    m2<-ggplot(merged_expression_table_normalized_stage_II, aes(Conections, T2, z = Expr))  + geom_point(aes(colour=Expr)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,   ": Stage II Expr. ",sep=""))+ geom_contour()  + xlim(0, 50) + ylim(0, 50)      #  + theme(legend.position="none")        a
    m3<-ggplot(merged_expression_table_normalized_stage_III, aes(Conections, T2, z = Expr))  + geom_point(aes(colour=Expr)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,  ": Stage III Expr. ",sep=""))+ geom_contour()  + xlim(0, 50) + ylim(0, 50)      #  + theme(legend.position="none")        a
    
    # FindClusters_resolution               
    png(filename=paste(output_dir,"countour_T2_Coonections_unique_melt_",normalization_scheme,"_Stage_all.png",sep=""), width = 25, height = 10, res=600, units = "cm")  
            ggarrange(m1,m2,m3,nrow = 1,ncol = 3, common.legend = TRUE, legend="bottom")
    dev.off()  
    #########################################################################################################################################    
    # Filter up Average expression greater than zero
    Interactomes_GC3_T2_selected_Stage_I  <-Interactomes_GC3_T2_selected_Stage_I[Interactomes_GC3_T2_selected_Stage_I$AveExp>0,]
    Interactomes_GC3_T2_selected_Stage_II <-Interactomes_GC3_T2_selected_Stage_II[Interactomes_GC3_T2_selected_Stage_II$AveExp>0,]
    Interactomes_GC3_T2_selected_Stage_III <-Interactomes_GC3_T2_selected_Stage_III[Interactomes_GC3_T2_selected_Stage_III$AveExp>0,]
        
    m1<-ggplot(Interactomes_GC3_T2_selected_Stage_I,   aes(Conections, T2, z = AveExp))  + geom_point(aes(colour=AveExp)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,  ": Stage I Average Expr.",sep=""))+ geom_contour()  + xlim(0, 50) + ylim(0, 50)      #  + theme(legend.position="none")        a
    m2<-ggplot(Interactomes_GC3_T2_selected_Stage_II,  aes(Conections, T2, z = AveExp))  + geom_point(aes(colour=AveExp)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,  ": Stage II Average Expr.",sep=""))+ geom_contour()  + xlim(0, 50) + ylim(0, 50)      #  + theme(legend.position="none")        a
    m3<-ggplot(Interactomes_GC3_T2_selected_Stage_III, aes(Conections, T2, z = AveExp))  + geom_point(aes(colour=AveExp)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,  ": Stage III Average Expr.",sep=""))+ geom_contour()  + xlim(0, 50) + ylim(0, 50)      #  + theme(legend.position="none")        a
    
    # FindClusters_resolution               
    png(filename=paste(output_dir,"countour_T2_Coonections_unique_avg_",normalization_scheme,"_Stage_all.png",sep=""), width = 25, height = 10, res=600, units = "cm")  
            ggarrange(m1,m2,m3,nrow = 1,ncol = 3, common.legend = TRUE, legend="bottom")
    dev.off()  
    #########################################################################################################################################

    #########################################################################################################################################
    # Filter up Average expression greater than zero
    merged_expression_table_normalized_stage_I  <-merged_expression_table_normalized_stage_I[merged_expression_table_normalized_stage_I$Exp_z_score>0,]
    merged_expression_table_normalized_stage_II <-merged_expression_table_normalized_stage_II[merged_expression_table_normalized_stage_II$Exp_z_score>0,]
    merged_expression_table_normalized_stage_III <-merged_expression_table_normalized_stage_III[merged_expression_table_normalized_stage_III$Exp_z_score>0,]
        
    m1<-ggplot(merged_expression_table_normalized_stage_I, aes(Conections_z_score, T2_z_score, z = Exp_z_score))  + geom_point(aes(colour=Exp_z_score)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,    ": Stage I Z-score",sep=""))+ geom_contour()     + theme(legend.position="none")        
    m2<-ggplot(merged_expression_table_normalized_stage_II, aes(Conections_z_score, T2_z_score, z = Exp_z_score))  + geom_point(aes(colour=Exp_z_score)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,    ": Stage II Z-score",sep=""))+ geom_contour()    + theme(legend.position="none")        
    m3<-ggplot(merged_expression_table_normalized_stage_III, aes(Conections_z_score, T2_z_score, z = Exp_z_score))  + geom_point(aes(colour=Exp_z_score)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,    ": Stage III Z-score",sep=""))+ geom_contour()    + theme(legend.position="none")        
    
    # FindClusters_resolution               
    png(filename=paste(output_dir,"countour_T2_Coonections_unique_zscore_",normalization_scheme,"_Stage_all.png",sep=""), width = 25, height = 10, res=600, units = "cm")  
            ggarrange(m1,m2,m3,nrow = 1,ncol = 3, common.legend = TRUE, legend="bottom")
    dev.off()  
    #########################################################################################################################################    
  
}








merged_expression_table_normalized_stage_I$Stage<-"Stage I"
merged_expression_table_normalized_stage_II$Stage<-"Stage II"	
merged_expression_table_normalized_stage_III$Stage<-"Stage III"

Interactomes_GC3_T2_selected_Stage_I$Stage<-"Stage I"
Interactomes_GC3_T2_selected_Stage_II$Stage<-"Stage II"	
Interactomes_GC3_T2_selected_Stage_III$Stage<-"Stage III"					

# Merge the three stages
merged_expression_table_normalized_all<-rbind(merged_expression_table_normalized_stage_I,merged_expression_table_normalized_stage_II,merged_expression_table_normalized_stage_III)   

# Merge the three stages
Interactomes_GC3_T2_selected_all<-rbind(Interactomes_GC3_T2_selected_Stage_I,Interactomes_GC3_T2_selected_Stage_II,Interactomes_GC3_T2_selected_Stage_III)
#########################################################################################################################################
# Three countour plots will be created
# One with the average expression
# Second with the expression per patient
# Third with the z-score 
#########################################################################################################################################
# Filter up Average expression greater than zero
merged_expression_table_normalized_all  <-merged_expression_table_normalized_all[merged_expression_table_normalized_all$Expr>0,]
Interactomes_GC3_T2_selected_all  <-Interactomes_GC3_T2_selected_all[Interactomes_GC3_T2_selected_all$AveExp>0,]
  

m1<-ggplot(merged_expression_table_normalized_all, aes(Conections, T2, z = Expr))  + geom_point(aes(colour=Expr, shape=Stage)) + geom_density_2d() + theme_bw() 
m2<-ggplot(Interactomes_GC3_T2_selected_all, aes(Conections, T2, z = AveExp))  + geom_point(aes(colour=AveExp, shape=Stage)) + geom_density_2d() + theme_bw() 
m3<-ggplot(merged_expression_table_normalized_all, aes(Conections_z_score, T2_z_score, z = T2_z_score))  + geom_point(aes(colour=T2_z_score,shape=Stage)) + geom_density_2d() + theme_bw() 	

# FindClusters_resolution               
png(filename=paste(output_dir,"countour_T2_Coonections_all_combined_",normalization_scheme,"_Stage_all.png",sep=""), width = 25, height = 10, res=600, units = "cm")  
    ggarrange(m1,m2,m3,nrow = 1,ncol = 3, common.legend = TRUE, legend="bottom")
dev.off()  




















####################################################################################################################################################
# Z-score normalization
# 15-10-2024
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
    # genes_stages_I
    genes_stages_I    <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_I",".tsv",sep=""), sep = '\t', header = TRUE)$gene #
    genes_stages_II   <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_II",".tsv",sep=""), sep = '\t', header = TRUE)$gene #
    genes_stages_III  <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_III",".tsv",sep=""), sep = '\t', header = TRUE)$gene #
  
    # Take also expression data from the normalization scheme set by "normalization_scheme"
    expression_table_normalized<-df_reads_count_all_projects[[normalization_scheme]]

    # Stage specific genes from each stage
    expression_table_normalized_stage_I  <-expression_table_normalized[genes_stages_I,]
    expression_table_normalized_stage_II <-expression_table_normalized[genes_stages_II,]
    expression_table_normalized_stage_III<-expression_table_normalized[genes_stages_III,]

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
    colnames(merged_expression_table_normalized_stage_I)[6]<-"Expr"
    colnames(merged_expression_table_normalized_stage_II)[6]<-"Expr"
    colnames(merged_expression_table_normalized_stage_III)[6]<-"Expr"

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
    merged_expression_table_normalized_stage_I$Exp_z_score        <-  (merged_expression_table_normalized_stage_I$Expr-mean(merged_expression_table_normalized_stage_I$Expr))/sd(merged_expression_table_normalized_stage_I$Expr)
    merged_expression_table_normalized_stage_I$Conections_z_score <-  (merged_expression_table_normalized_stage_I$Conections-mean(merged_expression_table_normalized_stage_I$Conections))/sd(merged_expression_table_normalized_stage_I$Conections) 
    merged_expression_table_normalized_stage_I$T2_z_score         <-  (merged_expression_table_normalized_stage_I$T2-mean(merged_expression_table_normalized_stage_I$T2))/sd(merged_expression_table_normalized_stage_I$T2) 

    merged_expression_table_normalized_stage_II$Exp_z_score        <-  (merged_expression_table_normalized_stage_II$Expr-mean(merged_expression_table_normalized_stage_II$Expr))/sd(merged_expression_table_normalized_stage_II$Expr)
    merged_expression_table_normalized_stage_II$Conections_z_score <-  (merged_expression_table_normalized_stage_II$Conections-mean(merged_expression_table_normalized_stage_II$Conections))/sd(merged_expression_table_normalized_stage_II$Conections) 
    merged_expression_table_normalized_stage_II$T2_z_score         <-  (merged_expression_table_normalized_stage_II$T2-mean(merged_expression_table_normalized_stage_II$T2))/sd(merged_expression_table_normalized_stage_II$T2) 
  
    merged_expression_table_normalized_stage_III$Exp_z_score        <-  (merged_expression_table_normalized_stage_III$Expr-mean(merged_expression_table_normalized_stage_III$Expr))/sd(merged_expression_table_normalized_stage_III$Expr)
    merged_expression_table_normalized_stage_III$Conections_z_score <-  (merged_expression_table_normalized_stage_III$Conections-mean(merged_expression_table_normalized_stage_III$Conections))/sd(merged_expression_table_normalized_stage_III$Conections) 
    merged_expression_table_normalized_stage_III$T2_z_score         <-  (merged_expression_table_normalized_stage_III$T2-mean(merged_expression_table_normalized_stage_III$T2))/sd(merged_expression_table_normalized_stage_III$T2)   
    #########################################################################################################################################
    # Three countour plots will be created
    # One with the average expression
    # Second with the expression per patient
    # Third with the z-score 
    #########################################################################################################################################
    # Filter up Average expression greater than zero
    merged_expression_table_normalized_stage_I  <-merged_expression_table_normalized_stage_I[merged_expression_table_normalized_stage_I$Expr>0,]
    merged_expression_table_normalized_stage_II <-merged_expression_table_normalized_stage_II[merged_expression_table_normalized_stage_II$Expr>0,]
    merged_expression_table_normalized_stage_III <-merged_expression_table_normalized_stage_III[merged_expression_table_normalized_stage_III$Expr>0,]
        
    m1<-ggplot(merged_expression_table_normalized_stage_I, aes(Conections, T2, z = Expr))  + geom_point(aes(colour=Expr)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,    ": Stage I Expr. ",sep=""))+ geom_contour()  + xlim(0, 50) + ylim(0, 50)      #  + theme(legend.position="none")        a
    m2<-ggplot(merged_expression_table_normalized_stage_II, aes(Conections, T2, z = Expr))  + geom_point(aes(colour=Expr)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,   ": Stage II Expr. ",sep=""))+ geom_contour()  + xlim(0, 50) + ylim(0, 50)      #  + theme(legend.position="none")        a
    m3<-ggplot(merged_expression_table_normalized_stage_III, aes(Conections, T2, z = Expr))  + geom_point(aes(colour=Expr)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,  ": Stage III Expr. ",sep=""))+ geom_contour()  + xlim(0, 50) + ylim(0, 50)      #  + theme(legend.position="none")        a
    
    # FindClusters_resolution               
    png(filename=paste(output_dir,"countour_T2_Coonections_unique_melt_",normalization_scheme,"_Stage_all.png",sep=""), width = 25, height = 10, res=600, units = "cm")  
            ggarrange(m1,m2,m3,nrow = 1,ncol = 3, common.legend = TRUE, legend="bottom")
    dev.off()  
    #########################################################################################################################################    
    # Filter up Average expression greater than zero
    Interactomes_GC3_T2_selected_Stage_I  <-Interactomes_GC3_T2_selected_Stage_I[Interactomes_GC3_T2_selected_Stage_I$AveExp>0,]
    Interactomes_GC3_T2_selected_Stage_II <-Interactomes_GC3_T2_selected_Stage_II[Interactomes_GC3_T2_selected_Stage_II$AveExp>0,]
    Interactomes_GC3_T2_selected_Stage_III <-Interactomes_GC3_T2_selected_Stage_III[Interactomes_GC3_T2_selected_Stage_III$AveExp>0,]
        
    m1<-ggplot(Interactomes_GC3_T2_selected_Stage_I,   aes(Conections, T2, z = AveExp))  + geom_point(aes(colour=AveExp)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,  ": Stage I Average Expr.",sep=""))+ geom_contour()  + xlim(0, 50) + ylim(0, 50)      #  + theme(legend.position="none")        a
    m2<-ggplot(Interactomes_GC3_T2_selected_Stage_II,  aes(Conections, T2, z = AveExp))  + geom_point(aes(colour=AveExp)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,  ": Stage II Average Expr.",sep=""))+ geom_contour()  + xlim(0, 50) + ylim(0, 50)      #  + theme(legend.position="none")        a
    m3<-ggplot(Interactomes_GC3_T2_selected_Stage_III, aes(Conections, T2, z = AveExp))  + geom_point(aes(colour=AveExp)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,  ": Stage III Average Expr.",sep=""))+ geom_contour()  + xlim(0, 50) + ylim(0, 50)      #  + theme(legend.position="none")        a
    
    # FindClusters_resolution               
    png(filename=paste(output_dir,"countour_T2_Coonections_unique_avg_",normalization_scheme,"_Stage_all.png",sep=""), width = 25, height = 10, res=600, units = "cm")  
            ggarrange(m1,m2,m3,nrow = 1,ncol = 3, common.legend = TRUE, legend="bottom")
    dev.off()  
    #########################################################################################################################################            
    m1<-ggplot(merged_expression_table_normalized_stage_I, aes(Conections_z_score, T2_z_score, z = Exp_z_score))  + geom_point(aes(colour=Exp_z_score)) + theme_bw() + ggtitle(paste(normalization_scheme,    ": Stage I Z-score (Value - mean(Value))/sd(Value)",sep=""))+ geom_contour()     + theme(legend.position="none")         + geom_vline(xintercept=mean(merged_expression_table_normalized_stage_I$T2_z_score), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_I$Conections_z_score), linetype="dashed", color = "red")   +   xlim(0, 10) + ylim(0, 30)
    m2<-ggplot(merged_expression_table_normalized_stage_II, aes(Conections_z_score, T2_z_score, z = Exp_z_score))  + geom_point(aes(colour=Exp_z_score)) + theme_bw() + ggtitle(paste(normalization_scheme,    ": Stage II Z-score (Value - mean(Value))/sd(Value)",sep=""))+ geom_contour()    + theme(legend.position="none")        + geom_vline(xintercept=mean(merged_expression_table_normalized_stage_II$T2_z_score), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_II$Conections_z_score), linetype="dashed", color = "red")   + xlim(0, 10) + ylim(0, 30)
    m3<-ggplot(merged_expression_table_normalized_stage_III, aes(Conections_z_score, T2_z_score, z = Exp_z_score))  + geom_point(aes(colour=Exp_z_score)) + theme_bw() + ggtitle(paste(normalization_scheme,    ": Stage III Z-score (Value - mean(Value))/sd(Value)",sep=""))+ geom_contour()    + theme(legend.position="none")      + geom_vline(xintercept=mean(merged_expression_table_normalized_stage_III$T2_z_score), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_III$Conections_z_score), linetype="dashed", color = "red") +  xlim(0, 10) + ylim(0, 30)     
      
    # FindClusters_resolution               
    png(filename=paste(output_dir,"countour_T2_Coonections_unique_zscore_",normalization_scheme,"_Stage_all.png",sep=""), width = 15, height = 25, res=600, units = "cm")  
            ggarrange(m1,m2,m3,nrow = 3,ncol = 1, common.legend = TRUE, legend="bottom")
    dev.off()  
    #########################################################################################################################################
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
    merged_expression_table_normalized_stage_I$Exp_z_score2        <-  scale(merged_expression_table_normalized_stage_I$Expr, center = TRUE, scale = FALSE)
    merged_expression_table_normalized_stage_I$Conections_z_score2 <-  scale(merged_expression_table_normalized_stage_I$Conections, center = TRUE, scale = FALSE)
    merged_expression_table_normalized_stage_I$T2_z_score2         <-  scale(merged_expression_table_normalized_stage_I$T2, center = TRUE, scale = FALSE)

    merged_expression_table_normalized_stage_II$Exp_z_score2         <-  scale(merged_expression_table_normalized_stage_II$Expr, center = TRUE, scale = FALSE)
    merged_expression_table_normalized_stage_II$Conections_z_score2  <-  scale(merged_expression_table_normalized_stage_II$Conections,  mean(merged_expression_table_normalized_stage_II$Conections, na.rm = TRUE),sd(merged_expression_table_normalized_stage_II$Conections, na.rm = TRUE))    
    merged_expression_table_normalized_stage_II$T2_z_score2          <-  scale(merged_expression_table_normalized_stage_II$T2,          mean(merged_expression_table_normalized_stage_II$T2, na.rm = TRUE),sd(merged_expression_table_normalized_stage_II$T2, na.rm = TRUE))      

    merged_expression_table_normalized_stage_III$Exp_z_score2        <-  scale(merged_expression_table_normalized_stage_III$Expr, center = TRUE, scale = FALSE)
    merged_expression_table_normalized_stage_III$Conections_z_score2 <-  scale(merged_expression_table_normalized_stage_III$Conections,  mean(merged_expression_table_normalized_stage_III$Conections, na.rm = TRUE),sd(merged_expression_table_normalized_stage_III$Conections, na.rm = TRUE))    
    merged_expression_table_normalized_stage_III$T2_z_score2         <-  scale(merged_expression_table_normalized_stage_III$T2,          mean(merged_expression_table_normalized_stage_III$T2, na.rm = TRUE),sd(merged_expression_table_normalized_stage_III$T2, na.rm = TRUE))          
    
    m1<-ggplot(merged_expression_table_normalized_stage_I, aes(Conections_z_score2, T2_z_score, z = Exp_z_score2))  + geom_point(aes(colour=Exp_z_score2)) + theme_bw() + ggtitle(paste(normalization_scheme,    ": Stage I Scale",sep=""))+ geom_contour()     + theme(legend.position="none")         + geom_vline(xintercept=mean(merged_expression_table_normalized_stage_I$T2_z_score), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_I$Conections_z_score), linetype="dashed", color = "red")   +   xlim(-1, 1) + ylim(-1, 1)
    m2<-ggplot(merged_expression_table_normalized_stage_II, aes(Conections_z_score2, T2_z_score, z = Exp_z_score2))  + geom_point(aes(colour=Exp_z_score2)) + theme_bw() + ggtitle(paste(normalization_scheme,    ": Stage II Scale",sep=""))+ geom_contour()    + theme(legend.position="none")        + geom_vline(xintercept=mean(merged_expression_table_normalized_stage_II$T2_z_score), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_II$Conections_z_score), linetype="dashed", color = "red")   + xlim(-1, 1) + ylim(-1, 1)
    m3<-ggplot(merged_expression_table_normalized_stage_III, aes(Conections_z_score2, T2_z_score, z = Exp_z_score2))  + geom_point(aes(colour=Exp_z_score2)) + theme_bw() + ggtitle(paste(normalization_scheme,    ": Stage III Scale",sep=""))+ geom_contour()    + theme(legend.position="none")      + geom_vline(xintercept=mean(merged_expression_table_normalized_stage_III$T2_z_score), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_III$Conections_z_score), linetype="dashed", color = "red") +  xlim(-1, 1) + ylim(-1, 1)     

      
    # FindClusters_resolution               
    png(filename=paste(output_dir,"countour_T2_Coonections_unique_zscore_",normalization_scheme,"_Stage_all.png",sep=""), width = 25, height = 10, res=600, units = "cm")  
            ggarrange(m1,m2,m3,nrow = 1,ncol = 3, common.legend = TRUE, legend="bottom")
    dev.off()

}










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
    # genes_stages_I                                                
    unique_genes_stages_I    <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","unique_stage_I",".tsv",sep=""), sep = '\t', header = TRUE)$gene
    unique_genes_stages_II   <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","unique_stage_II",".tsv",sep=""), sep = '\t', header = TRUE)$gene
    unique_genes_stages_III  <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","unique_stage_III",".tsv",sep=""), sep = '\t', header = TRUE)$gene

    # Stage specific genes from each stage
    expression_table_normalized_stage_I  <-expression_table_normalized[unique_genes_stages_I,]
    expression_table_normalized_stage_II <-expression_table_normalized[unique_genes_stages_II,]
    expression_table_normalized_stage_III<-expression_table_normalized[unique_genes_stages_III,]

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
    merged_expression_table_normalized_stage_I$Stage<-"Stage I"
    merged_expression_table_normalized_stage_II$Stage<-"Stage II"	
    merged_expression_table_normalized_stage_III$Stage<-"Stage III"
    
    Interactomes_GC3_T2_merged_Stage_I$Stage<-"Stage I"
    Interactomes_GC3_T2_merged_Stage_II$Stage<-"Stage II"	
    Interactomes_GC3_T2_merged_Stage_III$Stage<-"Stage III"					
    
    # Merge the three stages
    merged_expression_table_normalized_all<-rbind(merged_expression_table_normalized_stage_I,merged_expression_table_normalized_stage_II,merged_expression_table_normalized_stage_III)   
    
    # Merge the three stages
    Interactomes_GC3_T2_selected_all<-rbind(Interactomes_GC3_T2_merged_Stage_I,Interactomes_GC3_T2_merged_Stage_II,Interactomes_GC3_T2_merged_Stage_III)
    #########################################################################################################################################        
    merged_expression_pca_normalized_all <- prcomp(merged_expression_table_normalized_all[,c("tpm","T2","Conections")], scale. = TRUE)
    Interactomes_GC3_T2_pca_all <- prcomp(Interactomes_GC3_T2_selected_all[,c("AveExp","T2","Conections")], scale. = TRUE)

    p1<-autoplot(merged_expression_pca_normalized_all, data = merged_expression_table_normalized_all, colour = 'Stage')+ theme_bw() + ggtitle("All points")
    p2<-autoplot(Interactomes_GC3_T2_pca_all, data = Interactomes_GC3_T2_selected_all, colour = 'Stage')+ theme_bw() + ggtitle("Average per genes")

    # FindClusters_resolution               
    png(filename=paste(output_dir,"PCAs_merged_average_",normalization_scheme,"_Stage_all.png",sep=""), width = 20, height = 15, res=600, units = "cm")  
            ggarrange(p1,p2,nrow = 1,ncol = 3, common.legend = TRUE, legend="bottom")
    dev.off()

    # FindClusters_resolution               
    png(filename=paste(output_dir,"PCA_average_",normalization_scheme,"_Stage_all.png",sep=""), width = 15, height = 15, res=600, units = "cm")  
            autoplot(Interactomes_GC3_T2_pca_all, data = Interactomes_GC3_T2_selected_all, colour = 'Stage')+ theme_bw() + ggtitle("Average per genes")
    dev.off()
  
  }
