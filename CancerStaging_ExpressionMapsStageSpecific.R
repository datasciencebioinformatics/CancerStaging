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
Interactomes_GC3_T2_merged<-Interactomes_GC3_T2_merged[Interactomes_GC3_T2_merged$Conections <=1600,]
limit_expr=10000
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

    #######################################################################################################################################                                                                                                                                     #
    unique_stage_I  =intersect(setdiff(selected_genes_Stage_I_gene, c(selected_genes_Stage_II_gene,selected_genes_Stage_III_gene)),selected_genes_Stage_I_gene)
    unique_stage_II =intersect(setdiff(selected_genes_Stage_II_gene, c(selected_genes_Stage_I_gene,selected_genes_Stage_III_gene)),selected_genes_Stage_II_gene)
    unique_stage_III=intersect(setdiff(selected_genes_Stage_III_gene, c(selected_genes_Stage_I_gene,selected_genes_Stage_II_gene)),selected_genes_Stage_III_gene)
    #########################################################################################################################################        
    selected_genes_Stage_I_gene      <- genes_stages_I
    selected_genes_Stage_II_gene     <- genes_stages_II
    selected_genes_Stage_III_gene    <- genes_stages_III
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
 
    merged_expression_table_normalized_stage_I  <-merged_expression_table_normalized_stage_I[merged_expression_table_normalized_stage_I$Expr>0 & merged_expression_table_normalized_stage_I$Expr<limit_expr,]
    merged_expression_table_normalized_stage_II <-merged_expression_table_normalized_stage_II[merged_expression_table_normalized_stage_II$Expr>0 & merged_expression_table_normalized_stage_II$Expr<limit_expr,]
    merged_expression_table_normalized_stage_III <-merged_expression_table_normalized_stage_III[merged_expression_table_normalized_stage_III$Expr>0  & merged_expression_table_normalized_stage_III$Expr<limit_expr,]
    merged_expression_table_normalized_all_stages <-merged_expression_table_normalized_all_stages[merged_expression_table_normalized_all_stages$Expr>0  & merged_expression_table_normalized_all_stages$Expr<limit_expr,]    

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

    merged_expression_table_normalized_stage_overlapping$Exp_z_score        <-  calculate_z(merged_expression_table_normalized_stage_overlapping$Expr,        mean(merged_expression_table_normalized_stage_overlapping$Expr, na.rm = TRUE),sd(merged_expression_table_normalized_stage_overlapping$Expr, na.rm = TRUE))  
    merged_expression_table_normalized_stage_overlapping$Conections_z_score <-  calculate_z(merged_expression_table_normalized_stage_overlapping$Conections,  mean(merged_expression_table_normalized_stage_overlapping$Conections, na.rm = TRUE),sd(merged_expression_table_normalized_stage_overlapping$Conections, na.rm = TRUE))    
    merged_expression_table_normalized_stage_overlapping$T2_z_score         <-  calculate_z(merged_expression_table_normalized_stage_overlapping$T2,          mean(merged_expression_table_normalized_stage_overlapping$T2, na.rm = TRUE),sd(merged_expression_table_normalized_stage_overlapping$T2, na.rm = TRUE))        

    merged_expression_table_normalized_all_stages<-rbind(merged_expression_table_normalized_stage_I,merged_expression_table_normalized_stage_II,merged_expression_table_normalized_stage_III, merged_expression_table_normalized_stage_overlapping)  
  
        
    m1<-ggplot(merged_expression_table_normalized_stage_I, aes(Conections, T2, z = Expr))  + geom_point(aes(colour=Expr)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))+ geom_contour()  + xlim(0, 50) + ylim(0, 60)        + geom_vline(xintercept=mean(merged_expression_table_normalized_stage_I$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_I$T2), linetype="dashed", color = "red")     + geom_vline(xintercept=median(merged_expression_table_normalized_stage_I$Conections), linetype="dashed", color = "yellow")   +  geom_hline(yintercept=median(merged_expression_table_normalized_stage_I$T2), linetype="dashed", color =   "yellow")   
    m2<-ggplot(merged_expression_table_normalized_stage_II, aes(Conections, T2, z = Expr))  + geom_point(aes(colour=Expr)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,   ": All points Stage II Expr. ",sep=""))+ geom_contour()  + xlim(0, 50) + ylim(0, 60)       + geom_vline(xintercept=mean(merged_expression_table_normalized_stage_II$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_II$T2), linetype="dashed", color = "red")   + geom_vline(xintercept=median(merged_expression_table_normalized_stage_II$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(merged_expression_table_normalized_stage_II$T2), linetype="dashed", color =  "yellow")
    m3<-ggplot(merged_expression_table_normalized_stage_III, aes(Conections, T2, z = Expr))  + geom_point(aes(colour=Expr)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,  ": All points Stage III Expr. ",sep=""))+ geom_contour()  + xlim(0, 50) + ylim(0, 60)      + geom_vline(xintercept=mean(merged_expression_table_normalized_stage_III$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_III$T2), linetype="dashed", color = "red") + geom_vline(xintercept=median(merged_expression_table_normalized_stage_III$Conections), linetype="dashed", color = "yellow") +  geom_hline(yintercept=median(merged_expression_table_normalized_stage_III$T2), linetype="dashed", color = "yellow")    
    
    m4 <- ggplot(merged_expression_table_normalized_all_stages, aes(x=Stages, y=T2)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()           +  geom_hline(yintercept=median(merged_expression_table_normalized_all_stages$T2), linetype="dashed", color = "red")    + ggtitle(paste("T2 All points: ", normalization_scheme,sep="")) +  stat_compare_means(comparisons = my_comparisons, method = "t.test")
    m5 <- ggplot(merged_expression_table_normalized_all_stages, aes(x=Stages, y=Conections)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()   +  geom_hline(yintercept=median(merged_expression_table_normalized_all_stages$Conections), linetype="dashed", color = "red")    + ggtitle(paste("Connectivity All points: ", normalization_scheme,sep="")) +  stat_compare_means(comparisons = my_comparisons, method = "t.test")
    m6 <- ggplot(merged_expression_table_normalized_all_stages, aes(x=Stages, y=Expr)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()         +      geom_hline(yintercept=median(merged_expression_table_normalized_all_stages$Expr), linetype="dashed", color = "red")    + ggtitle(paste("Expr. All points: ", normalization_scheme,sep="")) +  stat_compare_means(comparisons = my_comparisons, method = "t.test")

    m7<-ggplot(merged_expression_table_normalized_stage_I, aes(Conections, T2, z = Expr))  + geom_point(aes(colour=Expr))   +  theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))+ geom_contour()  + xlim(0, 50) + ylim(0, 60)        + geom_vline(xintercept=mean(merged_expression_table_normalized_stage_I$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_I$T2), linetype="dashed", color = "red")     + geom_vline(xintercept=median(merged_expression_table_normalized_stage_I$Conections), linetype="dashed", color = "yellow")   +  geom_hline(yintercept=median(merged_expression_table_normalized_stage_I$T2), linetype="dashed", color =   "yellow")   
    m8<-ggplot(merged_expression_table_normalized_stage_II, aes(Conections, T2, z = Expr))  + geom_point(aes(colour=Expr))  +  theme_bw() + ggtitle(paste(normalization_scheme,   ": All points Stage II Expr. ",sep=""))+ geom_contour()  + xlim(0, 50) + ylim(0, 60)       + geom_vline(xintercept=mean(merged_expression_table_normalized_stage_II$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_II$T2), linetype="dashed", color = "red")   + geom_vline(xintercept=median(merged_expression_table_normalized_stage_II$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(merged_expression_table_normalized_stage_II$T2), linetype="dashed", color =  "yellow")
    m9<-ggplot(merged_expression_table_normalized_stage_III, aes(Conections, T2, z = Expr))  + geom_point(aes(colour=Expr)) + theme_bw() + ggtitle(paste(normalization_scheme,  ": All points Stage III Expr. ",sep=""))+ geom_contour()  + xlim(0, 50) + ylim(0, 60)      + geom_vline(xintercept=mean(merged_expression_table_normalized_stage_III$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_III$T2), linetype="dashed", color = "red") + geom_vline(xintercept=median(merged_expression_table_normalized_stage_III$Conections), linetype="dashed", color = "yellow") +  geom_hline(yintercept=median(merged_expression_table_normalized_stage_III$T2), linetype="dashed", color = "yellow")      
  
    # FindClusters_resolution               
    png(filename=paste(output_dir,"countour_T2_Coonections_melt_",normalization_scheme,"_Stage_all.png",sep=""), width = 30, height = 50, res=600, units = "cm")  
            ggarrange(m7, m8, m9, m1,m2,m3,m4,m5,m6, nrow = 3,ncol = 3, common.legend = TRUE, legend="bottom")
    dev.off()  
      
    # Filter up Average expression greater than zero        
    m1<-ggplot(merged_expression_table_normalized_stage_I, aes(Conections_z_score, T2_z_score, z = Exp_z_score))  + geom_point(aes(colour=Exp_z_score)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,    ": Stage I Z-score",sep=""))+ geom_contour()        + geom_vline(xintercept=mean(merged_expression_table_normalized_stage_I$Conections_z_score), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_I$T2_z_score), linetype="dashed", color = "red")        + geom_vline(xintercept=median(merged_expression_table_normalized_stage_I$T2_z_score), linetype="dashed", color = "yellow") +  geom_hline(yintercept=median(merged_expression_table_normalized_stage_I$Conections_z_score), linetype="dashed", color = "yellow")+xlim(0,4)+ylim(-2,2)
    m2<-ggplot(merged_expression_table_normalized_stage_II, aes(Conections_z_score, T2_z_score, z = Exp_z_score))  + geom_point(aes(colour=Exp_z_score)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,    ": Stage II Z-score",sep=""))+ geom_contour()      + geom_vline(xintercept=mean(merged_expression_table_normalized_stage_II$Conections_z_score), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_II$T2_z_score), linetype="dashed", color = "red")      + geom_vline(xintercept=median(merged_expression_table_normalized_stage_II$T2_z_score), linetype="dashed", color = "yellow") +  geom_hline(yintercept=median(merged_expression_table_normalized_stage_II$Conections_z_score), linetype="dashed", color = "yellow")+xlim(0,4)+ylim(-2,2)      
    m3<-ggplot(merged_expression_table_normalized_stage_III, aes(Conections_z_score, T2_z_score, z = Exp_z_score))  + geom_point(aes(colour=Exp_z_score)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,    ": Stage III Z-score",sep=""))+ geom_contour()    + geom_vline(xintercept=mean(merged_expression_table_normalized_stage_III$Conections_z_score), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(merged_expression_table_normalized_stage_III$T2_z_score), linetype="dashed", color = "red")     + geom_vline(xintercept=median(merged_expression_table_normalized_stage_III$T2_z_score), linetype="dashed", color = "yellow") +  geom_hline(yintercept=median(merged_expression_table_normalized_stage_III$Conections_z_score), linetype="dashed", color = "yellow")+xlim(0,4)+ylim(-2,2)  

    m4 <- ggplot(merged_expression_table_normalized_all_stages, aes(x=Stages, y=T2_z_score)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()           +  geom_hline(yintercept=mean(merged_expression_table_normalized_all_stages$T2_z_score), linetype="dashed", color = "red")    + ggtitle(paste("T2 z-score: ", normalization_scheme,sep="")) 
    m5 <- ggplot(merged_expression_table_normalized_all_stages, aes(x=Stages, y=Conections_z_score)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()   +  geom_hline(yintercept=mean(merged_expression_table_normalized_all_stages$Conections_z_score), linetype="dashed", color = "red")    + ggtitle(paste("Connectivity z-score: ", normalization_scheme,sep="")) 
    m6 <- ggplot(merged_expression_table_normalized_all_stages, aes(x=Stages, y=Exp_z_score)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()          +  geom_hline(yintercept=mean(merged_expression_table_normalized_all_stages$Exp_z_score), linetype="dashed", color = "red")    + ggtitle(paste("Expr. Z-score: ", normalization_scheme,sep="")) 
      
    # FindClusters_resolution               
    png(filename=paste(output_dir,"countour_T2_Coonections_zscore_",normalization_scheme,"_Stage_all.png",sep=""), width = 30, height = 25, res=600, units = "cm")  
            ggarrange(m1,m2,m3,m4,m5,m6, nrow = 2,ncol = 3, common.legend = TRUE, legend="bottom")
    dev.off()  
    #########################################################################################################################################    
}
