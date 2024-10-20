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
################################################################################################################################
saveRDS(object = Interactomes_GC3_T2_merged, file = paste(output_dir,"Interactomes_GC3_T2_merged.rds",sep=""))                 #
saveRDS(object = normalization_schemes    , file = paste(output_dir,"normalization_schemes.rds",sep=""))                       #
saveRDS(object = expression_table_normalized    , file = paste(output_dir,"expression_table_normalized.rds",sep=""))           # 
################################################################################################################################
# Interactomes_GC3_T2.csv file has 15650 entries. The number of annotated genes with gene length geneLength_ENTREZID_ENSEMBL is 14609. Among these, 14726 are common to Interactomes_GC3_T2 and geneLength_ENTREZID_ENSEMBL and will be used to create the maps. 
# Consitency - check filters meticulously.
# FPKM, TPM  - take these as robust.
# Paramter to set the normalization_scheme

# Interactomes_GC3_T2.csv file has 15650 entries. The number of annotated genes with gene length geneLength_ENTREZID_ENSEMBL is 14609. Among these, 14726 are common to Interactomes_GC3_T2 and geneLength_ENTREZID_ENSEMBL and will be used to create the maps. 
# Consitency - check filters meticulously.
# FPKM, TPM  - take these as robust.
# Paramter to set the normalization_scheme
normalization_schemes<-c("tpm","fpkm","tmm","rpkm")
normalization_schemes<-c("tpm","tmm")





# For each normlization normalization_scheme
for (normalization_scheme in normalization_schemes)
{     
    # genes_stages_I
    genes_stages_I    <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_I",".tsv",sep=""), sep = '\t', header = TRUE)$gene #
    genes_stages_II   <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_II",".tsv",sep=""), sep = '\t', header = TRUE)$gene #
    genes_stages_III  <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_III",".tsv",sep=""), sep = '\t', header = TRUE)$gene #
  
    # Take also expression data from the normalization scheme set by "normalization_scheme"
    expression_table_normalized<-expression_table_normalized[[normalization_scheme]]

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
    m1<-ggplot(Interactomes_GC3_T2_merged_Stage_I, aes(Conections, T2, z = AveExp))   + geom_density_2d_filled(alpha = 0.5) + theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))+  geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_I$T2), linetype="dashed", color = "red")   +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_I$T2), linetype="dashed", color = "red")     + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_I$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_I$T2), linetype="dashed", color =   "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE)) + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))
    m2<-ggplot(Interactomes_GC3_T2_merged_Stage_II, aes(Conections, T2, z = AveExp))  + geom_density_2d_filled(alpha = 0.5) + theme_bw() + ggtitle(paste(normalization_scheme,   ": All points Stage II Expr. ",sep=""))+  geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_II$T2), linetype="dashed", color = "red")  +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_II$T2), linetype="dashed", color = "red")   + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_II$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_II$T2), linetype="dashed", color =  "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))
    m3<-ggplot(Interactomes_GC3_T2_merged_Stage_III, aes(Conections, T2, z = AveExp)) + geom_density_2d_filled(alpha = 0.5) + theme_bw() + ggtitle(paste(normalization_scheme,  ": All points Stage III Expr. ",sep=""))+  geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_III$T2), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_III$T2), linetype="dashed", color = "red") + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_III$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_III$T2), linetype="dashed", color = "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))
  
    m4 <- ggplot(Interactomes_GC3_T2_merged_all, aes(x=Stages, y=T2)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()           +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_all$T2), linetype="dashed", color = "red")    + ggtitle(paste("T2 All points: ", normalization_scheme,sep=""))                   +  stat_compare_means(comparisons = my_comparisons, method = "t.test") 
    m5 <- ggplot(Interactomes_GC3_T2_merged_all, aes(x=Stages, y=Conections)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()   +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_all$Conections), linetype="dashed", color = "red")    + ggtitle(paste("Connectivity All points: ", normalization_scheme,sep="")) +  stat_compare_means(comparisons = my_comparisons, method = "t.test") 
    m6 <- ggplot(Interactomes_GC3_T2_merged_all, aes(x=Stages, y=AveExp)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()         +      geom_hline(yintercept=median(Interactomes_GC3_T2_merged_all$Expr), linetype="dashed", color = "red")    + ggtitle(paste("Expr. All points: ", normalization_scheme,sep=""))          +  stat_compare_means(comparisons = my_comparisons, method = "t.test")       
  
    m7<-ggplot(Interactomes_GC3_T2_merged_Stage_I, aes(Conections, T2, z = AveExp))    + geom_point(aes(colour=AveExp))   +  theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))+  xlim(0, 50)       + geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_I$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_I$T2), linetype="dashed", color = "red")     + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_I$Conections), linetype="dashed", color = "yellow") +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_I$T2), linetype="dashed", color =   "yellow")   + ylim(10,40) + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(0, 40, by = 10), limits = c(10, 40))
    m8<-ggplot(Interactomes_GC3_T2_merged_Stage_II, aes(Conections, T2, z = AveExp))   + geom_point(aes(colour=AveExp))  +  theme_bw() + ggtitle(paste(normalization_scheme,   ": All points Stage II Expr. ",sep=""))+ xlim(0, 50) + geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_II$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_II$T2), linetype="dashed", color = "red")   + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_II$Conections), linetype="dashed", color = "yellow") +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_II$T2), linetype="dashed", color =  "yellow")   + ylim(10,40) + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)       + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(0, 40, by = 10), limits = c(10, 40))
    m9<-ggplot(Interactomes_GC3_T2_merged_Stage_III, aes(Conections, T2, z = AveExp))  + geom_point(aes(colour=AveExp)) + theme_bw() + ggtitle(paste(normalization_scheme,  ": All points Stage III Expr. ",sep=""))+ xlim(0, 50) +  geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_III$Conections), linetype="dashed", color = "red")  +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_III$T2), linetype="dashed", color = "red") + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_III$Conections), linetype="dashed", color = "yellow") +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_III$T2), linetype="dashed", color = "yellow")   + ylim(10,40) +  guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE)) + ylim(10,40)  + xlim(0, 50)     + scale_y_continuous(minor_breaks = seq(10, 40, by = 1), breaks = seq(0, 40, by = 10), limits = c(10, 40))
  
    m10<-ggplot(Interactomes_GC3_T2_merged_Stage_I, aes(Conections, T2, z = AveExp))   +  theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))   + geom_density_2d(bin=10)         + geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_I$Conections), linetype="dashed", color = "red")   +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_I$T2), linetype="dashed", color = "red")     + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_I$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_I$T2), linetype="dashed", color =   "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE)) + ylim(10,40)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(0, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))
    m11<-ggplot(Interactomes_GC3_T2_merged_Stage_II, aes(Conections, T2, z = AveExp))  +  theme_bw() + ggtitle(paste(normalization_scheme,   ": All points Stage II Expr. ",sep=""))   + geom_density_2d(bin=10)        + geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_II$Conections), linetype="dashed", color = "red")  +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_II$T2), linetype="dashed", color = "red")   + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_II$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_II$T2), linetype="dashed", color =  "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(0, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))
    m12<-ggplot(Interactomes_GC3_T2_merged_Stage_III, aes(Conections, T2, z = AveExp)) +  theme_bw() + ggtitle(paste(normalization_scheme,  ": All points Stage III Expr. ",sep=""))   + geom_density_2d(bin=10)       + geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_III$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_III$T2), linetype="dashed", color = "red") + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_III$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_III$T2), linetype="dashed", color = "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(10,40)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(0, 40, by = 1), breaks = seq(10, 40, by = 10), limits = c(10, 40))  
      
    # FindClusters_resolution               
    png(filename=paste(output_dir,"countour_T2_Coonections_melt_",normalization_scheme,"_",TCGA_project,"_Stage_all_T2.png",sep=""), width = 30, height = 50, res=600, units = "cm")  
            plot<-ggarrange(m4, m5, m6, m7, m8, m9, m10, m11,m12, m1,m2,m3, nrow = 4,ncol = 3, common.legend = TRUE, legend="bottom")
            print(annotate_figure(plot, top = text_grob(TCGA_project, face = "bold", size = 14)))
    dev.off()

    #########################################################################################################################################
    m1<-ggplot(Interactomes_GC3_T2_merged_Stage_I, aes(Conections, GC3, z = AveExp))   + geom_density_2d_filled(alpha = 0.5) + theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))+  geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_I$Conections), linetype="dashed", color = "red")   +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_I$GC3), linetype="dashed", color = "red")     + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_I$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_I$GC3), linetype="dashed", color =   "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE)) + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85))
    m2<-ggplot(Interactomes_GC3_T2_merged_Stage_II, aes(Conections, GC3, z = AveExp))  + geom_density_2d_filled(alpha = 0.5) + theme_bw() + ggtitle(paste(normalization_scheme,   ": All points Stage II Expr. ",sep=""))+  geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_II$Conections), linetype="dashed", color = "red")  +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_II$GC3), linetype="dashed", color = "red")   + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_II$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_II$GC3), linetype="dashed", color =  "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(25,85)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85))
    m3<-ggplot(Interactomes_GC3_T2_merged_Stage_III, aes(Conections, GC3, z = AveExp)) + geom_density_2d_filled(alpha = 0.5) + theme_bw() + ggtitle(paste(normalization_scheme,  ": All points Stage III Expr. ",sep=""))+  geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_III$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_III$GC3), linetype="dashed", color = "red") + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_III$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_III$GC3), linetype="dashed", color = "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(25,85)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85))
  
    m4 <- ggplot(Interactomes_GC3_T2_merged_all, aes(x=Stages, y=GC3)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()           +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_all$GC3), linetype="dashed", color = "red")    + ggtitle(paste("T2 All points: ", normalization_scheme,sep=""))                   +  stat_compare_means(comparisons = my_comparisons, method = "t.test") 
    m5 <- ggplot(Interactomes_GC3_T2_merged_all, aes(x=Stages, y=Conections)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()   +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_all$Conections), linetype="dashed", color = "red")    + ggtitle(paste("Connectivity All points: ", normalization_scheme,sep="")) +  stat_compare_means(comparisons = my_comparisons, method = "t.test") 
    m6 <- ggplot(Interactomes_GC3_T2_merged_all, aes(x=Stages, y=AveExp)) +  geom_violin(trim=FALSE) + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2,color="red" )+ theme_bw()         +      geom_hline(yintercept=median(Interactomes_GC3_T2_merged_all$Expr), linetype="dashed", color = "red")    + ggtitle(paste("Expr. All points: ", normalization_scheme,sep=""))          +  stat_compare_means(comparisons = my_comparisons, method = "t.test")       
  
    m7<-ggplot(Interactomes_GC3_T2_merged_Stage_I, aes(Conections, GC3, z = AveExp))    + geom_point(aes(colour=AveExp))   +  theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))+  xlim(0, 50)       + geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_I$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_I$GC3), linetype="dashed", color = "red")     + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_I$Conections), linetype="dashed", color = "yellow") +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_I$GC3), linetype="dashed", color =   "yellow")   + ylim(25,85) + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(25,85)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(0, 40, by = 10), limits = c(25, 85))
    m8<-ggplot(Interactomes_GC3_T2_merged_Stage_II, aes(Conections, GC3, z = AveExp))   + geom_point(aes(colour=AveExp))  +  theme_bw() + ggtitle(paste(normalization_scheme,   ": All points Stage II Expr. ",sep=""))+ xlim(0, 50) + geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_II$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_II$GC3), linetype="dashed", color = "red")   + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_II$Conections), linetype="dashed", color = "yellow") +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_II$GC3), linetype="dashed", color =  "yellow")   + ylim(25,85) + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(25,85)  + xlim(0, 50)       + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(0, 40, by = 10), limits = c(25, 85))
    m9<-ggplot(Interactomes_GC3_T2_merged_Stage_III, aes(Conections, GC3, z = AveExp))  + geom_point(aes(colour=AveExp)) + theme_bw() + ggtitle(paste(normalization_scheme,  ": All points Stage III Expr. ",sep=""))+ xlim(0, 50) +  geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_III$Conections), linetype="dashed", color = "red")  +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_III$GC3), linetype="dashed", color = "red") + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_III$Conections), linetype="dashed", color = "yellow") +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_III$GC3), linetype="dashed", color = "yellow")   + ylim(25,85) +  guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE)) + ylim(25,85)  + xlim(0, 50)     + scale_y_continuous(minor_breaks = seq(25, 85, by = 1), breaks = seq(0, 40, by = 10), limits = c(25, 85))
  
    m10<-ggplot(Interactomes_GC3_T2_merged_Stage_I, aes(Conections, GC3, z = AveExp))   +  theme_bw() + ggtitle(paste(normalization_scheme,    ": All points Stage I Expr. ",sep=""))   + geom_density_2d(bin=10)         + geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_I$Conections), linetype="dashed", color = "red")   +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_I$GC3), linetype="dashed", color = "red")     + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_I$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_I$GC3), linetype="dashed", color =   "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE)) + ylim(25,85)  + xlim(0, 50)   + scale_y_continuous(minor_breaks = seq(0, 40, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85))
    m11<-ggplot(Interactomes_GC3_T2_merged_Stage_II, aes(Conections, GC3, z = AveExp))  +  theme_bw() + ggtitle(paste(normalization_scheme,   ": All points Stage II Expr. ",sep=""))   + geom_density_2d(bin=10)        + geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_II$Conections), linetype="dashed", color = "red")  +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_II$GC3), linetype="dashed", color = "red")   + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_II$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_II$GC3), linetype="dashed", color =  "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(25,85)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(0, 40, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85))
    m12<-ggplot(Interactomes_GC3_T2_merged_Stage_III, aes(Conections, GC3, z = AveExp)) +  theme_bw() + ggtitle(paste(normalization_scheme,  ": All points Stage III Expr. ",sep=""))   + geom_density_2d(bin=10)       + geom_vline(xintercept=mean(Interactomes_GC3_T2_merged_Stage_III$Conections), linetype="dashed", color = "red") +  geom_hline(yintercept=mean(Interactomes_GC3_T2_merged_Stage_III$GC3), linetype="dashed", color = "red") + geom_vline(xintercept=median(Interactomes_GC3_T2_merged_Stage_III$Conections), linetype="dashed", color = "yellow")  +  geom_hline(yintercept=median(Interactomes_GC3_T2_merged_Stage_III$GC3), linetype="dashed", color = "yellow")   + guides(x = guide_axis(minor.ticks = TRUE),y = guide_axis(minor.ticks = TRUE))  + ylim(25,85)  + xlim(0, 50)  + scale_y_continuous(minor_breaks = seq(0, 40, by = 1), breaks = seq(25, 85, by = 10), limits = c(25, 85))  
      
    # FindClusters_resolution               
    png(filename=paste(output_dir,"countour_T2_Coonections_melt_",normalization_scheme,"_",TCGA_project,"_Stage_all_GC3.png",sep=""), width = 30, height = 50, res=600, units = "cm")  
            plot<-ggarrange(m4, m5, m6, m7, m8, m9, m10, m11,m12, m1,m2,m3, nrow = 4,ncol = 3, common.legend = TRUE, legend="bottom")
            print(annotate_figure(plot, top = text_grob(TCGA_project, face = "bold", size = 14)))
    dev.off()  

    h1<-ggplot(unique(Interactomes_GC3_T2_merged_all[,c("ENSEMBL","T2","GC3","Stages")]), aes(x=T2, color=Stages)) +  geom_histogram(fill="white", alpha=0.5, position="identity", bins=20) + theme_bw()     + xlim(0, 50)
    h2<-ggplot(unique(Interactomes_GC3_T2_merged_all[,c("ENSEMBL","T2","GC3","Stages")]), aes(x=GC3, color=Stages)) +  geom_histogram(fill="white", alpha=0.5, position="identity", bins=20) + theme_bw()    + xlim(0, 100)
    h3<-ggplot(unique(Interactomes_GC3_T2_merged_all[,c("ENSEMBL","T2","AveExp","Stages")]), aes(x=AveExp, color=Stages)) +  geom_histogram(fill="white", alpha=0.5, position="identity", bins=20) + theme_bw()  + xlim(0, 10000)
    
    # FindClusters_resolution               
    png(filename=paste(output_dir,"countour_T2_Coonections_melt_",normalization_scheme,"_",TCGA_project,"_histogram.png",sep=""), width = 25, height = 15, res=600, units = "cm")            
          plot<-ggarrange(h1, h2, h3, nrow = 1, common.legend = TRUE, legend="bottom") 
          annotate_figure(plot, top = text_grob(TCGA_project, face = "bold", size = 14))  
    dev.off()    
}



