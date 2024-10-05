#######################################################################################################################
source("/home/felipe/Documents/github/CancerStaging/CancerStaging_SetupAllParamters.R")                               #
source("/home/felipe/Documents/github/CancerStaging/CancerStaging_LoadRPackages.R")                                   #
#######################################################################################################################
# saveRDS                                                                                                             #
#saveRDS(object = df_reads_count_all_projects, file = paste(output_dir,"df_reads_count_all_projects.rds",sep=""))      #
#saveRDS(object = geneLength_ENTREZID_ENSEMBL, file = paste(output_dir,"geneLength_ENTREZID_ENSEMBL.rds",sep=""))      #
#saveRDS(object = Interactomes_GC3_T2_merged, file = paste(output_dir,"Interactomes_GC3_T2_merged.rds",sep=""))        #
                                                                                                                      #
# Restore the object                                                                                                  #
geneLength_ENTREZID_ENSEMBL<-readRDS(file = paste(output_dir,"geneLength_ENTREZID_ENSEMBL.rds",sep=""))               #
Interactomes_GC3_T2_merged <-readRDS(file = paste(output_dir,"Interactomes_GC3_T2_merged.rds",sep=""))                #
df_reads_count_all_projects<-readRDS(file = paste(output_dir,"df_reads_count_all_projects.rds",sep=""))               ###########################
merged_data_patient_info   <-read.table(file = "/home/felipe/Documents/Cancer_staging/merged_data_patient_info.tsv", sep = '\t', header = TRUE) #
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
  
####################################################################################################################################################
# Interactomes_GC3_T2.csv file has 15650 entries. The number of annotated genes with gene length geneLength_ENTREZID_ENSEMBL is 14609. Among these, 14726 are common to Interactomes_GC3_T2 and geneLength_ENTREZID_ENSEMBL and will be used to create the maps. 
# Consitency - check filters meticulously.
# FPKM, TPM  - take these as robust.
# Paramter to set the normalization_scheme
normalization_schemes<-c("tpm","fpkm","tmm","rpkm","tpm_calc")

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
    
    # Set colnames
    colnames(merged_expression_table_normalized_stage_I)[6]<-normalization_scheme
    colnames(merged_expression_table_normalized_stage_II)[6]<-normalization_scheme
    colnames(merged_expression_table_normalized_stage_III)[6]<-normalization_scheme
    ####################################################################################################################################################
    # I have two tables to be use.                                                                                                                     #
    # First , a table with with "T2", "GC3", "Conections", "ENSEMBL" "AveExp" per gene                                                                 #
    # Interactomes_GC3_T2_merged                                                                                                                       #
    # Second, a table with with "T2", "GC3", "Conections", "ENSEMBL" "Exp"    per patient                                                              #
    # melt_expression_interactomes                                                                                                                     #
    #############################################################################################################################################################
    Interactomes_GC3_T2_selected_Stage_I                       <-merged_expression_table_normalized_stage_I[,c("T2",normalization_scheme,"Conections","GC3")]   #
    Interactomes_GC3_T2_selected_Stage_II                      <-merged_expression_table_normalized_stage_II[,c("T2",normalization_scheme,"Conections","GC3")]  #
    Interactomes_GC3_T2_selected_Stage_III                     <-merged_expression_table_normalized_stage_III[,c("T2",normalization_scheme,"Conections","GC3")] #
                                                                                                                                                                #
    # FindClusters_resolution          #                                                                                                                        #
    png(filename=paste(output_dir,"scatterplot3d_avg_",normalization_scheme,"_stage_I.png",sep=""), width = 24, height = 24, res=600, units = "cm")             #
            scatterplot3d(Interactomes_GC3_T2_selected_Stage_I[,c("T2","Conections",normalization_scheme)], pch = 16)                                           #
    dev.off()                                                                                                                                                   #
    png(filename=paste(output_dir,"scatterplot3d_avg_",normalization_scheme,"_stage_II.png",sep=""), width = 24, height = 24, res=600, units = "cm")            #
            scatterplot3d(Interactomes_GC3_T2_selected_Stage_II[,c("T2","Conections",normalization_scheme)], pch = 16)                                          #
    dev.off()                                                                                                                                                   #
    png(filename=paste(output_dir,"scatterplot3d_avg_",normalization_scheme,"_stage_III.png",sep=""), width = 24, height = 24, res=600, units = "cm")           #
            scatterplot3d(Interactomes_GC3_T2_selected_Stage_III[,c("T2","Conections",normalization_scheme)], pch = 16)                                         #
    dev.off()                                                                                                                                                   #
  ###############################################################################################################################################################
    # FindClusters_resolution          #          
    png(filename=paste(output_dir,"scatterplot3d_melt_",normalization_scheme,"_stage_I.png",sep=""), width = 24, height = 24, res=600, units = "cm")                                      #
            scatterplot3d(merged_expression_table_normalized_stage_I[,c("T2","Conections",normalization_scheme)], pch = 16)
    dev.off()
    # FindClusters_resolution          #               
    png(filename=paste(output_dir,"scatterplot3d_melt_",normalization_scheme,"_stage_II.png",sep=""), width = 24, height = 24, res=600, units = "cm")                                      #
            scatterplot3d(merged_expression_table_normalized_stage_II[,c("T2","Conections",normalization_scheme)], pch = 16)
    dev.off()  
      # FindClusters_resolution          #               
    png(filename=paste(output_dir,"scatterplot3d_melt_",normalization_scheme,"_stage_III.png",sep=""), width = 24, height = 24, res=600, units = "cm")                                      #
            scatterplot3d(merged_expression_table_normalized_stage_II[,c("T2","Conections",normalization_scheme)], pch = 16)
    dev.off()     
    ####################################################################################################################################################
    # Only Variable Labels on the outside (no axis labels)
    Interactomes_GC3_T2_mean <- ggpairs(Interactomes_GC3_T2_selected[,c("T2","GC3","AveExp","Conections")], axisLabels = "none")
    
    # FindClusters_resolution               
    png(filename=paste(output_dir,"correaltion_matrix_GC3_T2_mean_",normalization_scheme,".png",sep=""), width = 20, height = 20, res=600, units = "cm")  
            Interactomes_GC3_T2_mean
    dev.off()
    ####################################################################################################################################################
    # Only Variable Labels on the outside (no axis labels)
    Interactomes_GC3_T2_melt_Stage_I   <- ggpairs(merged_expression_table_normalized_stage_I[,c("T2","GC3",normalization_scheme,"Conections")], axisLabels = "none")
    Interactomes_GC3_T2_melt_Stage_II  <- ggpairs(merged_expression_table_normalized_stage_II[,c("T2","GC3",normalization_scheme,"Conections")], axisLabels = "none")
    Interactomes_GC3_T2_melt_Stage_III <- ggpairs(merged_expression_table_normalized_stage_III[,c("T2","GC3",normalization_scheme,"Conections")], axisLabels = "none")
    
    # FindClusters_resolution
    png(filename=paste(output_dir,"correaltion_matrix_melt_",normalization_scheme,"_stage_I.png",sep=""), width = 20, height = 20, res=600, units = "cm")  
            Interactomes_GC3_T2_melt_Stage_I
    dev.off()
    # FindClusters_resolution
    png(filename=paste(output_dir,"correaltion_matrix_melt_",normalization_scheme,"_stage_II.png",sep=""), width = 20, height = 20, res=600, units = "cm")  
            Interactomes_GC3_T2_melt_Stage_II
    dev.off()
    # FindClusters_resolution
    png(filename=paste(output_dir,"correaltion_matrix_melt_",normalization_scheme,"_stage_III.png",sep=""), width = 20, height = 20, res=600, units = "cm")  
            Interactomes_GC3_T2_melt_Stage_III
    dev.off()  
    ###########################################################################################################################################################
    # Only Variable Labels on the outside (no axis labels)
    Interactomes_GC3_T2_avg_Stage_I   <- ggpairs(Interactomes_GC3_T2_merged_Stage_I[,c("T2","GC3","AveExp","Conections")], axisLabels = "none")
    Interactomes_GC3_T2_avg_Stage_II  <- ggpairs(Interactomes_GC3_T2_merged_Stage_II[,c("T2","GC3","AveExp","Conections")], axisLabels = "none")
    Interactomes_GC3_T2_avg_Stage_III <- ggpairs(Interactomes_GC3_T2_merged_Stage_III[,c("T2","GC3","AveExp","Conections")], axisLabels = "none")
    
    # FindClusters_resolution
    png(filename=paste(output_dir,"correaltion_matrix_melt_",normalization_scheme,"_stage_I.png",sep=""), width = 20, height = 20, res=600, units = "cm")  
            Interactomes_GC3_T2_avg_Stage_I
    dev.off()
    # FindClusters_resolution
    png(filename=paste(output_dir,"correaltion_matrix_melt_",normalization_scheme,"_stage_II.png",sep=""), width = 20, height = 20, res=600, units = "cm")  
            Interactomes_GC3_T2_avg_Stage_II
    dev.off()
    # FindClusters_resolution
    png(filename=paste(output_dir,"correaltion_matrix_melt_",normalization_scheme,"_stage_III.png",sep=""), width = 20, height = 20, res=600, units = "cm")  
            Interactomes_GC3_T2_avg_Stage_III
    dev.off()  
    ###########################################################################################################################################################
  
    # FindClusters_resolution               
    png(filename=paste(output_dir,"geom_contour_melt",normalization_scheme,".png",sep=""), width = 14, height = 14, res=600, units = "cm")  
            ggplot(melt_expression_interactomes, aes(T2, AveExp, z = Conections))  + geom_point(aes(colour=Conections))
    dev.off()
    ###########################################################################################################################################################
    # melt_expression_interactomes
    melt_expression_interactomes$Sample.Type<-merged_data_patient_info[match(melt_expression_interactomes$variable, merged_data_patient_info$sample_id, nomatch = NA_integer_, incomparables = NULL),"Sample.Type"]

    # Select collumns
    tp53_expresion<-melt_expression_interactomes[which(melt_expression_interactomes$ENSEMBL=="ENSG00000141510"),c("Sample.Type","T2","GC3",normalization_scheme,"variable","Conections")] #
    ###########################################################################################################################################################
    # Basic box plot
    p1_tp53 <- ggplot(tp53_expresion, aes(x=Sample.Type, y=normalization_scheme)) + geom_boxplot(notch = TRUE)
    p1_tp53 <- p1_tp53 + theme_bw() + stat_summary(fun.y=mean, geom="point", shape=18,size=3, color="red")
    
    # Basic box plot
    p2_tp53 <- ggplot(tp53_expresion, aes(x=Sample.Type, y=T2)) + geom_boxplot(notch = TRUE)
    p2_tp53 <- p2_tp53 + theme_bw() + stat_summary(fun.y=mean, geom="point", shape=18,size=3, color="red")
    
    # Basic box plot
    p3_tp53 <- ggplot(tp53_expresion, aes(x=Sample.Type, y=GC3)) + geom_boxplot(notch = TRUE)
    p3_tp53 <- p3_tp53 + theme_bw() + stat_summary(fun.y=mean, geom="point", shape=18,size=3, color="red")

    # Basic box plot
    p4_tp53 <- ggplot(tp53_expresion, aes(x=Sample.Type, y=Conections)) + geom_boxplot(notch = TRUE)
    p4_tp53 <- p4_tp53 + theme_bw() + stat_summary(fun.y=mean, geom="point", shape=18,size=3, color="red")  
        
    # Basic box plot
    p1_all <- ggplot(melt_expression_interactomes, aes(x=Sample.Type, y=normalization_scheme)) + geom_boxplot(notch = TRUE)
    p1_all <- p1_all + theme_bw() + stat_summary(fun.y=mean, geom="point", shape=18,size=3, color="red")
    
    # Basic box plot
    p2_all <- ggplot(melt_expression_interactomes, aes(x=Sample.Type, y=T2)) + geom_boxplot(notch = TRUE)
    p2_all <- p2_all + theme_bw() + stat_summary(fun.y=mean, geom="point", shape=18,size=3, color="red")
    
    # Basic box plot
    p3_all <- ggplot(melt_expression_interactomes, aes(x=Sample.Type, y=GC3)) + geom_boxplot(notch = TRUE)
    p3_all <- p3_all + theme_bw() + stat_summary(fun.y=mean, geom="point", shape=18,size=3, color="red")

    # Basic box plot
    p4_all <- ggplot(melt_expression_interactomes, aes(x=Sample.Type, y=Conections)) + geom_boxplot(notch = TRUE)
    p4_all <- p4_all + theme_bw() + stat_summary(fun.y=mean, geom="point", shape=18,size=3, color="red")  
    
    grid_arrange_tp53<-grid.arrange(p1_tp53, p2_tp53, p3_tp53,p4_tp53, nrow = 1, top = "tp53 only")
    grid_arrange_all <-grid.arrange(p1_all, p2_all, p3_all, p4_all, nrow = 1, top = "all genes")    
    ###########################################################################################################################################################
    # FindClusters_resolution
    png(filename=paste(output_dir,"boxplot_GC3_T2_tp53_",normalization_scheme,".png",sep=""), width = 30, height = 20, res=600, units = "cm")  
      grid.arrange(grid_arrange_tp53, grid_arrange_all,  nrow = 2)
    dev.off()
    ###########################################################################################################################################################
    # Conections, T2, AvgExpression
    # Combine AvgExpression, Conections, T2
    # harmonic mean
    # z-core : Composite Scores
    # Implementing the Z score formula in R is quite straightforward. 
    # To reuse code, we will create a function called calculate_z using the mean and sd base functions to calculate Z. 
    # sd calculates the standard deviation in R.
    # weighted average
    # Z-score for AveExp_expression
    Interactomes_GC3_T2_selected$AveExp_z_score <- calculate_z(Interactomes_GC3_T2_selected$AveExp, 
                             mean(Interactomes_GC3_T2_selected$AveExp, na.rm = TRUE),
                             sd(Interactomes_GC3_T2_selected$AveExp, na.rm = TRUE))  
    # Z-score for AveExp_expression
    Interactomes_GC3_T2_selected$Conections_z_score <- calculate_z(Interactomes_GC3_T2_selected$Conections, 
                             mean(Interactomes_GC3_T2_selected$Conections, na.rm = TRUE),
                             sd(Interactomes_GC3_T2_selected$Conections, na.rm = TRUE)) 
    # Z-score for AveExp_expression
    Interactomes_GC3_T2_selected$T2_z_score <- calculate_z(Interactomes_GC3_T2_selected$T2, 
                             mean(Interactomes_GC3_T2_selected$T2, na.rm = TRUE),
                             sd(Interactomes_GC3_T2_selected$T2, na.rm = TRUE))   	    
    ###########################################################################################################################################################
    colnames(melt_expression_interactomes)[6]<-"Expr"

    # Conections, T2, AvgExpression
    # Combine AvgExpression, Conections, T2
    # harmonic mean
    # z-core : Composite Scores
    # Implementing the Z score formula in R is quite straightforward. 
    # To reuse code, we will create a function called calculate_z using the mean and sd base functions to calculate Z. 
    # sd calculates the standard deviation in R.
    # weighted average
    # Z-score for AveExp_expression
    melt_expression_interactomes$AveExp_z_score <- calculate_z(melt_expression_interactomes$AveExp, 
                             mean(melt_expression_interactomes$AveExp, na.rm = TRUE),
                             sd(melt_expression_interactomes$AveExp, na.rm = TRUE))  
    # Z-score for AveExp_expression
    melt_expression_interactomes$Conections_z_score <- calculate_z(melt_expression_interactomes$Conections, 
                             mean(melt_expression_interactomes$Conections, na.rm = TRUE),
                             sd(melt_expression_interactomes$Conections, na.rm = TRUE)) 
    # Z-score for AveExp_expression
    melt_expression_interactomes$T2_z_score <- calculate_z(melt_expression_interactomes$T2, 
                             mean(melt_expression_interactomes$T2, na.rm = TRUE),
                             sd(melt_expression_interactomes$T2, na.rm = TRUE))     
    #########################################################################################################################################
    # Filter up Average expression greater than zero
    Interactomes_GC3_T2_selected<-Interactomes_GC3_T2_selected[Interactomes_GC3_T2_selected$AveExp>0,]
    Interactomes_GC3_T2_selected<-Interactomes_GC3_T2_selected[Interactomes_GC3_T2_selected$Conections>0,]
    Interactomes_GC3_T2_selected<-Interactomes_GC3_T2_selected[Interactomes_GC3_T2_selected$T2>0,]
 
    m1<-ggplot(Interactomes_GC3_T2_selected, aes(Conections, T2, z = AveExp))  + geom_point(aes(colour=AveExp)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,  ": AveExpv vs .Conections : All points",sep=""))+ geom_contour()+ theme(legend.position="none")  # + theme(legend.position="none")
    m2<-ggplot(Interactomes_GC3_T2_selected, aes(Conections, T2, z = AveExp))  + geom_point(aes(colour=AveExp)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,  ": AveExpv vs .Conections : 0-100",sep=""))+ geom_contour()  + xlim(0, 100) + ylim(0, 100)     # + theme(legend.position="none")
    m3<-ggplot(Interactomes_GC3_T2_selected, aes(Conections, T2, z = AveExp))  + geom_point(aes(colour=AveExp)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,  ": AveExpv vs .Conections : 0-50",sep=""))+ geom_contour()  + xlim(0, 50) + ylim(0, 50)      #  + theme(legend.position="none")        
    m4<-ggplot(Interactomes_GC3_T2_selected, aes(Conections, T2, z = AveExp))  + geom_point(aes(colour=AveExp)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,  ": AveExpv vs .Conections : 0-25",sep=""))+ geom_contour()  + xlim(0, 25) + ylim(0, 25)      #  + theme(legend.position="none")          
    #########################################################################################################################################    
    # FindClusters_resolution               
    png(filename=paste(output_dir,"countour_T2_Coonections_Expression_",normalization_scheme,".png",sep=""), width = 25, height = 25, res=600, units = "cm")  
            ggarrange(m1,m2,m3,m4,nrow = 2,ncol = 2, common.legend = TRUE, legend="bottom")
    dev.off()  
    #########################################################################################################################################    
    #########################################################################################################################################
    # Change collumn id
    colnames(melt_expression_interactomes)[6]<-"Expr"
  
    # Filter up Average expression greater than zero
    melt_expression_interactomes<-melt_expression_interactomes[melt_expression_interactomes$Expr>0,]
    melt_expression_interactomes<-melt_expression_interactomes[melt_expression_interactomes$Conections>0,]
    melt_expression_interactomes<-melt_expression_interactomes[melt_expression_interactomes$T2>0,]

    m1<-ggplot(melt_expression_interactomes, aes(Conections, T2, z = Expr))  + geom_point(aes(colour=Expr)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,  ": Expv vs .Conections : All points",sep=""))+ geom_contour()+ theme(legend.position="none")  # + theme(legend.position="none")
    m2<-ggplot(melt_expression_interactomes, aes(Conections, T2, z = Expr))  + geom_point(aes(colour=Expr)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,  ": Expv vs .Conections : 0-100",sep=""))+ geom_contour()     + xlim(0, 100) + ylim(0, 100)     # + theme(legend.position="none")
    m3<-ggplot(melt_expression_interactomes, aes(Conections, T2, z = Expr))  + geom_point(aes(colour=Expr)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,  ": Expv vs .Conections : 0-50",sep=""))+ geom_contour()      + xlim(0, 50) + ylim(0, 50)      #  + theme(legend.position="none")        
    m4<-ggplot(melt_expression_interactomes, aes(Conections, T2, z = Expr))  + geom_point(aes(colour=Expr)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,  ": Expv vs .Conections : 0-25",sep=""))+ geom_contour()      + xlim(0, 25) + ylim(0, 25)      #  + theme(legend.position="none")          
    #########################################################################################################################################    
    # FindClusters_resolution               
    png(filename=paste(output_dir,"countour_T2_Coonections_AvgExpression_",normalization_scheme,".png",sep=""), width = 20, height = 20, res=600, units = "cm")  
            ggarrange(m1,m2,m3,m4,nrow = 2,ncol = 2, common.legend = TRUE, legend="bottom")
    dev.off()  
    #########################################################################################################################################      
    #########################################################################################################################################
    m1<-ggplot(Interactomes_GC3_T2_selected, aes(Conections_z_score, T2_z_score, z = AveExp_z_score))  + geom_point(aes(colour=AveExp_z_score)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,  ": z_score AveExpv vs. Conections : All points",sep=""))+ geom_contour()+ theme(legend.position="none")  # + theme(legend.position="none")
    m2<-ggplot(Interactomes_GC3_T2_selected, aes(Conections_z_score, T2_z_score, z = AveExp_z_score))  + geom_point(aes(colour=AveExp_z_score)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,  ": z_score AveExpv vs .Conections : 0-100",sep=""))+ geom_contour()  + xlim(0, 100) + ylim(0, 100)     # + theme(legend.position="none")
    m3<-ggplot(Interactomes_GC3_T2_selected, aes(Conections_z_score, T2_z_score, z = AveExp_z_score))  + geom_point(aes(colour=AveExp_z_score)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,  ": z_score AveExpv vs .Conections : 0-50",sep=""))+ geom_contour()  + xlim(0, 50) + ylim(0, 50)      #  + theme(legend.position="none")        
    m4<-ggplot(Interactomes_GC3_T2_selected, aes(Conections_z_score, T2_z_score, z = AveExp_z_score))  + geom_point(aes(colour=AveExp_z_score)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,  ": z_score AveExpv vs .Conections : 0-2.5 - 0-30",sep=""))+ geom_contour()  + xlim(0, 2.5) + ylim(0, 30)      #  + theme(legend.position="none")          
    #########################################################################################################################################    
    # FindClusters_resolution               
    png(filename=paste(output_dir,"countour_T2_Coonections_zscore_",normalization_scheme,".png",sep=""), width = 25, height = 25, res=600, units = "cm")  
            ggarrange(m1,m2,m3,m4,nrow = 2,ncol = 2, common.legend = TRUE, legend="bottom")
    dev.off()  
    #########################################################################################################################################    
    #########################################################################################################################################
    # Change collumn id
    colnames(melt_expression_interactomes)[6]<-"Expr"
  
    # Filter up Average expression greater than zero
    melt_expression_interactomes<-melt_expression_interactomes[melt_expression_interactomes$Expr>0,]
    melt_expression_interactomes<-melt_expression_interactomes[melt_expression_interactomes$Conections>0,]
    melt_expression_interactomes<-melt_expression_interactomes[melt_expression_interactomes$T2>0,]

    m1<-ggplot(melt_expression_interactomes, aes(Conections_z_score, T2_z_score, z = AveExp_z_score))  + geom_point(aes(colour=AveExp_z_score)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,  ": Expv vs .Conections : All points",sep=""))+ geom_contour()+ theme(legend.position="none")  # + theme(legend.position="none")
    m2<-ggplot(melt_expression_interactomes, aes(Conections_z_score, T2_z_score, z = AveExp_z_score))  + geom_point(aes(colour=AveExp_z_score)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,  ": Expv vs .Conections : 0-100",sep=""))+ geom_contour()     + xlim(0, 100) + ylim(0, 100)     # + theme(legend.position="none")
    m3<-ggplot(melt_expression_interactomes, aes(Conections_z_score, T2_z_score, z = AveExp_z_score))  + geom_point(aes(colour=AveExp_z_score)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,  ": Expv vs .Conections : 0-50",sep=""))+ geom_contour()      + xlim(0, 50) + ylim(0, 50)      #  + theme(legend.position="none")        
    m4<-ggplot(melt_expression_interactomes, aes(Conections_z_score, T2_z_score, z = AveExp_z_score))  + geom_point(aes(colour=AveExp_z_score)) + geom_density_2d_filled() + theme_bw() + ggtitle(paste(normalization_scheme,  ": Expv vs .Conections : 0-25",sep=""))+ geom_contour()      + xlim(0, 25) + ylim(0, 25)      #  + theme(legend.position="none")          
    #########################################################################################################################################    
    # FindClusters_resolution               
    png(filename=paste(output_dir,"countour_T2_Coonections_zscore_",normalization_scheme,".png",sep=""), width = 20, height = 20, res=600, units = "cm")  
            ggarrange(m1,m2,m3,m4,nrow = 2,ncol = 2, common.legend = TRUE, legend="bottom")
    dev.off()  
    #########################################################################################################################################      				   
}
