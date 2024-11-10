normalization_scheme<-"tpm"
#######################################################################################################################################
# Path to files of selected_genes                                                                                                             # 
# genes_stages_I
selected_genes_Stage_I_data    <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_I",".tsv",sep=""), sep = '\t', header = TRUE) #
selected_genes_Stage_II_data   <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_II",".tsv",sep=""), sep = '\t', header = TRUE) #
selected_genes_Stage_III_data  <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_III",".tsv",sep=""), sep = '\t', header = TRUE) #

rownames(selected_genes_Stage_I_data)<-selected_genes_Stage_I_data$gene
rownames(selected_genes_Stage_II_data)<-selected_genes_Stage_II_data$gene
rownames(selected_genes_Stage_III_data)<-selected_genes_Stage_III_data$gene

selected_genes_Stage_I_gene      <- selected_genes_Stage_I_data$gene
selected_genes_Stage_II_gene     <- selected_genes_Stage_II_data$gene
selected_genes_Stage_III_gene    <- selected_genes_Stage_III_data$gene  #######################################################################################################################################                                                                                                                                     #
  unique_stage_I  =intersect(setdiff(selected_genes_Stage_I_gene, c(selected_genes_Stage_II_gene,selected_genes_Stage_III_gene)),selected_genes_Stage_I_gene)
  unique_stage_II =intersect(setdiff(selected_genes_Stage_II_gene, c(selected_genes_Stage_I_gene,selected_genes_Stage_III_gene)),selected_genes_Stage_II_gene)
  unique_stage_III=intersect(setdiff(selected_genes_Stage_III_gene, c(selected_genes_Stage_I_gene,selected_genes_Stage_II_gene)),selected_genes_Stage_III_gene)

# Considering that the average connection score of up-regulated genes across the 1722 of the three stages combined was 50.06 by reference to the IntAct reactome,
intersect_conectivity<-connectivity[names(connectivity) %in% unique(c(selected_genes_Stage_I_gene,selected_genes_Stage_II_gene,selected_genes_Stage_III_gene))]

mean(intersect_conectivity)

# we chose 50 as a threshold above which to consider a protein with a higher connection score as a hub. 
intersect_conectivity<-intersect_conectivity[which(intersect_conectivity>75)]

# Because therapy is intended to maximize the patient comfort we also filtered out genes whose expression was larger than ~10 in the control since targeting drugs could affect the healthy tissue in case of basal expression.
selected_genes<-df_FC[df_FC$Mean_normal<=10,]


# Store conectivity
df_selected_conectivity<-data.frame(intersect_conectivity[names(intersect_conectivity) %in% selected_genes$Gene])

# rownames(df_selected_conectivity)
rownames(df_selected_conectivity)<-df_selected_conectivity$Var1

# Add stage
df_selected_conectivity_Stage_I<-na.omit(df_selected_conectivity[unique_stage_I,])
df_selected_conectivity_Stage_II<-na.omit(df_selected_conectivity[unique_stage_II,])
df_selected_conectivity_Stage_III<-na.omit(df_selected_conectivity[unique_stage_III,])

# Add stage
df_selected_conectivity_Stage_I$Stage<-"Stage I"
df_selected_conectivity_Stage_II$Stage<-"Stage II"
df_selected_conectivity_Stage_III$Stage<-"Stage III"

# Merge tables
df_selected_conectivity<-rbind(df_selected_conectivity_Stage_I,df_selected_conectivity_Stage_II,df_selected_conectivity_Stage_III)

# Set colnames
colnames(df_selected_conectivity)<-c("Gene","Connectivity","Stage")

# Mer data.frame
df_selected_conectivity<-merge(df_selected_conectivity,df_FC,by="Gene")[,c("Gene", "Connectivity", "Stage", "Mean_normal", "sd_normal", "Mean_tumor","sd_tumor","FC")]

# Save TSV file with genes from Stage3
write_tsv(df_selected_conectivity, paste(output_dir,"/selected_genes.tsv",sep=""))
