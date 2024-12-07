normalization_scheme<-"tpm"
#######################################################################################################################################
rownames(genes_rankData_stage_all_genes_merged)<-genes_rankData_stage_all_genes_merged$ensembl_gene_id
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
selected_genes_Stage_III_gene    <- selected_genes_Stage_III_data$gene  
#######################################################################################################################################                                                                                                                                     #
unique_stage_I  =intersect(setdiff(selected_genes_Stage_I_gene, c(selected_genes_Stage_II_gene,selected_genes_Stage_III_gene)),selected_genes_Stage_I_gene)
unique_stage_II =intersect(setdiff(selected_genes_Stage_II_gene, c(selected_genes_Stage_I_gene,selected_genes_Stage_III_gene)),selected_genes_Stage_II_gene)
unique_stage_III=intersect(setdiff(selected_genes_Stage_III_gene, c(selected_genes_Stage_I_gene,selected_genes_Stage_II_gene)),selected_genes_Stage_III_gene)

# Considering that the average connection score of up-regulated genes across the 1722 of the three stages combined was 50.06 by reference to the IntAct reactome,
intersect_conectivity<-genes_rankData_stage_all_genes_merged[unique(c(selected_genes_Stage_I_gene,selected_genes_Stage_II_gene,selected_genes_Stage_III_gene)),]

# we chose 50 as a threshold above which to consider a protein with a higher connection score as a hub. 
intersect_conectivity<-na.omit(intersect_conectivity[intersect_conectivity$Conections>75,])

# Because therapy is intended to maximize the patient comfort we also filtered out genes whose expression was larger than ~10 in the control since targeting drugs could affect the healthy tissue in case of basal expression.
selected_genes<-intersect_conectivity[intersect_conectivity$Mean_normal<=10,]

selected_genes_Stage_I<-na.omit(selected_genes[unique_stage_I,])
selected_genes_Stage_II<-na.omit(selected_genes[unique_stage_II,])
selected_genes_Stage_III<-na.omit(selected_genes[unique_stage_III,])

selected_genes_Stage_I$FC<-selected_genes_Stage_I$FC_Stage_I
selected_genes_Stage_II$FC<-selected_genes_Stage_II$FC_Stage_II
selected_genes_Stage_III$FC<-selected_genes_Stage_II$FC_Stage_III

selected_genes_Stage_I$Mean_stage<-selected_genes_Stage_I$Mean_Stage_I
selected_genes_Stage_II$Mean_stage<-selected_genes_Stage_II$Mean_Stage_II
selected_genes_Stage_III$Mean_stage<-selected_genes_Stage_III$Mean_Stage_III

selected_genes_Stage_I$SD_stage<-selected_genes_Stage_I$sd_Stage_I
selected_genes_Stage_II$SD_stage<-selected_genes_Stage_II$sd_Stage_II
selected_genes_Stage_III$SD_stage<-selected_genes_Stage_III$sd_Stage_III

selected_genes_Stage_I$Cnx<-selected_genes_Stage_I$Conections
selected_genes_Stage_II$Cnx<-selected_genes_Stage_II$Conections
selected_genes_Stage_III$Cnx<-selected_genes_Stage_III$Conections


# Selected all sategs
selected_genes_ALL_STAGES<-rbind(selected_genes_Stage_I,selected_genes_Stage_II,selected_genes_Stage_III)

# Select all tables
selected_genes_ALL_STAGES<-selected_genes_ALL_STAGES[,c("hgnc_symbol","Mean_normal","sd_normal","Mean_stage","SD_stage","Cnx")]

# Save TSV file with genes from Stage3
write_tsv(selected_genes_ALL_STAGES, paste(output_dir,"/Table4.tsv",sep=""))


#######################################################################################################################################
stage_I_II    =unique(intersect(selected_genes_Stage_I_gene,selected_genes_Stage_II_gene))
stage_I_III   =unique(intersect(selected_genes_Stage_I_gene,selected_genes_Stage_III_gene))
stage_II_III  =unique(intersect(selected_genes_Stage_II_gene,selected_genes_Stage_III_gene))


# Considering that the average connection score of up-regulated genes across the 1722 of the three stages combined was 50.06 by reference to the IntAct reactome,
intersect_conectivity_stage_I_II   <- stage_I_II[which(stage_I_II %in% names(connectivity))]
intersect_conectivity_stage_I_III  <- stage_I_III[which(stage_I_III %in% names(connectivity))]
intersect_conectivity_stage_II_III <- stage_II_III[which(stage_II_III %in% names(connectivity))]

# we chose 50 as a threshold above which to consider a protein with a higher connection score as a hub. 
intersect_conectivity_stage_I_II <-intersect_conectivity_stage_I_II[which(intersect_conectivity_stage_I_II>55)]
intersect_conectivity_stage_I_III<-intersect_conectivity_stage_I_III[which(intersect_conectivity_stage_I_III>55)]
intersect_conectivity_stage_II_III<-intersect_conectivity_stage_II_III[which(intersect_conectivity_stage_II_III>55)]

# Because therapy is intended to maximize the patient comfort we also filtered out genes whose expression was larger than ~10 in the control since targeting drugs could affect the healthy tissue in case of basal expression.
selected_genes<-genes_rankData_stage_all_genes[genes_rankData_stage_all_genes$Mean_normal<=10.0,]

# Store conectivity
df_selected_conectivity_stage_I_II<-data.frame(Gene=intersect_conectivity_stage_I_II[intersect_conectivity_stage_I_II %in% selected_genes$ensembl_gene_id])
df_selected_conectivity_stage_I_III<-data.frame(Gene=intersect_conectivity_stage_I_III[intersect_conectivity_stage_I_III %in% selected_genes$ensembl_gene_id])
df_selected_conectivity_stage_II_III<-data.frame(Gene=intersect_conectivity_stage_II_III[intersect_conectivity_stage_II_III %in% selected_genes$ensembl_gene_id])

# Add stage
df_selected_conectivity_stage_I_II$Stage  <-"Stages_I_II"
df_selected_conectivity_stage_I_III$Stage <-"Stages_I_III"
df_selected_conectivity_stage_II_III$Stage<-"Stages_II_III"

colnames(df_selected_conectivity_stage_I_II)<-c("Gene","Stage")
colnames(df_selected_conectivity_stage_I_III)<-c("Gene","Stage")
colnames(df_selected_conectivity_stage_II_III)<-c("Gene","Stage")

# Merge data.frame
df_selected_conectivity_stage_all<-unique(rbind(df_selected_conectivity_stage_I_II,df_selected_conectivity_stage_I_III,df_selected_conectivity_stage_II_III))

# Mer data.frame
# Re set colnames
colnames(genes_rankData_stage_all_genes)[1]<-"Gene"

# Merge tables
df_selected_conectivity<-merge(df_selected_conectivity_stage_all,genes_rankData_stage_all_genes,by="Gene")


# Save TSV file with genes from Stage3
write_tsv(df_selected_conectivity, paste(output_dir,"/selected_conectivity_2.tsv",sep=""))












#######################################################################################################################################                                                                                                                                     #
unique_stage=unique(c(selected_genes_Stage_I_gene,selected_genes_Stage_II_gene,selected_genes_Stage_III_gene))

# Considering that the average connection score of up-regulated genes across the 1722 of the three stages combined was 50.06 by reference to the IntAct reactome,
intersect_conectivity<-connectivity[names(connectivity) %in% unique_stage]

# we chose 50 as a threshold above which to consider a protein with a higher connection score as a hub. 
intersect_conectivity<-intersect_conectivity[which(intersect_conectivity>50)]

# Take DE genes with normal mean <= 10.0
selected_genes<-genes_rankData_stage_all_genes[rownames(genes_rankData_stage_all_genes) %in% unique_stage,]

# Because therapy is intended to maximize the patient comfort we also filtered out genes whose expression was larger than ~10 in the control since targeting drugs could affect the healthy tissue in case of basal expression.
selected_genes<-genes_rankData_stage_all_genes[which(selected_genes$Mean_normal<=10.0),]

# Filter also by FC
selected_genes<-selected_genes[which(selected_genes$FC>=2.0),]

# Among the 1808 up-regulated tumor genes, 167 genes have mean expression in normal samples TPM<=10 and FC>=2. Amnong these, 19 biomarkers have also conectivity greater thabn 50.

# Store conectivity
df_selected_conectivity<-data.frame(intersect_conectivity[names(intersect_conectivity) %in% rownames(selected_genes)])

# rownames(df_selected_conectivity)
rownames(df_selected_conectivity)<-df_selected_conectivity$Var1

# Set colnames
colnames(df_selected_conectivity)<-c("ensembl_gene_id","Connectivity")

# Mer data.frame
df_selected_conectivity<-merge(df_selected_conectivity,df_FC,by="ensembl_gene_id")

# Slecte collumns
df_selected_conectivity<-df_selected_conectivity[,c("ensembl_gene_id","Connectivity","Mean_normal","sd_normal","Mean_tumor","sd_tumor","FC")]

# Set colnames
colnames(df_selected_conectivity)<-c("Gene", "Cnx", "Avg.normal", "sd.normal", "Avg.tumor","sd.tumor","FC")

# Re-order table
df_selected_conectivity<-df_selected_conectivity[,c("Gene","Avg.normal", "sd.normal","Avg.tumor","sd.tumor","FC", "Cnx")]

# Change to gene symbols
df_selected_conectivity$Gene<-genes_df_FC[genes_df_FC$ensembl_gene_id %in% df_selected_conectivity$Gene,"hgnc_symbol"]

  

# Save TSV file with genes from Stage3
write_tsv(df_selected_conectivity, paste(output_dir,"/Table_3.tsv",sep=""))
