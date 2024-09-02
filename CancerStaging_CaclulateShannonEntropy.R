#######################################################################################################################################
# A script to caluclate entropy from lists of genes from each stage
#######################################################################################################################################
df_reads_count_all_projects[[1]]

# A filter to keep only genes that are positivelly regulated
genes_ids<-c()

# For each gene in stage I
for (gene_id in rownames(unstranded_data_filter))
{
  # Store gene id in the vector
  genes_ids<-c(genes_ids,strsplit(gene_id, split = "\\.")[[1]][1])
}
# Gene_ids
genes_ids<-unique(genes_ids)
########################################################################################################################################
# Interactome data
interactome_data

# for each  normalization scheme
for (normalization_scheme in normalization_schemes)
{	
  # Take the expression data

  # A filter to keep only genes that are positivelly regulated
  
  
  # for each  normalization scheme
  for (comparisson_index in rownames(df_table_comparisson))
  {
    # Store log2change_Stage_i
    log2change_Stage_i<-list_stage_specific_genes[[paste(normalization_scheme,"_",substring(Stage_i,8,20),sep="")]]

    # Store stage-specif genes
    stage_specifi_genes<-log2change_Stage_i[log2change_Stage_i$Category=="Per stage genes","gene"]
  }  
}


# Gene table
genes_Stage_I       <-read.table(file = file_genes_Stage_I, sep = '\t', header = TRUE,fill=TRUE)         
genes_Stage_II      <-read.table(file = file_genes_Stage_II, sep = '\t', header = TRUE,fill=TRUE)
genes_Stage_III     <-read.table(file = file_genes_Stage_III, sep = '\t', header = TRUE,fill=TRUE)   
########################################################################################################################################

#######################################################################################################      
# Take all genes from interactom
# store both, gene in pair one and gene in pair two in a same vectors

# The gene in lists of genes per stage must be filterd to keep only entris that are present in the interactome
# ~99% of genes in the selected lists are in the interactome
genes_interactome_stage_I  <-genes_id_vector_stage_I[genes_id_vector_stage_I %in% genes_ids]
genes_interactome_stage_II <-genes_id_vector_stage_II[genes_id_vector_stage_II %in% genes_ids]
genes_interactome_stage_III<-genes_id_vector_stage_III[genes_id_vector_stage_III %in% genes_ids]

# A vector with all genes of the interactome,full_interactome<-unique(c(genes_interactome_stage_I$Gene1,genes_interactome_stage_I$Gene2))
# Calculate all pairwise combinations of genes, without redunctancy
full_interactome_stage_I<- data.frame(expand.grid.unique(x = genes_interactome_stage_I, y = genes_interactome_stage_I,include.equals=FALSE))
full_interactome_stage_II<- data.frame(expand.grid.unique(x = genes_interactome_stage_II, y = genes_interactome_stage_II,include.equals=FALSE))
full_interactome_stage_III<- data.frame(expand.grid.unique(x = genes_interactome_stage_III, y = genes_interactome_stage_III,include.equals=FALSE))

# set colnames
colnames(full_interactome_stage_I)<-c("Gene1","Gene2")
colnames(full_interactome_stage_II)<-c("Gene1","Gene2")
colnames(full_interactome_stage_III)<-c("Gene1","Gene2")

rownames(full_interactome_stage_I)<-paste(full_interactome_stage_I$Gene1,full_interactome_stage_I$Gene2,sep="-")
rownames(full_interactome_stage_II)<-paste(full_interactome_stage_II$Gene1,full_interactome_stage_II$Gene2,sep="-")
rownames(full_interactome_stage_III)<-paste(full_interactome_stage_III$Gene1,full_interactome_stage_III$Gene2,sep="-")
#######################################################################################################
interactome_stage_I  <-interactome_data_inv[which(rownames(interactome_data_inv) %in% rownames(full_interactome_stage_I)),]
interactome_stage_II <-interactome_data_inv[which(rownames(interactome_data_inv) %in% rownames(full_interactome_stage_II)),]
interactome_stage_III<-interactome_data_inv[which(rownames(interactome_data_inv) %in% rownames(full_interactome_stage_III)),]
########################################################################################################################################
stage_I_genes_factor  <-factor(c(interactome_stage_I$Gene1,interactome_stage_I$Gene2),level=unique(c(interactome_stage_I$Gene1,interactome_stage_I$Gene2)))
stage_II_genes_factor <-factor(c(interactome_stage_II$Gene1,interactome_stage_II$Gene2),level=unique(c(interactome_stage_II$Gene1,interactome_stage_II$Gene2)))
stage_III_genes_factor<-factor(c(interactome_stage_III$Gene1,interactome_stage_III$Gene2),level=unique(c(interactome_stage_III$Gene1,interactome_stage_III$Gene2)))
   
df_stageI_connectivity   <-unique(data.frame(Conectivity=table(stage_I_genes_factor)))
df_stageII_connectivity  <-unique(data.frame(Conectivity=table(stage_II_genes_factor)))
df_stageIII_connectivity <-unique(data.frame(Conectivity=table(stage_III_genes_factor)))

#df_stageI_connectivity   <-unique(data.frame(Conectivity=table(c(interactome_stage_I$Gene1,interactome_stage_I$Gene2))))
#df_stageII_connectivity  <-unique(data.frame(Conectivity=table(c(interactome_stage_II$Gene1,interactome_stage_II$Gene2))))
#df_stageIII_connectivity <-unique(data.frame(Conectivity=table(c(interactome_stage_III$Gene1,interactome_stage_III$Gene2))))
########################################################################################################################################
colnames(df_stageI_connectivity)<-c("Gene","Conectivity")
colnames(df_stageII_connectivity)<-c("Gene","Conectivity")
colnames(df_stageIII_connectivity)<-c("Gene","Conectivity")
########################################################################################################################################
# Table for the calculation of entropy
df_entropy_calulation_I   <-data.frame(table(df_stageI_connectivity$Conectivity),p_k=0,log2_pk=0,p_k_mult_log2_pk=0)
df_entropy_calulation_II  <-data.frame(table(df_stageII_connectivity$Conectivity),p_k=0,log2_pk=0,p_k_mult_log2_pk=0)
df_entropy_calulation_III <-data.frame(table(df_stageIII_connectivity$Conectivity),p_k=0,log2_pk=0,p_k_mult_log2_pk=0)

# Rename colnames
colnames(df_entropy_calulation_I)<-c("k","count","p_k","log2_pk","p_k_mult_log2_pk")
colnames(df_entropy_calulation_II)<-c("k","count","p_k","log2_pk","p_k_mult_log2_pk")
colnames(df_entropy_calulation_III)<-c("k","count","p_k","log2_pk","p_k_mult_log2_pk")

# Calculate p(k)
df_entropy_calulation_I$p_k<-df_entropy_calulation_I$count/sum(df_entropy_calulation_I$count)
df_entropy_calulation_II$p_k<-df_entropy_calulation_II$count/sum(df_entropy_calulation_II$count)
df_entropy_calulation_III$p_k<-df_entropy_calulation_III$count/sum(df_entropy_calulation_III$count)

# Calculate log2(p(k))
df_entropy_calulation_I$log2_pk<-log(df_entropy_calulation_I$p_k,2)
df_entropy_calulation_II$log2_pk<-log(df_entropy_calulation_II$p_k,2)
df_entropy_calulation_III$log2_pk<-log(df_entropy_calulation_III$p_k,2)

# Calculate p(k)*log2(p(k))
df_entropy_calulation_I$p_k_mult_log2_pk<-df_entropy_calulation_I$p_k*df_entropy_calulation_I$log2_pk
df_entropy_calulation_II$p_k_mult_log2_pk<-df_entropy_calulation_II$p_k*df_entropy_calulation_II$log2_pk
df_entropy_calulation_III$p_k_mult_log2_pk<-df_entropy_calulation_III$p_k*df_entropy_calulation_III$log2_pk

# Caclulate entropy value
Entropy_stage_I_value_Carels  <-abs(sum(df_entropy_calulation_I$p_k_mult_log2_pk))
Entropy_stage_II_value_Carels <-abs(sum(df_entropy_calulation_II$p_k_mult_log2_pk))
Entropy_stage_III_value_Carels<-abs(sum(df_entropy_calulation_III$p_k_mult_log2_pk))
########################################################################################################################################
# Save TSV file with genes from Stage1
write_tsv(df_stageI_connectivity, paste(output_dir,"df_stageI_connectivity_I_interactome",".tsv",sep=""))
write_tsv(df_stageII_connectivity, paste(output_dir,"df_stageII_connectivity_II_interactome",".tsv",sep=""))
write_tsv(df_stageIII_connectivity, paste(output_dir,"df_stageIII_connectivity_III_interactome",".tsv",sep=""))

# Save TSV file with genes from Stage1
write_tsv(interactome_stage_I, paste(output_dir,"df_stageI_interactome_interactome",".tsv",sep=""))
write_tsv(interactome_stage_II, paste(output_dir,"df_stageII_interactome_interactome",".tsv",sep=""))
write_tsv(interactome_stage_III, paste(output_dir,"df_stageIII_interactome_interactome",".tsv",sep=""))
########################################################################################################################################
cat(print(paste("\nNº of vertex/Nº/Entropy of edges, sub-interactome network for Stage I: ",  paste(length(df_stageI_connectivity$Gene),dim(unique(interactome_stage_I))[1],round(Entropy_stage_I_value_Carels,4),sep="/"),sep="")),file=paste(output_dir,"outfile.txt",sep="/"),append=TRUE)
cat(print(paste("\nNº of vertex/Nº/Entropy of edges, sub-interactome network for Stage II: ", paste(length(df_stageII_connectivity$Gene),dim(unique(interactome_stage_II))[1],round(Entropy_stage_II_value_Carels,4),sep="/"),sep="")),file=paste(output_dir,"outfile.txt",sep="/"),append=TRUE)
cat(print(paste("\nNº of vertex/Nº/Entropy of edges, sub-interactome network for Stage III: ",paste(length(df_stageIII_connectivity$Gene),dim(unique(interactome_stage_III))[1],round(Entropy_stage_III_value_Carels,4),sep="/"),sep="")),file=paste(output_dir,"outfile.txt",sep="/"),append=TRUE)
