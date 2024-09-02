#######################################################################################################################################
# A script to caluclate entropy from lists of genes from each stage
#######################################################################################################################################
# Interactome data
# for each  normalization scheme
for (normalization_scheme in normalization_schemes)
{	
  # Take the expression data
  # First, I will load the expression table   	
  normalized_expression_table<-normalized_expression_table_list[[normalization_scheme]]

  # A filter to keep only genes that are positivelly regulated
  genes_ids<-c()
  
  # For each gene
  for (gene_id in rownames(normalized_expression_table))
  {
    # Store gene id in the vector
    genes_ids<-c(genes_ids,strsplit(gene_id, split = "\\.")[[1]][1])
  }
  # Gene_ids
  genes_ids<-unique(genes_ids)

  # for each  normalization scheme
  for (comparisson_index in rownames(df_table_comparisson))
  {
	# Store log2change_Stage_i
	log2change_Stage_i<-list_stage_specific_genes[[paste(normalization_scheme,"_",substring(Stage_i,8,20),sep="")]]
	
	# Store stage-specif genes
	DE_genes<-log2change_Stage_i[log2change_Stage_i$Category=="Per stage genes","gene"]
	
	# Stages
	Stage_i          <-df_table_comparisson[comparisson_index,"Stage_i"]
		
	# The gene in lists of genes per stage must be filterd to keep only entris that are present in the interactome
	# ~99% of genes in the selected lists are in the interactome
	genes_interactome_stage  <-DE_genes[DE_genes %in% genes_ids]
	
	# A vector with all genes of the interactome,full_interactome<-unique(c(genes_interactome_stage_I$Gene1,genes_interactome_stage_I$Gene2))
	# Calculate all pairwise combinations of genes, without redunctancy
	full_interactome_stage<- data.frame(expand.grid.unique(x = genes_interactome_stage, y = genes_interactome_stage,include.equals=FALSE))    

	# set colnames
	colnames(full_interactome_stage)<-c("Gene1","Gene2")    
	rownames(full_interactome_stage)<-paste(full_interactome_stage$Gene1,full_interactome_stage$Gene2,sep="-")

	# Check overlap between downloaded interactome 		  
	interactome_stage  <-interactome_data[which(rownames(interactome_data) %in% rownames(full_interactome_stage)),]
	########################################################################################################################################
	stage_genes_factor  <-factor(c(interactome_stage$Gene1,interactome_stage$Gene2),level=unique(c(interactome_stage$Gene1,interactome_stage$Gene2)))      
	df_stage_connectivity   <-unique(data.frame(Conectivity=table(stage_genes_factor)))    
	########################################################################################################################################
	colnames(df_stage_connectivity)<-c("Gene","Conectivity")    
	########################################################################################################################################
	# Table for the calculation of entropy
	df_entropy_calulation   <-data.frame(table(df_stage_connectivity$Conectivity),p_k=0,log2_pk=0,p_k_mult_log2_pk=0)
	
	# Rename colnames
	colnames(df_entropy_calulation)<-c("k","count","p_k","log2_pk","p_k_mult_log2_pk")
	
	# Calculate p(k)
	df_entropy_calulation$p_k<-df_entropy_calulation$count/sum(df_entropy_calulation$count)
	
	# Calculate log2(p(k))
	df_entropy_calulation$log2_pk<-log(df_entropy_calulation$p_k,2)
	
	# Calculate p(k)*log2(p(k))
	df_entropy_calulation$p_k_mult_log2_pk<-df_entropy_calulation$p_k*df_entropy_calulation$log2_pk
	
	# Caclulate entropy value
	Entropy_stage_value  <-abs(sum(df_entropy_calulation$p_k_mult_log2_pk))
	
	paste(dim)
	print(normalization_scheme," : ", Entropy_stage_value )
}
