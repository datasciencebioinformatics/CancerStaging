#######################################################################################################################################
# A script to caluclate entropy from lists of genes from each stage
#######################################################################################################################################
# Interactome data
# for each  normalization scheme
for (normalization_scheme in normalization_schemes)
{	
	# Take the expression data
	# First, I will load the expression table   	
	normalized_expression_table<-df_reads_count_all_projects[[normalization_scheme]]

	# genes_ids
	genes_ids <- rownames(normalized_expression_table)	
		
	# Store log2change_Stage_i
	log2change_Stage_i   <-list_stage_specific_genes[[paste(normalization_scheme,"_","stage_I",sep="")]]
	log2change_Stage_ii  <-list_stage_specific_genes[[paste(normalization_scheme,"_","stage_II",sep="")]]
	log2change_Stage_iii <-list_stage_specific_genes[[paste(normalization_scheme,"_","stage_III",sep="")]]	
	
	# Vectors to store gene ids from each stage
	genes_id_vector_stage_I<-log2change_Stage_i[log2change_Stage_i$Category=="Per stage genes","gene"]
	genes_id_vector_stage_II<-log2change_Stage_ii[log2change_Stage_ii$Category=="Per stage genes","gene"]
	genes_id_vector_stage_III<-log2change_Stage_iii[log2change_Stage_iii$Category=="Per stage genes","gene"]
	
	# Take all genes from interactom
	# store both, gene in pair one and gene in pair two in a same vectors
	
	# The gene in lists of genes per stage must be filterd to keep only entris that are present in the interactome
	# ~99% of genes in the selected lists are in the interactome
	genes_interactome_stage_I  <-genes_id_vector_stage_I[genes_id_vector_stage_I %in% genes_ids]
	genes_interactome_stage_II <-genes_id_vector_stage_II[genes_id_vector_stage_II %in% genes_ids]
	genes_interactome_stage_III<-genes_id_vector_stage_III[genes_id_vector_stage_III %in% genes_ids]

	print(paste(length(genes_interactome_stage_I),length(genes_interactome_stage_II),length(genes_interactome_stage_III), sep="-"))
	
	# A vector with all genes of the interactome,full_interactome<-unique(c(genes_interactome_stage_I$Gene1,genes_interactome_stage_I$Gene2))
	# Calculate all pairwise combinations of genes, without redunctancy
	full_interactome_stage_I<- data.frame(expand.grid.unique(x = genes_interactome_stage_I, y = genes_interactome_stage_I,include.equals=FALSE))
	full_interactome_stage_II<- data.frame(expand.grid.unique(x = genes_interactome_stage_II, y = genes_interactome_stage_II,include.equals=FALSE))
	full_interactome_stage_III<- data.frame(expand.grid.unique(x = genes_interactome_stage_III, y = genes_interactome_stage_III,include.equals=FALSE))

	# If at least one overlapping interaction
	if(dim(full_interactome_stage_I)[1]>0)
	{		
		# set colnames
		colnames(full_interactome_stage_I)<-c("Gene1","Gene2")
	}
	# If at least one overlapping interaction
	if(dim(full_interactome_stage_II)[1]>0)
	{		
		# set colnames
		colnames(full_interactome_stage_II)<-c("Gene1","Gene2")
	}
	# If at least one overlapping interaction
	if(dim(full_interactome_stage_III)[1]>0)
	{		
		# set colnames
		colnames(full_interactome_stage_III)<-c("Gene1","Gene2")
	}	  
	  	  	
	rownames(full_interactome_stage_I)  <- paste(full_interactome_stage_I$Gene1,full_interactome_stage_I$Gene2,sep="-")
	rownames(full_interactome_stage_II) <-paste(full_interactome_stage_II$Gene1,full_interactome_stage_II$Gene2,sep="-")
	rownames(full_interactome_stage_III)<-paste(full_interactome_stage_III$Gene1,full_interactome_stage_III$Gene2,sep="-")
	#######################################################################################################
	interactome_stage_I  <-interactome_data[which(rownames(interactome_data) %in% rownames(full_interactome_stage_I)),]
	interactome_stage_II <-interactome_data[which(rownames(interactome_data) %in% rownames(full_interactome_stage_II)),]
	interactome_stage_III<-interactome_data[which(rownames(interactome_data) %in% rownames(full_interactome_stage_III)),]
	########################################################################################################################################
	stage_I_genes_factor  <-factor(c(interactome_stage_I$Gene1,interactome_stage_I$Gene2),level=unique(c(interactome_stage_I$Gene1,interactome_stage_I$Gene2)))
	stage_II_genes_factor <-factor(c(interactome_stage_II$Gene1,interactome_stage_II$Gene2),level=unique(c(interactome_stage_II$Gene1,interactome_stage_II$Gene2)))
	stage_III_genes_factor<-factor(c(interactome_stage_III$Gene1,interactome_stage_III$Gene2),level=unique(c(interactome_stage_III$Gene1,interactome_stage_III$Gene2)))

	print(paste(length(unique(stage_I_genes_factor)),length(unique(stage_II_genes_factor)),length(unique(stage_III_genes_factor)), sep="-"))

	   
	df_stageI_connectivity   <-unique(data.frame(Conectivity=table(stage_I_genes_factor)))
	df_stageII_connectivity  <-unique(data.frame(Conectivity=table(stage_II_genes_factor)))
	df_stageIII_connectivity <-unique(data.frame(Conectivity=table(stage_III_genes_factor)))	
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
	print(paste(normalization_scheme," : Stage I :", round(Entropy_stage_I_value_Carels,3), " : Stage II :", round(Entropy_stage_II_value_Carels,3), " : Stage III :",  round(Entropy_stage_III_value_Carels,3), sep="") )
	cat(print(paste("\n",normalization_scheme," : Stage I :",  round(Entropy_stage_I_value_Carels,3), " : Stage II :",  round(Entropy_stage_II_value_Carels,3), " : Stage III :",  round(Entropy_stage_III_value_Carels,3),"\n", sep="") ),file=results_files,append=TRUE)
	
}
