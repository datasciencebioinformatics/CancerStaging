# All tumor and control samples
colData_tumor  <-unique(merged_data_patient_info_count[merged_data_patient_info_count$Sample.Type=="Primary Tumor",])
colData_normal <-unique(merged_data_patient_info_count[merged_data_patient_info_count$Sample.Type=="Solid Tissue Normal",])
                                                                                                                                     #
# Vector with each stage                                                                                                              #
stages_str<-c("stage_I","stage_II","stage_III")                                                                                       #
                                                                                                                                      #
# Samples of each stage stored in colData                                                                                             #
sample_stage_I  <-unique(colData_tumor[colData_tumor$stages=="Stage I","sample_id"])                                                          #
sample_stage_II <-unique(colData_tumor[colData_tumor$stages=="Stage II","sample_id"])                                                         #
sample_stage_III<-unique(colData_tumor[colData_tumor$stages=="Stage III","sample_id"])                                                        #
sample_normal   <-unique(colData_normal[,"sample_id"])      
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
	#######################################################################################################################################
	# Select stages 
	stages_I_II_III_unique<-ggVennDiagram(list(Stage_I=unique_stage_I,Stage_II=unique_stage_II,Stage_III=unique_stage_III), label_alpha = 0.9,set_color = c("grey50","grey50","grey50")) +  scale_fill_gradient(low = "white", high = "white") + theme_bw() + ggtitle("Stages I, II and III")+ guides(fill="none")
	stages_I_II_III_unique<-ggVennDiagram(list(Stage_I=selected_genes_Stage_I_gene,Stage_II=selected_genes_Stage_II_gene,Stage_III=selected_genes_Stage_III_gene), label_alpha = 0.9,set_color = c("grey50","grey50","grey50")) +  scale_fill_gradient(low = "white", high = "white") + theme_bw() + ggtitle("Stages I, II and III")+ guides(fill="none") + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme_bw() +theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())  + ggtitle("B") 
		
	# FindClusters_resolution
	png(filename=paste(output_dir,"stages_I_II_III_unique.png",sep=""), width = 12.0, height = 12, res=1600, units = "cm")
		stages_I_II_III_unique
	dev.off()	
	
	#############################################################################################################################################################################  
	# List of used genes
	genes<-unique(c(selected_genes_Stage_I_gene,selected_genes_Stage_II_gene,selected_genes_Stage_III_gene))
	
	# List of considered samples
	samples<-c(sample_stage_I,sample_stage_II,sample_stage_III,sample_normal)
	
	# Store normalized table
	normalized_table      <-df_reads_count_all_projects[[normalized_table_names]][genes,samples]
	normalized_table_paired<-df_reads_count_all_projects[[normalized_table_names]][genes,c(paired_sample_df$normal,paired_sample_df$tumor)]
	
	# Merge merged_data_patient_sel
	merged_data_patient_sel<-unique(merged_data_patient_info[merged_data_patient_info$sample_id %in% colnames(normalized_table),c("Sample.Type","stages","sample_id")])  
	
	# Set rownames
	rownames(merged_data_patient_sel)<-merged_data_patient_sel$sample_id
	
	# Tanspose RPKM table                                                                                                                                                               #	
	transporse_normalized_table<-data.frame(t(normalized_table))                                                                                                                        #
	transporse_normalized_table_paired<-data.frame(t(normalized_table_paired))                                                                                                                        #
																		    #
	# Calculate prcomp for stage                                                                                                                                                        #
	pca_res_tumor_normal   <- prcomp(transporse_normalized_table, scale. = TRUE)                                                                                                              #
	pca_res_tumor_normal_paired   <- prcomp(transporse_normalized_table_paired, scale. = TRUE)                                                                                                              #
	
	# Rename collumns
	colnames(merged_data_patient_sel)[]<-"tumor_normal"
	
	# Plot PCA tumor versus normal                                                                                                                                                      #
	plot_res_tumor_normal        <- autoplot(pca_res_tumor_normal, data = merged_data_patient_sel[rownames(transporse_normalized_table),], colour = 'tumor_normal')+ theme_bw()  + theme(legend.position="bottom") + ggtitle("A")               + scale_color_manual(values=c("#E69F00", "#56B4E9")) + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16),, axis.title.x = element_text(size = 16)) + theme(legend.text=element_text(size=16))
	plot_res_tumor_normal_2      <- fviz_pca_ind(pca_res_tumor_normal, col.ind = merged_data_patient_sel[rownames(transporse_normalized_table),"tumor_normal"],palette = c("#E69F00", "#56B4E9"), col.var = c("#E69F00", "#56B4E9"),addEllipses = TRUE, ellipse.level = 0.95,labels="none",labelsize=0) + theme(legend.position="bottom") + ggtitle("A")               + scale_color_manual(values=c("#E69F00", "#56B4E9")) + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16),, axis.title.x = element_text(size = 16)) + theme(legend.text=element_text(size=16))

	
	# FindClusters_resolution
	png(filename=paste(output_dir,"plot_res_tumor_normal.png",sep=""), width = 12.0, height = 12, res=1600, units = "cm")
		plot_res_tumor_normal
	dev.off()

	# FindClusters_resolution
	png(filename=paste(output_dir,"plot_res_tumor_normal_2.png",sep=""), width = 12.0, height = 12, res=1600, units = "cm")
		plot_res_tumor_normal_2
	dev.off()
	
	
	
	#############################################################################################################################
	write_tsv(selected_genes_Stage_I_data[unique_stage_I,],     paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","unique_stage_I",".tsv",sep=""))
	write_tsv(selected_genes_Stage_II_data[unique_stage_II,],   paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","unique_stage_II",".tsv",sep=""))
	write_tsv(selected_genes_Stage_III_data[unique_stage_III,], paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","unique_stage_III",".tsv",sep=""))
	#############################################################################################################################
	cat(print(paste("\n", normalization_scheme, " Number of tumor genes per stage for Stage I : ",paste(length(rownames(selected_genes_Stage_I_data)),"/",length(unique_stage_I),sep=""))),file=paste(output_dir,"outfile.txt",sep="/"),append=TRUE)
	cat(print(paste("\n", normalization_scheme, " Number of tumor genes per stage for Stage II : ",paste(length(rownames(selected_genes_Stage_II_data)),"/",length(unique_stage_II),sep=""))),file=paste(output_dir,"outfile.txt",sep="/"),append=TRUE)
	cat(print(paste("\n", normalization_scheme, " Number of tumor genes per stage for Stage III : ",paste(length(rownames(selected_genes_Stage_III_data)),"/",length(unique_stage_III),sep=""))),file=paste(output_dir,"outfile.txt",sep="/"),append=TRUE)
}
