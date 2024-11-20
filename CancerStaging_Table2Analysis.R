library("xlsx")


# Set Table4.xlsx
Table4 <-data.frame(read_excel("/home/felipe/Downloads/Table4.xlsx", skip = 1))

# Set rownames
rownames(Table4)<-Table4$ENSEMBL

# Table 2 analysis
normalization_scheme<-"tpm"

genes_stages_I    <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_I",".tsv",sep=""), sep = '\t', header = TRUE)$gene #
genes_stages_II   <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_II",".tsv",sep=""), sep = '\t', header = TRUE)$gene #
genes_stages_III  <-read.table(file = paste(output_dir,"/FindStageSpecificGenes_",normalization_scheme,"_","sample_stage_III",".tsv",sep=""), sep = '\t', header = TRUE)$gene 


tumor_genes<-c("ENSG00000141338", "ENSG00000134917", "ENSG00000171094 ", "ENSG00000134982", "ENSG00000196914", "ENSG00000108381", "ENSG00000087586", "ENSG00000178999", "ENSG00000089685", "ENSG00000157764", "ENSG00000012048", "ENSG00000129993", "ENSG00000134057", "ENSG00000157456", "ENSG00000170312", "ENSG00000129757", "ENSG00000138180", "ENSG00000149554", "ENSG00000157404", "ENSG00000105976", "ENSG00000136848", "ENSG00000187323", "ENSG00000162733", "ENSG00000108654", "ENSG00000146648", "ENSG00000141736", "ENSG00000012061", "ENSG00000182197", "ENSG00000066468", "ENSG00000022267", "ENSG00000114861", "ENSG00000147257", "ENSG00000141736", "ENSG00000096968", "ENSG00000198553", "ENSG00000129451", "ENSG00000133703", "ENSG00000150457", "ENSG00000144791", "ENSG00000171444", "ENSG00000105976", "ENSG00000158747", "ENSG00000198400", "ENSG00000141510", "ENSG00000132646", "ENSG00000188389", "ENSG00000120217", "ENSG00000197646", "ENSG00000121879", "ENSG00000185920", "ENSG00000171862", "ENSG00000183010", "ENSG00000122679", "ENSG00000116473", "ENSG00000101265", "ENSG00000139687", "ENSG00000122707", "ENSG00000165731", "ENSG00000143878", "ENSG00000047936", "ENSG00000167325", "ENSG00000111961", "ENSG00000133121", "ENSG00000154144", "ENSG00000146648", "ENSG00000131747", "ENSG00000141510", "ENSG00000088325", "ENSG00000204977", "ENSG00000176890", "ENSG00000175063") 

unique_stage_I    =intersect(setdiff(genes_stages_I, c(genes_stages_II,genes_stages_III)),genes_stages_I)
unique_stage_II   =intersect(setdiff(genes_stages_II, c(genes_stages_I,genes_stages_III)),genes_stages_II)
unique_stage_III  =intersect(setdiff(genes_stages_III, c(genes_stages_I,genes_stages_II)),genes_stages_III)

rownames(df_reads_count_all_projects) %in% tumor_genes
