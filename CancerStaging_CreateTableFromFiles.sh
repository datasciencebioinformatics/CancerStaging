# This bash script scan all downloaded .tsv files and create a table containing information for all samples.
# Output folder  /home/felipe/googledrive/Cancer_staging/tables/:
# unstranded.rna_seq.augmented_star_gene_counts.tsv          
# stranded_first.rna_seq.augmented_star_gene_counts.tsv
# stranded_second.rna_seq.augmented_star_gene_counts.tsv
# tpm_unstranded.rna_seq.augmented_star_gene_counts.tsv
# fpkm_unstranded.rna_seq.augmented_star_gene_counts.tsv
# header.txt ENSG gene id
# gene_name.txt gene symbol

# A script to obtain all expression data into a data table
input_folder=/home/felipe/googledrive/Cancer_staging/tables/

# gene_id	            #1
# gene_name	            #2  
# gene_typec	        #3
# unstranded	        #4	
# stranded_first	    #5
# stranded_second	    #6	
# tpm_unstranded	    #7	
#fpkm_unstranded	    #6
# List of files fo gene id
# For each file
ls  /home/felipe/googledrive/Cancer_staging/tables/*/*.rna_seq.augmented_star_gene_counts.tsv  | while read file
do
    folder_name=$(echo $file  | sed 's/\// /g' | awk '{print "/"$1"/"$2"/"$3"/"$4"/"$5"/"$6"/"}')
    file_name_2=$(echo $file  | sed 's/\// /g' |  sed 's/\./ /g' | awk '{print $7}')    
    file_name=$(echo $file  | sed 's/\// /g' | awk '{print $7}')    
    file_path=$(echo $folder_name"/"$file_name)    
    echo $file_name
    echo $folder_name$file_name_2".gene_id.txt"
    cat $file_path | grep -v GENCODE | grep -v gene_id | grep -v N_unmapped | grep -v N_multimapping | grep -v N_noFeature | grep -v N_ambiguous | awk '{print $1}' > $folder_name$file_name_2".gene_id.txt"
    cat $file_path | grep -v GENCODE | grep -v gene_id | grep -v N_unmapped | grep -v N_multimapping | grep -v N_noFeature | grep -v N_ambiguous | awk '{print $2}' > $folder_name$file_name_2".gene_name.txt"
    cat $file_path | grep -v GENCODE | grep -v gene_id | grep -v N_unmapped | grep -v N_multimapping | grep -v N_noFeature | grep -v N_ambiguous | awk '{print $3}' > $folder_name$file_name_2".gene_typec.txt"
    cat $file_path | grep -v GENCODE | grep -v gene_id | grep -v N_unmapped | grep -v N_multimapping | grep -v N_noFeature | grep -v N_ambiguous | awk '{print $4}' > $folder_name$file_name_2".unstranded.txt"
    cat $file_path | grep -v GENCODE | grep -v gene_id | grep -v N_unmapped | grep -v N_multimapping | grep -v N_noFeature | grep -v N_ambiguous | awk '{print $5}' > $folder_name$file_name_2".stranded_first.txt"
    cat $file_path | grep -v GENCODE | grep -v gene_id | grep -v N_unmapped | grep -v N_multimapping | grep -v N_noFeature | grep -v N_ambiguous | awk '{print $6}' > $folder_name$file_name_2".stranded_second.txt"
    cat $file_path | grep -v GENCODE | grep -v gene_id | grep -v N_unmapped | grep -v N_multimapping | grep -v N_noFeature | grep -v N_ambiguous | awk '{print $7}' > $folder_name$file_name_2".tpm_unstranded.txt"
    cat $file_path | grep -v GENCODE | grep -v gene_id | grep -v N_unmapped | grep -v N_multimapping | grep -v N_noFeature | grep -v N_ambiguous | awk '{print $8}' > $folder_name$file_name_2".fpkm_unstranded.txt"    
done
# Finish here, to do
# Add collumn to the table
# Check order of gene_ids
paste /home/felipe/googledrive/Cancer_staging/tables/*/*."gene_id.txt"  > /home/felipe/googledrive/Cancer_staging/tables/gene_id.order.txt

# I stopped here 
ls /home/felipe/googledrive/Cancer_staging/tables/*/*.rna_seq.augmented_star_gene_counts.tsv | sed 's/\// /g' | sed 's/.rna_seq.augmented_star_gene_counts.tsv/ /g' | awk '{print $7}'  | sed -e 's/ \t/\t/g'  > /home/felipe/googledrive/Cancer_staging/tables/header.txt
paste /home/felipe/googledrive/Cancer_staging/tables/*/*."gene_id.txt" | awk '{print $1}' >   /home/felipe/googledrive/Cancer_staging/tables/gene_ids.txt
paste /home/felipe/googledrive/Cancer_staging/tables/*/*."gene_name.txt" | awk '{print $1}' > /home/felipe/googledrive/Cancer_staging/tables/gene_name.txt

# Save all rnaseq counts
paste /home/felipe/googledrive/Cancer_staging/tables/*/*."unstranded.txt"      > /home/felipe/googledrive/Cancer_staging/tables/unstranded.rna_seq.augmented_star_gene_counts.tsv
paste /home/felipe/googledrive/Cancer_staging/tables/*/*."stranded_first.txt"  > /home/felipe/googledrive/Cancer_staging/tables/stranded_first.rna_seq.augmented_star_gene_counts.tsv
paste /home/felipe/googledrive/Cancer_staging/tables/*/*."stranded_second.txt" > /home/felipe/googledrive/Cancer_staging/tables/stranded_second.rna_seq.augmented_star_gene_counts.tsv
paste /home/felipe/googledrive/Cancer_staging/tables/*/*."tpm_unstranded.txt"  > /home/felipe/googledrive/Cancer_staging/tables/tpm_unstranded.rna_seq.augmented_star_gene_counts.tsv
paste /home/felipe/googledrive/Cancer_staging/tables/*/*."fpkm_unstranded.txt" > /home/felipe/googledrive/Cancer_staging/tables/fpkm_unstranded.rna_seq.augmented_star_gene_counts.tsv
