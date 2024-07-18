# Manifest files
# Filter r
head /home/felipe/googledrive/Cancer_staging/gdc_manifest.2024-07-18.txt -n 1 > /home/felipe/googledrive/Cancer_staging/gdc_manifest.2024-07-18.filtered.txt

# Keep only rna_seq.augmented_star_gene_counts
cat /home/felipe/googledrive/Cancer_staging/gdc_manifest.2024-07-18.txt	| grep -f /home/felipe/googledrive/Cancer_staging/used_file_names.tsv >> /home/felipe/googledrive/Cancer_staging/gdc_manifest.2024-07-18.filtered.txt

# Create sub-sub folder for Prostate sub-type 
mkdir /home/felipe/googledrive/Cancer_staging/tables/

# Download all files
gdc-client download -m /home/felipe/googledrive/Cancer_staging/gdc_manifest.2024-07-18.filtered.txt -d /home/felipe/googledrive/Cancer_staging/tables/
