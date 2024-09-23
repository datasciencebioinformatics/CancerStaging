# Consitency - check filters meticulously.
# FPKM, TPM  - take these as robust.

# Construnction of 3d coordinates.
# X=T2           : Thymine composition in second codon position (T2)
# Y=connectivity : interactome_data                 # Carels checked    .
# Z=expression   : reads_count_all_projects         # Use the fpkm
# Reading the contents of TSV file using read_tsv() method
Interactomes_GC3_T2_file <-"/home/felipe/Documents/github/CancerStaging/Interactomes_GC3_T2.csv"

# Read data
Interactomes_GC3_T2_data <-read.table(file = Interactomes_GC3_T2_file, sep = '\t', header = TRUE,fill=TRUE) 

# Accross cancer types.
# The hypohtoses is to put the clusters ordered by entropy in the 3d maps                   .
# The hypotheseis is to order the clusters by entropy in the 3d map                         .
# Amont the 2d coordinates this can have low signal, in a three map can have a stroger sinal.


# October, paper article.
# October.



