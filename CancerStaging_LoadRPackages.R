# Function to expand.grid.unique without redundancy
expand.grid.unique <- function(x, y, include.equals=FALSE)
{
    x <- unique(x)
    y <- unique(y)
    g <- function(i)
    {
        z <- setdiff(y, x[seq_len(i-include.equals)])

        if(length(z)) cbind(x[i], z, deparse.level=0)
    }
    do.call(rbind, lapply(seq_along(x), g))
}
# Calculating Z score in R
calculate_z <- function(X, X_mean, S)
{
  return((X-X_mean)/S)
}

# Packages and data use throught 
library(ggplot2)
library("readr")
library("reshape2")
library(ggfortify)
#library(AnnotationHub)
#library (edgeR)
#library (EDASeq)
#library("biomaRt")
#library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#library("ggplot2")
#library("GGally")
#library("goseq")
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#library("clusterProfiler")
#library("org.Hs.eg.db")
#library("ggplot2")
#library("data.table")
#library("plotly")
#library(GGally)
#library("scatterplot3d") # load
library("ggpubr")
library("gridExtra")
#library("ggVennDiagram")



#library(metR)
#library("ggpairs")
#library(AnnotationHub)
#library("biomaRt")
#library(data.table)
#library(DescTools)
#library("DESeq2")
#library (EDASeq)
#library("enrichplot")
#library (edgeR)
#library(ggfortify) 
#library(ggplot2)
#library(ggplot2)  
#library(ggVennDiagram)
#library(gridExtra)
#library("enrichplot")
#library("GWENA")
#library(igraph)
#library("readr")
#library(readr)
#library(readr)                                                                                                            #
#library("readxl")
#library(viridis)
#library("xlsx")
#library("xlsx")
#library("gtools")
#library(tidyr)
#library(dplyr)
#library("clusterProfiler")
#library("org.Hs.eg.db")
#library("org.Hs.eg.db")
#library(RCy3)
#library(RColorBrewer)
#library("clusterProfiler")
#library(ReactomePA)
#library(stringr)
#library("amap")
#library("ggpubr")
