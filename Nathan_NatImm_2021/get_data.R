# Gabriel Hoffman
# 
# September 6, 2022
#
# Download and format data from Nathan, et al. 2021

library(data.table)
library(SingleCellExperiment)
library(zellkonverter)
library(stringr)
source("read.sparseMatrix.R")

# Download data from GEO
#-----------------------
cmd1 = "wget --no-check-certificate https://ftp.ncbi.nlm.nih.gov/geo/series/GSE158nnn/GSE158769/suppl/GSE158769%5Fmeta%5Fdata%2Etxt%2Egz"

cmd2 = "wget --no-check-certificate https://ftp.ncbi.nlm.nih.gov/geo/series/GSE158nnn/GSE158769/suppl/GSE158769%5Fexprs%5Fraw%2Etsv%2Egz"

system(cmd1)
system(cmd2)

# Read data
#----------

setwd("/sc/arion/projects/CommonMind/hoffman/scRNAseq_data/Nathan_NatImm_2021")

# read metadata
df_meta = as.data.frame(fread('GSE158769_meta_data.txt.gz'))
rownames(df_meta) = df_meta$cell_id

# read counts as sparseMatrix
# This uses > 50Gb memory. A more efficient implemenation
#    can dramatically reduce this
counts = read.sparseMatrix("GSE158769_exprs_raw.tsv.gz", 100)

# free memory
gc() 

# Convert to SingleCellExperiment and save as H5AD
#-------------------------------------------------

# create SingleCellExperiment
sce = SingleCellExperiment( assays = list(counts = counts),
                            colData = df_meta[colnames(counts),])

# write to h5ad
outpath = "/sc/arion/projects/CommonMind/hoffman/scRNAseq_data/Nathan_NatImm_2021/"
outfile = paste0(outpath, '/Nathan_NatImm_2021.h5ad')
writeH5AD(sce, outfile, compression="lzf", verbose=TRUE)
