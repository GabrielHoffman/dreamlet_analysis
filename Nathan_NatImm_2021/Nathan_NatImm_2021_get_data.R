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

# set working directory for input and output files
wd = "/sc/arion/projects/CommonMind/hoffman/scRNAseq_data/Nathan_NatImm_2021"
setwd(wd)

# Download data from GEO
#-----------------------
cmd1 = "wget --no-check-certificate https://ftp.ncbi.nlm.nih.gov/geo/series/GSE158nnn/GSE158769/suppl/GSE158769%5Fmeta%5Fdata%2Etxt%2Egz"

cmd2 = "wget --no-check-certificate https://ftp.ncbi.nlm.nih.gov/geo/series/GSE158nnn/GSE158769/suppl/GSE158769%5Fexprs%5Fraw%2Etsv%2Egz"

system(cmd1)
system(cmd2)

# Read data
#----------

# read metadata
df_meta = as.data.frame(fread('GSE158769_meta_data.txt.gz'))
rownames(df_meta) = df_meta$cell_id

# read counts as sparseMatrix
# This uses > 50Gb memory. A more efficient implemenation
#    can dramatically reduce this
counts = read.sparseMatrix("GSE158769_exprs_raw.tsv.gz", 1000)

gc() # free memory

# Convert to SingleCellExperiment and save as H5AD
#-------------------------------------------------

# create SingleCellExperiment
sce = SingleCellExperiment( assays = list(counts = counts),
                            colData = df_meta[colnames(counts),])

# write to h5ad
outfile = paste0(wd, '/Nathan_NatImm_2021.h5ad')
writeH5AD(sce, outfile, compression="lzf", verbose=TRUE)
