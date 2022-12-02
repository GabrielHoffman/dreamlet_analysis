
# Download data
###############
cd /sc/arion/projects/CommonMind/hoffman/crumblr_analysis/Velmeshev_2019/data

wget --no-check-certificate https://cells.ucsc.edu/autism/rawMatrix.zip
wget --no-check-certificate https://cells.ucsc.edu/autism/meta.tsv
wget --no-check-certificate https://cells.ucsc.edu/autism/tMinusSNE.coords.tsv.gz

unzip rawMatrix.zip
gzip matrix.mtx

# Load data into R
##################

library(SingleCellExperiment)
library(zellkonverter)
library(Matrix)

counts = readMM("matrix.mtx")
barcodes = read.table("barcodes.tsv")

genes = read.table("genes.tsv")
colnames(genes) = c("ENSEMBL", "SYMBOL")
rownames(genes) = gene$ENSEMBL

colnames(counts) = barcodes$V1
rownames(counts) = genes$ENSEMBL

df_tsne = read.table("tMinusSNE.coords.tsv.gz", row.names=1)
colnames(df_tsne) = c( "tsne1", "tsne2")
identical(colnames(counts), rownames(df_tsne))

meta = read.table("meta.tsv", sep="\t", header=TRUE)
rownames(meta) = meta$cell
meta$cell = c()
meta$individual = factor(meta$individual)
identical(colnames(counts), rownames(meta))

sce = SingleCellExperiment( assays = list(counts = counts),
							colData = meta,
							rowData = genes,
							reducedDims = list(tSNE = df_tsne)  )

sce$diagnosis = factor(sce$diagnosis, c("Control", "ASD"))

writeH5AD(sce, file="Velmeshev_2019.h5ad", compression="lzf", verbose=TRUE)




