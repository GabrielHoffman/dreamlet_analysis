
# Compute GC content for each gene

# Compute GC content
# adapted from https://github.com/dpryan79/Answers/blob/master/SEQanswers_42420/GTF2LengthGC.R
library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)
library(R.utils)

# Specify GTF and FASTA files
GTFfile = "/sc/arion/projects/psychAD/ref_GRCh38/Homo_sapiens.GRCh38.104.gtf"
FASTAfile = "/sc/arion/projects/psychAD/ref_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

GTF <- import.gff(GTFfile, format="gtf", genome="GRCh38.104", feature.type="exon")
grl <- GenomicRanges::reduce(split(GTF, elementMetadata(GTF)$gene_id))
reducedGTF <- unlist(grl, use.names=T)
elementMetadata(reducedGTF)$gene_id <- rep(names(grl), elementNROWS(grl))

# Open the fasta file
FASTA <- FaFile(FASTAfile)
open(FASTA)

# Add the GC numbers
elementMetadata(reducedGTF)$nGCs <- letterFrequency(getSeq(FASTA, reducedGTF), "GC")[,1]
elementMetadata(reducedGTF)$widths <- width(reducedGTF)

# Create a list of the ensembl_id/GC/length
calc_GC_length <- function(x) {
    nGCs = sum(elementMetadata(x)$nGCs)
    width = sum(elementMetadata(x)$widths)
    c(width, nGCs/width)
}
df_GC <- t(sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_id), calc_GC_length))
colnames(df_GC) <- c("Length", "GC")
df_GC = as.data.frame(df_GC)
df_GC$ENSEMBL = rownames(df_GC)

i = match(df_GC$ENSEMBL, GTF$gene_id)
df_GC$SYMBOL = GTF$gene_name[i]
rownames(df_GC) = c()

# write data.frame to file and gzip
file = "GRCh38.104_gc_content.tsv"
write.table(df_GC, file, quote=FALSE, sep="\t", row.names=FALSE)
gzip(file, overwrite=TRUE)


