# Gabriel Hoffman
#
# Sept 11, 2023
#
# How many genes are retained after subsampling cells

# cd /sc/arion/projects/CommonMind/hoffman/dreamlet_analysis/PsychAD_r0
# ml python git
# git pull
# R --vanilla
# source("subsample_cells.R")

library(SingleCellExperiment)
library(zellkonverter)
library(dreamlet)
library(ggplot2)
library(DelayedArray)

# update block size for reading h5ad file from disk
setAutoBlockSize(1e9)

# Read H5AD
folder = "/sc/arion/projects/psychAD/NPS-AD/public_release_0/" 
file = paste0(folder, "PsychAD_r0_Dec_28_2022.h5ad")
sce = readH5AD(file, use_hdf5=TRUE, verbose=TRUE)
assayNames(sce)[1] = "counts"

df_grd = expand.grid(ncells = c(3, 5, 8, 10, 20, 25, 40, 50, 75, 100, 150, 200, 300),
						CT = "Astro")

df_counts = lapply( seq(nrow(df_grd)), function(i){

	message(i)

	# subsample 
	keep = sapply( levels(sce$Channel), function(sid){

		keep = which((sce$Channel == sid) & (sce$subclass == as.character(df_grd$CT[i])))

		# subsample
		if( length(keep) > df_grd$ncells[i]){
			keep = sample(keep, df_grd$ncells[i])
		}
		keep
		})
	keep = sort(unlist(keep))

	# Compute pseudobulk from subsetted data
	pb <- aggregateToPseudoBulk(sce[,keep],
	  assay = "counts",
	  cluster_id = "subclass",
	  sample_id = "Channel",
	  verbose = FALSE)

	# apply filtering
	res.proc <- processAssays(pb, ~ 1)

	# number of genes retained
	ngenes = nrow(assay(res.proc, 1))

	data.frame(ncells = df_grd$ncells[i], 
				CT = df_grd$CT[i], 
				ngenes = ngenes)
})
df_counts = do.call(rbind, df_counts)

# save results to file
write.table(df_counts, file="df_counts.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# plot
#######

ctorder = c('EN_L2_3_IT', 'EN_L3_5_IT_1', 'EN_L3_5_IT_2', 'EN_L3_5_IT_3', 'EN_L5_6_NP', 'EN_L6_CT', 'EN_L6_IT', 'EN_NF', 'IN_ADARB2', 'IN_LAMP5', 'IN_PVALB', 'IN_PVALB_CHC', 'IN_SST', 'IN_VIP', 'Oligo', 'OPC', 'Astro', 'Micro_PVM', 'CD8_T', 'PC', 'VLMC','Endo')

df_counts$CT = factor(df_counts$CT, ctorder)

ymax = max(df_counts$ngenes) * 1.03

fig = ggplot(df_counts, aes(ncells, ngenes, color=CT)) +
	geom_line() +
	geom_point() +
	theme_classic() +
	theme(aspect.ratio=1, legend.position="none") +
	facet_wrap(~ CT) +
	scale_y_continuous(expand=c(0,0), limits=c(0, ymax)) +
	scale_x_log10() +
	xlab("# of cells") +
	ylab("# of genes passing filter")

ggsave(fig, file="ngenes_subsampling.pdf")










