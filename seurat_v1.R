#### initialize seurat object and perform analysis to tsne visualization ####
## uses Seurat v1
## used for v1 datasets: E12.5 Batch 2, E14.5 Batch 1, E14.5 Batch 2, and E17.5 Batch 2

library(Seurat)
library(dplyr)
library(methods)

args <- commandArgs(trailingOnly=T)

input_name <- args[1]
output_dir <- args[2]

## Load in cellranger generated folder with original matrix, cell, and gene information
## includes matrix.mtx, genes.tsv, barcodes.tsv files in outs cellranger output.

my.pancreas.data <- Read10X(input_name)

## Initialize seurat object
## Filter cells, normalize
my.pancreas <- new("seurat", raw.data = sparse_counts)
my.pancreas <- Setup(my.pancreas, min.cells = 3, min.genes = 200, do.logNormalize = T, project = "PANCREAS", do.center = F, do.scale=F)

## find variable genes
my.pancreas <- MeanVarPlot(my.pancreas, x.low.cutoff = 0.1)

## regress out nUMI and/or batch effects
my.pancreas <- RegressOut(my.pancreas, latent.vars = "nUMI", genes.regress = my.pancreas@var.genes)

print(length(my.pancreas@var.genes))

## run PCA
my.pancreas <- PCAFast(my.pancreas, pcs.compute = 40)

## choose significant PCs based on scree plot as described in tutorial
png(paste(output_dir,".png",sep=""))
PCElbowPlot(my.pancreas, num.pc=40)
dev.off()

## clustering and tsne with sig PCs
## dims.use should match between RunTSNE and FindClusters.
my.pancreas <- RunTSNE(my.pancreas, dims.use = 1:15, do.fast = T)
my.pancreas <- FindClusters(my.pancreas, pc.use = 1:15, save.SNN = T, do.sparse = T)

## save object
## for loading into downstream scipts
save(my.pancreas, file=paste(output_dir,".Rdata",sep=""))





