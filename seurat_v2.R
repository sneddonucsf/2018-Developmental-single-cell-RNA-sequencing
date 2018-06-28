#### initialize seurat object and perform analysis to tsne visualization ####
## uses Seurat v2
## used for v2 datasets: E12.5, E14.5, E17.5 and E14.5 Fev-Cre;mTmG

library(Seurat)
library(dplyr)
library(Matrix)

args=commandArgs(TRUE)

## output directory, filename info
filename = args[1]
output=args[2]

#read in data
pancreas.data <- Read10X(filename)

## make seurat object
seur_ob <- CreateSeuratObject(raw.data = pancreas.data, min.cells = 3, min.genes = 200, 
                           project = filename)
## add mito data to meta.data
mito.genes <- grep(pattern = "^mt-", x = rownames(x = seur_ob@data), value = TRUE)
percent.mito <- Matrix::colSums(seur_ob@raw.data[mito.genes, ])/Matrix::colSums(seur_ob@raw.data)

seur_ob <- AddMetaData(object = seur_ob, metadata = percent.mito, col.name = "percent.mito")

png(paste0(filename,"_nGene_nUMI_mito_plot.png"))
VlnPlot(object = seur_ob, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()

## visualize nUMI, genes, % mito

par(mfrow = c(1, 2))
png(paste0(filename,"_QCplots.png"))
GenePlot(object = seur_ob, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = seur_ob, gene1 = "nUMI", gene2 = "nGene")
dev.off()

#normalize data
seur_ob <- NormalizeData(object = seur_ob, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

##variable gene identification
## use cut-offs suggested by online tutorial
par(mfrow = c(1, 1))
seur_ob <- FindVariableGenes(object = seur_ob, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(seur_ob@var.genes)

#scaling and regression
seur_ob <- ScaleData(object = seur_ob, vars.to.regress = c("nUMI"))

## dimension reduction
seur_ob <- RunPCA(object = seur_ob, pc.genes = seur_ob@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5,pcs.compute=40)

## choose sig PCs based on scree plot as outlined in tutorial
png(paste0(filename,"_PCelbow.png"))
PCElbowPlot(object = seur_ob)
dev.off()

save(seur_ob, file=paste0(filename,".Rdata"))
