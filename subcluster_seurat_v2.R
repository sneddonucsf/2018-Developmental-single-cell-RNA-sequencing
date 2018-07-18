### Subcluster seurat object to only contain clusters of interest ####
### Seurat V2 ####

library(Seurat)
library(dplyr)
library(Matrix)

setwd("~/Documents/")

load("path_to_seur_ob.Rdata")

seur_ob = SubsetData(seur_ob,ident.use = c(0,1))

seur_ob <- FindVariableGenes(object = seur_ob, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

seur_ob <- ScaleData(object = seur_ob, vars.to.regress = c("nUMI"))

seur_ob <- RunPCA(object = seur_ob, pc.genes = seur_ob@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5,pcs.compute=40)

save(seur_ob,file="subclustered.Rdata")
