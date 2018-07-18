#### Seurat Package ####
#### Subcluster out particular groups of cells for reanalysis ####
### Seurat V1 ####

library(Seurat)
library(dplyr)
library(methods)

### Load in Seurat Object for Subclustering ###

load("path_to_seurat_ob.Rdata")
my.data <- seurat_ob

### Subset out cells of interest ###
clusts <- c(11,1,18,12,13,2,4,5,3,14)
my.pancreas <- SubsetData(my.data, id=clusts)

### Reanalyze: Var genes, regress out, run PCA
### add "Labels" to RegressOut function when combining batches
analyze <- function(x){	
        x <- MeanVarPlot(x, x.low.cutoff = 0.1)
        x <- RegressOut(x, latent.vars = c("nUMI"), genes.regress = x@var.genes)
        x <- PCAFast(x, pcs.compute = 40)
        return(x)
}

my.pancreas <- analyze(my.pancreas)

### Choose PCAs to use in downstream analysis
PCElbowPlot(my.pancreas, num.pc=40)

### Run TSNE and find clusters with determined number of PCs (TSNE and FindClusters should match)
my.pancreas <- RunTSNE(my.pancreas, dims.use = 1:10, do.fast = T)
my.pancreas <- FindClusters(my.pancreas ,pc.use = 1:10,save.SNN = T)

## Add old IDs from clustering done with all cells
Clusters_all_cell <- data.frame(Clusters_all_cell = my.data@ident)
my.pancreas <- AddMetaData(my.pancreas,metadata = Clusters_all_cell)

save(my.pancreas, file="subclustered.Rdata")
