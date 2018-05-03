#### Merge two or more Seurat objects. Used for v1 datasets with seurat v1. Linear regression used to remove batch effects between seurat objects. #####

library(Seurat) ## v1.4
library(dplyr)
library(methods)

output="~/Documents/"
filename = "MyData"


################### Merging ###########################

### Load in Objects for Merging ###

load("path_to_object_1.Rdata")
# rename object
E14_batch1 <- my.pancreas

load("path_to_object_2.Rdata")
E14_batch2 <- my.pancreas

### Merge objects

my.pancreas <- MergeSeurat(E14_batch1,E14_batch2, do.scale=F,do.center=F) #will be scaled and centering in RegressOut function below.

### make Labels. If named each seruat object something different in "project" parameterin Setup function (seurat_v1.R), can also use that column for regression rather than making separate column.

E14_B1_labels <- rep("E14_B1", dim(E14_batch1@data)[2])
E14_B2_labels <- rep("E14_B2",dim(E14_batch2@data)[2])

# name of column will be used in RegressOut
labels <- data.frame(Labels = unlist(c(E14_B1_labels,E14_B2_labels)))
rownames(labels) <- rownames(my.pancreas@data.info)

my.pancreas <- AddMetaData(my.pancreas, metadata=labels)

### Reanalyze merged dataset. Regress on nUMI and name given above ("Labels")  ###

analyze <- function(x){	
  x <- MeanVarPlot(x, x.low.cutoff = 0.1)
  x <- RegressOut(x, latent.vars = c("nUMI", "Labels"), genes.regress = x@var.genes)
  x <- PCAFast(x, pcs.compute = 40)
  return(x)
}

my.pancreas=analyze(my.pancreas)

### Choose PCAs to use in downstream analysis by scree plot
png(paste(output,filename,"_pcaplot.png",sep=""))
PCElbowPlot(my.pancreas, num.pc=40)
dev.off()

### Run TSNE and find clusters with determined number of PCs (TSNE and FindClusters should match)
my.pancreas <- RunTSNE(my.pancreas, dims.use = 1:10, do.fast = T)
my.pancreas <- FindClusters(my.pancreas, pc.use = 1:10, save.SNN = T, do.sparse = T)

## TSNE Plot of reanalyzed Data
TSNEPlot(my.pancreas,do.label = T,group.by="res.0.8")

## Write out object for downstream analysis

save(my.pancreas, file=paste(output,filename,".Rdata",sep=""))

