### Script for clustering and finding differentially expressed markers ###

library(Seurat) ## Seurat v2.2.1
library(dplyr)
library(Matrix)

args=commandArgs(TRUE)

filename=args[1]
output="~/Documents/"
load(paste0(filename,".Rdata"))

### save old cluster identities, if necessary (i.e if subclustering and want to track clusters)
seur_ob <- StashIdent(object = seur_ob, save.name = "Old_Names")

## find clusters using signficant PCs deteremined in seurat_v2.R script
dim1 = 1
dim2= 15
seur_ob <- FindClusters(object =seur_ob, reduction.type = "pca", dims.use = dim1:dim2, 
                     resolution = 0.8, print.output = 0, save.SNN = TRUE)

seur_ob <- RunTSNE(object = seur_ob, dims.use = dim1:dim2, do.fast = TRUE)

## print new TSNE
png(paste0(output,filename,"res08_tsne.png"))
TSNEPlot(object = seur_ob,do.label=T)
dev.off()

## print TSNE labeled by previous names, if used
png(paste0(output,filename,"_oldnames_tsne.png"))
TSNEPlot(object = seur_ob,do.label=T,group.by="Old_Names")
dev.off()

## save seurat object
save(seur_ob, file=paste0(output,filename,".Rdata"))

## Find differentially expressed genes
get_markers = function(srt){
  seur.markers <- FindAllMarkers(object = srt, only.pos = TRUE, min.pct = 0.25, 
                               thresh.use = 0.25)
  df <- NULL

  for(x in unique(seur.markers$cluster)){clusx <- seur.markers[which(seur.markers$cluster==x),]
  sorted_clus <- clusx[order(clusx$avg_logFC,decreasing=TRUE),]
  write.csv(sorted_clus, file=paste0(output,filename,"_clus.",x,"_res08_markers.csv"))
  df <- rbind(df, sorted_clus)}
  write.csv(df,file=paste0(output,filename,"res08_all_markers.csv"))
}

get_markers(seur_ob)

