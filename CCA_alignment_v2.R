## script for aligning v2 datasets based on metagenebicor plot created in cca_merging_v2.R

library(Seurat) ## version v2.2
library(dplyr)
library(Matrix)

dir = "~/Document/" #output directory

load("path_to_object_from_cca_merging_v2.R")
seur_ob = named_object
dims= 1:10 # select number of dims from metagenebicor plot
### select CCs for downstream alignment
# Run rare non-overlapping filtering
seur_ob <- CalcVarExpRatio(object = seur_ob, reduction.type = "pca",
                                       grouping.var = "Batch", dims.use = dims)

seur_ob <- SubsetData(seur_ob, subset.name = "var.ratio.pca",
                                  accept.low = 0.5)

# Alignment
seur_ob <- AlignSubspace(seur_ob, reduction.type = "cca",grouping.var = "Batch",
                                     dims.align = dims)

# t-SNE and Clustering
seur_ob <- FindClusters(seur_ob, reduction.type = "cca.aligned",
                                    dims.use = dims, save.SNN = T, resolution = 0.8)
seur_ob <- RunTSNE(seur_ob,
                               reduction.use = "cca.aligned",
                               dims.use = dims)

png(paste0(dir,"tsne_batchlabels.png"))
TSNEPlot(seur_ob, do.label=F,group.by="Batch")
dev.off()

png(paste0(dir,"res08_tsne.png"))
TSNEPlot(seur_ob, do.label=T)
dev.off()

save(seur_ob, file=paste0(dir,"cca_aligned.Rdata"))

get_markers = function(srt){
  seur.markers <- FindAllMarkers(object = srt, only.pos = TRUE, min.pct = 0.25,
                               thresh.use = 0.25)

  df <- NULL

  for(x in unique(seur.markers$cluster)){clusx <- seur.markers[which(seur.markers$cluster==x),]
  sorted_clus <- clusx[order(clusx$avg_logFC,decreasing=TRUE),]
  write.csv(sorted_clus, file=paste0(dir,"cca_aligned_clus",x,"_res08_markers.csv"))
  df <- rbind(df, sorted_clus)}
  write.csv(df,file=paste0(dir,"cca_aligned_res08_all_markers.csv"))
}

get_markers(seur_ob)
