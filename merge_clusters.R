### Determine the number of differentially expressed genes for each cluster. Merge clusters that do not meet >=9 genes with >2 fold change threshold ###

### used for both v1 and v2 datasets (need to change name of fold change column in differentially expressed genes if using seurat v1 - avg_diff instead of avg_logFC. Noted in script)

library(Seurat)
library(dplyr)

args=commandArgs(TRUE)

filename = args[1]
output=args[2]

### need seurat object and differentially expressed genes for all clusters (from FindAllMarkers)
load(filename)
markers <- read.csv("seurat_object_markers_res08.csv", header=T)

## rename seurat object 
obj <- seurat_object_name

### start seurat object off at desired resolution
TSNEPlot(obj, do.label=T)
obj <- FindClusters(obj, resolution =2, print.output = F)

# sort by avg_diff

sorted_markers <- markers[order(markers$cluster,markers$avg_logFC,decreasing=TRUE),] ## for seurat v1 - change avg_logFC to avg_diff

# filter by > log

filtered_markers <- sorted_markers %>% filter(avg_logFC > log(2)) ## for seurat v1 - change avg_logFC to avg_diff

# make table of total differentially expressed genes

marker_table <- table(filtered_markers$cluster)
df <-as.data.frame(marker_table)
df
write.table(df,file=paste(output,"_first_diff_genes_table.txt",sep=""))

# build phylogenetic tree to choose which clusters to merge with which
pdf(paste(output,"_first_tree.pdf",sep=""))
obj <- BuildClusterTree(obj)
dev.off()

# make column in data.info of seurat_ob with merged clusters with less than 10 differentially expressed genes
merged_clus <- obj@ident
merged_clus[merged_clus == c(2)] <- 1

# place back in seurat obj
merged_df <- data.frame(merged = merged_clus)
rownames(merged_df) <- rownames(obj@data.info)
obj <- AddMetaData(obj,metadata=merged_df )

# view TSNE Plot
pdf(paste(output, "_first_merge.pdf",sep=""),width=7,height=5)
TSNEPlot(obj, do.label=T, group.by="merged")
dev.off()

# print new differentially expressed genes

obj <- SetAllIdent(obj, id = "merged")
new_markers <- FindAllMarkers(obj,only.pos=F, min.pct=0.25, thresh.use=0.25)
write.csv(new_markers, file=paste(output,"_first_merge.csv",sep=""))

markers <- new_markers

### check pairwise markers for those with similar markers

pairwise=FindMarkers(obj,1,7,pos.only=F)
pairwise=pairwise[order(pairwise$avg_logFC,decreasing=T),] ## for seurat v1 - change avg_logFC to avg_diff
write.csv(pairwise, file=paste(output,"_pairwise_1_vs_7_markers.csv",sep=""))

save(obj, file=paste(output,"_first_merge.Rdata",sep=""))

## can now repeat from top of script with next clusters if needed