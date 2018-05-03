##### Analyze Seurat Object (or Subclustered Seurat Object) ####
### sections do not necessarily need to be run in order

library(Seurat) ## version 2.3.0. Can still be used with seurat v1, some default settings change 
library(dplyr)

#### Inputs for Program ####

input_file <- "path_to_seurat_object.Rdata"
output="~/Documents/"

load(input_file)
## rename seurat object
seurat_ob <- seurat_object ## if using object created in seurat v1, can update object with seurat_ob=UpdateSeuratObject(seurat_object)
#######################################################################################

########## Choose Resolution for downstream analysis, if required #########

seurat_ob <- SetAllIdent(seurat_ob, id = "res.0.8")

######### PCA Analysis ############

png(paste(output,"pca1_2.png",sep=""))
PCAPlot(seurat_ob,1,2)
dev.off()

pdf(paste(output,"pcaheatmap_top10.pdf",sep=""))
PCHeatmap(seurat_ob, pc.use = 1:10, cells.use=500,do.balanced = TRUE, cexRow=0.5)
dev.off()

ProjectPCA(seurat_ob)
PrintPCA(seurat_ob, pcs.print = 1:2, use.full=T)

png(paste(root_dir,file_name, "_vizpca_1_2.png",sep=""))
VizPCA(seurat_ob, 2:3)
dev.off()

########## Plot TSNE results ##########

png(paste(output, "_res08.png",sep=""), width=695, height=538)
TSNEPlot(seurat_ob,do.label=T,pt.size=3,do.return=T, label.size=8)
dev.off()


########## Find Markers of Clusters ################

## Find Cluster Markers ##
### seurat v1 - default is "bimod"
### seurat v2 - default is "wilcox", "MAST" is also available
### can alter with parameter test.use

all_markers <- FindAllMarkers(seurat_ob,only.pos=T, min.pct=0.25, thresh.use=0.25)
write.csv(all_markers, file=paste(output, "_res08_markers.csv",sep=""))

df <- NULL
for(x in unique(all_markers$cluster)){clusx <- all_markers[which(all_markers$cluster==x),]
sorted_clus <- clusx[order(clusx$avg_diff,decreasing=TRUE),]
df <- rbind(df, sorted_clus)}

write.csv(df, file=paste(output,"_res08_sorted_markers_avgdiff.csv",sep=""))

# Find Pairwise Markers #

## between two clusters
markers <- FindMarkers(seurat_ob,2,3,min.pct=0.25, thresh.use=0.25, only.pos=F)
sorted_markers <- markers[order(markers$avg_diff,decreasing=TRUE),]

## 1 cluster vs. 2 others
markers <- FindMarkers(seurat_ob,1,c(2:3),min.pct=0, thresh.use=0, only.pos=F)
write.csv(sorted_markers, file=paste(output, "_clus1_clus2_3_avgdiff_markers.csv",sep=""))

######## Feature and Violin Plots for Differentially Expressed Genes #######

### load in already calculated markers if needed ###
all_markers <- read.csv("path_to_markers.csv") ## from FindAllMarkers

sorted <- all_markers %>% group_by(cluster) %>% top_n(50, avg_diff)

genes <- unique(sorted$gene)

endo_genes=c("Gcg","Ins1","Fev","Neurog3")

feature_plots <- function(genes){
  for(x in genes){
    pdf(paste(output, "feature_plot_",x,".pdf",sep=""), width=7, height=5,useDingbats = F)
    FeaturePlot(seurat_ob,x,cols.use = c("gray","red"), pt.size=2,no.legend=T,no.axes=T)
    dev.off()
  }
}

feature_plots(endo_genes)
feature_plots(genes)

violin_plots <- function(genes){
  for(x in genes){
    png(paste(output, "_violin_",x,".png",sep=""),width=1000,height=200)
    VlnPlot(seurat_ob, x,size.x.use = 12,size.y.use=10,size.title.use=20)
    dev.off()
  }
}

violin_plots(endo_genes)
violin_plots(genes)

######## Dot Plots ###########

pdf(paste(output, "_dotplot.pdf"),width=15, height=4)
DotPlot(seurat_ob,genes,cols.use=myPalette(low = "blue",high="red"), cex.use = 2)
dev.off()

########################### Change TSNE Resolution ###############################

seurat_ob <- StashIdent(seurat_ob, save.name = "orig.res")
seurat_ob <- FindClusters(seurat_ob,resolution=1,print.output = F)

png(paste(output, "_res1.png",sep=""),width=695, height=538)
	TSNEPlot(seurat_ob,
		do.label=T,
		do.return=T, 
		pt.size=2, 
		label.size=8,
		no.legend=F)
	dev.off()

	
### get cell names of particular cluster ###
	
new_ob = SubsetData(seurat_ob,ident.use = 15)

new_cell_names = colnames(new_ob@data)

########################## Adding Metadata ##################################

current.cluster.ids <- levels(seurat_ob@ident)
new.cluster.ids <- c("Acinar","Mature Acinar","Prolif. Acinar","Prolif. Ductal","Ductal","Ngn3","Fev","Beta","Alpha","Epsilon")
seurat_ob@ident <- plyr::mapvalues(seurat_ob@ident, from = current.cluster.ids, to = new.cluster.ids)

meta <- data.frame(Cluster_Names=seurat_ob@ident)

seurat_ob <- AddMetaData(seurat_ob,metadata=meta)


save(seurat_ob, file=paste(output,"analysis.Rdata",sep=""))















