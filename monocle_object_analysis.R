### analysis on HSMM object created in monocle scripts ###

library(monocle)
library(Seurat)

filename ="path_to_HSMM_object.Rdata"

load(filname)

#rename to HSMM
HSMM <- HSMM_object

### accessing pheno and feature data ###

pheno.data <- HSMM@phenoData@data

gene.data<-HSMM@featureData@data

### change state/end of pseudotime ###

HSMM <- orderCells(HSMM, root_state = 4)

#### Plot cell trajectory

plot_cell_trajectory(HSMM,cell_size=4,color_by="State",cell_name_size=3,show_branch_points=T)

###### branch analysis ######
## follow online tutorial

BEAM_res <- BEAM(HSMM, branch_point=6, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
write.csv(BEAM_res, file="BEAM_res_branch6.csv")

## subset based on qval cut-off
subsetted <- subset(BEAM_res, qval < 0.001)

## plot branched heatmap - choose number of gene clusters ##
x<- plot_genes_branched_heatmap(HSMM[row.names(subsetted),],
                            branch_point=2,
                            num_clusters =8,
                            cores=1,
                            use_gene_short_name = T,
                            show_rownames = F,
                            return_heatmap = T,
                            add_annotation_row = df
                            )


## heatmap sometimes has too many genes to easily visualize individual genes
## Find genes of particular cluster and plot separately, for better visualization of trend
## print out individual cluster heatmaps and list of genes in same order as those in heatmap

## make z = 1: number of clusters chosen above in branched heatmap call
## k = number of clusters chosen above
for(z in 1:8){
  hsmm.clust <-  data.frame(cluster = cutree(x$ph$tree_row, k = 8)) # choose same cluster # as above
  hsmm_clust_ordered <- hsmm.clust[order(hsmm.clust$cluster),,drop=FALSE]
  hsmm_clust_ordered <- cbind(hsmm_clust_ordered, gene_name=rownames(hsmm_clust_ordered))

  clus <- hsmm_clust_ordered[hsmm_clust_ordered$cluster==z,]
  clus_of_interest <- as.character(clus$gene_name)

  png(paste("monocle_clus_",z,"_branch2.png",sep=""))
  hmap <- plot_genes_branched_heatmap(HSMM[clus_of_interest,],
                                    branch_point=2,
                                    num_clusters =1,
                                    cores=1,
                                    use_gene_short_name = T,
                                    show_rownames = T,
                                    return_heatmap = T)
  dev.off()



  heatmap_genes = rownames(hmap$annotation_row)
  ordered_genes= heatmap_genes[hmap$ph$tree_row$order]

  final_genes <- BEAM_res[match(ordered_genes,rownames(BEAM_res)),]

  write.csv(final_genes, file=paste("monocle_clus_",z,"_branch2_heatmap_genes.csv",sep=""))
}

###### plotting individual gene trends #####

genes =c("Ins1","Gcg","Fev")
plot_genes_in_pseudotime(HSMM[genes,],color_by="State")

## chose end states
plot_multiple_branches_pseudotime(HSMM[genes,],branches = c(3,5)) 

## chose branch point
plot_genes_branched_pseudotime(HSMM[genes,],branch_point = 2)

