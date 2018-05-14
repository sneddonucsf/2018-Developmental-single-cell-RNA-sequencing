library(Seurat) #version 2.2
library(dplyr)
library(Matrix)
library(methods)

### make seurat objects from raw data in order to add unique cell IDs. Use seurat object to pull out cells of interest from raw matrix file.

filename1= "example1"
filename2="example2"
filename3 ="example3"

dir= "~/Documents/" #output directory
dir_with_raw ="path_to_raw" # path to where raw data is stored (matrix.mtx, barcodes.tsv, gene.tsv files from cellranger output)

load("path_to_dataset1")
dataset1=seur_ob

load("path_to_dataset2")
dataset2=seur_ob

load("path_to_dataset3")
dataset3=seur_ob

prepare_matrices = function(seurat_object, filename){
  ## subset out cells of interest from original matrix.mtx and barcodes.tsv files
  cell_names = seurat_object@cell.names
  barcodes = read.table(paste0(dir_with_raw,filename,"/barcodes.tsv"),sep = '\t', header = F)
 	print(head(barcodes))
  old_mat <- readMM(paste0(dir_with_raw,filename,"/matrix.mtx"))
  genes = read.table(paste0(dir_with_raw,filename,"/genes.tsv"))

  ## prepare barcodes
  if (all(grepl(pattern = "\\-1$", x = barcodes$V1))) {
    new_barcodes <- as.vector(x = as.character(x = sapply(X = as.character(barcodes$V1), 
                                                      FUN = ExtractField, field = 1, delim = "-")))
  } else{new_barcodes=as.vector(as.character(barcodes$V1))}
  
	print(which(is.na(new_barcodes)))
	print("prepared barcodes")
  ## subset barcodes of interest from seurat object
  sub_barcodes = new_barcodes[match(cell_names,new_barcodes)]
  named_barcodes =as.character(sapply(sub_barcodes,function(x) paste(filename,x,sep="_")))
  ## subset matrix from seurat object
  new_mat = old_mat[,match(cell_names,new_barcodes)]
  colnames(new_mat)=named_barcodes
  rownames(x = new_mat) <- make.unique(as.character(genes$V2))


  new_object=CreateSeuratObject(raw.data=as(as.matrix(new_mat),"dgCMatrix")) ## already selected cells of interest from seurat object so more filtering is not needed
  new_object <- NormalizeData(new_object,normalization.method = "LogNormalize",
                      scale.factor = 10000)
  new_object <- FindVariableGenes(new_object,mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
  new_object <- ScaleData(new_object,vars.to.regress = c("nUMI"))
  new_object@meta.data$Batch <- filename
  return(new_object)
}

object1=prepare_matrices(dataset1,filename1)
print("finished object1 prep")
object2=prepare_matrices(dataset2,filename2)
print("finished object2 prep")
object3=prepare_matrices(dataset3,filename3)
print("finished object3 prep")

# Determine genes to use for CCA, must be highly variable in at least 2 datasets. As shown in online tutorial.
seurat_list =c(object1,object2,object3)

genes.use = c()
for (i in 1:length(seurat_list)) {
  genes.use <- c(genes.use, head(rownames(seurat_list[[i]]@hvg.info), 1000))
}
genes.use <- names(which(table(genes.use) > 1))
for (i in 1:length(seurat_list)) {
  genes.use <- genes.use[genes.use %in% rownames(seurat_list[[i]]@scale.data)]
}

seur_ob <- RunMultiCCA(object.list = seurat_list, genes.use = genes.use, num.ccs = 40)

write.csv(genes.use, file=paste0(dir,"merged_cca_variable_genes.csv"))

# CC Selection
png(paste0(dir,"merged_Bicorplot.png"))
MetageneBicorPlot(seur_ob, grouping.var = "Batch", dims.eval = 1:40)
dev.off()

save(seur_ob,file=paste0(dir,"merged.Rdata"))


