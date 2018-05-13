library(monocle) ## version 2.6
library(Seurat) ## version 2.3.0

output= "~/Documents/"
filename = "example"

## rename seurat object to seur_ob

seur_ob = seurat_object

seur_ob = UpdateSeuratObject(seur_ob)

## create CDS from seurat object data
# expression Matrix - select only cells of interest from seurat_ob@data
exprs <- seur_ob@raw.data
cells <- colnames(seur_ob@data)
exprs_mat <- exprs[,match(cells, colnames(exprs))]

# check number of cells from both
print("raw data dim:")
print(dim(exprs))
print("filtered data dim:")
print(dim(exprs_mat))

## create phenodata from meta data (@data.info for seurat v1)
pheno.data = seur_ob@meta.data

## create feature data
genes <- data.frame(gene_short_name = rownames(exprs))
rownames(genes) <- rownames(exprs)

# Make CellDataSet object
pd <- new("AnnotatedDataFrame", data=pheno.data)
fd <- new("AnnotatedDataFrame", data=genes)
HSMM_expr_matrix <- exprs_mat

cds <- newCellDataSet(as(HSMM_expr_matrix,"sparseMatrix"),
  phenoData=pd,
  featureData=fd,
  lowerDetectionLimit=0.5,
  expressionFamily=negbinomial.size())


HSMM <- estimateSizeFactors(cds)
HSMM <- estimateDispersions(HSMM)

## use seurat determined variable genes for ordering
seurat_var_genes = seur_ob@var.genes

HSMM_seur_var = setOrderingFilter(HSMM, seurat_var_genes)

pdf(paste(output,filename,"_seur_var_genes.pdf",sep=""))
plot_ordering_genes(HSMM_seur_var)
dev.off()

## reduce dimensionality with DDRTree, regression on nUMI
HSMM_seur_var <- reduceDimension(HSMM_seur_var, reduction_method="DDRTree",max_components = 2,residualModelFormulaStr = "~nUMI")

## order cells
HSMM_seur_var <- orderCells(HSMM_seur_var)

## print trajectories
pdf(paste(output,filename,"_monocle_traj.pdf",sep=""))
plot_cell_trajectory(HSMM_seur_var, color_by = "Pseudotime",cell_size = 4,show_branch_points = T)
dev.off()

save(HSMM_seur_var, file=paste(output,filename,"_monocle.Rdata",sep=""))


