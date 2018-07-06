### monocle for timecourse data
### use seurat normalized and transformed data as input into monocle generated in prepare_for_monocle_v1.R

library(Seurat)
library(monocle)

filename = "path_to_seurat_object.Rdata"
output = "~/Documents/"

load(filename)

####### Expression matrix #########
# rename seurat object
obj <-my.pancreas

#expression matrix - use scaled data, which is regressed, normalized
exprs <- obj@scale.data

## use variable genes determined by seurat##
var_genes <- obj@var.genes

############# PhenoData ################

## AnnotatedDataFrame object, where rows are cells, and columns are cell attributes 

pheno.data <- obj@data.info 

######### FeatureData #################
# AnnotatedDataFrame object,  where  rows  are  features  (e.g.   genes),  and  columns  are  gene attributes, such as biotype, gc content, etc

genes <- data.frame(gene_short_name = rownames(exprs))
rownames(genes) <- rownames(exprs)

###########################################################
###########################################################

# Make CellDataSet object
pd <- new("AnnotatedDataFrame", data=pheno.data)
fd <- new("AnnotatedDataFrame", data=genes)
HSMM_expr_matrix <- as.matrix(exprs)

## use gaussianff family since inputting normalized, transformed values
HSMM <- newCellDataSet(
  HSMM_expr_matrix,
  phenoData=pd,
  featureData=fd,
 # lowerDetectionLimit=0.5,
  expressionFamily=gaussianff())

ordering_genes <- var_genes

HSMM <- setOrderingFilter(HSMM, ordering_genes)
print(dim(exprs(HSMM)))

## reduce dimension - do not normalize or include pseudo count. do use monocle scaling
HSMM <- reduceDimension(HSMM,norm_method="none", reduction_method="DDRTree",max_components=2,scaling=TRUE,verbose=TRUE,pseudo_expr=0)

## order cells
HSMM <- orderCells(HSMM)

pdf(output_trajectory)
plot_cell_trajectory(HSMM, color_by="Labels")
dev.off()

save(HSMM, file=output)

