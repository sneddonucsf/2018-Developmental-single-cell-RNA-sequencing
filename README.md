# 2018-Developmental-single-cell-RNA-sequencing

Scripts used for analysis of data in manuscript under revision at Nature Communications. "Lineage Dynamics of Pancreatic Development at Single-Cell Resolution." Byrnes L. & Wong D. et al.

# Scripts description

seurat_v1.R - initialize seurat object from 10X Genomics cellranger outputs. Includes filtering, normalization, regression, variable gene identification, PCA analysis, clustering, tSNE visualization. Used for v1 datasets. <br /> <br />
merge_seurat.R - merge two or more seurat objects into one seurat object. Perform linear regression to remove batch effects from separate objects. Used for v1 datasets. <br /> <br />
seurat_v2.R - initialize seurat object from 10X Genomics cellranger outputs. Includes filtering, normalization, regression, variable gene identification, and PCA analysis. Used for v2 datasets. <br /><br />
clustering_markers_v2.R - clustering and tSNE visualization for v2 datasets. <br /><br />
seurat_object_analysis_v1_and_v2.R - downstream analysis and plotting functions for seurat object created by seurat_v1.R or seurat_v2.R. <br /><br />
merge_clusters.R - merge clusters that do not meet gene threshold. Used for both v1 and v2 datasets. <br /><br />
monocle_v1_seurat_input.R - monocle script using seurat batch corrected values as input for v1 merged timecourse datasets. <br /><br />
monocle_object_analysis.R - downstream analysis for monocle object - BEAM and plotting. <br /><br />
