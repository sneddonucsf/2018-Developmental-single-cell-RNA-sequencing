# 2018-Developmental-single-cell-RNA-sequencing

Scripts used for analysis of data in manuscript under revision at Nature Communications. "Lineage Dynamics of Pancreatic Development at Single-Cell Resolution." Byrnes L. & Wong D. et al.

# Scripts description

seurat_v1.R - initialize seurat object from 10X Genomics cellranger outputs. Used for v1 datasets. <br />
seurat_v2.R - initialize seurat object from 10X Genomics cellranger outputs. Used for v2 datasets.
seurat_object_analysis_v1_and_v2.R - downstream analysis for seurat object created by seurat_v1.R or seurat_v2.R. 
merge_clusters.R - merge clusters that do not meet gene threshold. Used for both v1 and b2 datasets.
monocle_v1_seurat_input.R - monocle script using seurat batch corrected values as input for v1 merged timecourse datasets.
monocle_object_analysis.R - downstream analysis for monocle object - BEAM and plotting
