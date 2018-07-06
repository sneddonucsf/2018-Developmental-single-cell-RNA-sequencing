### Generate seurat objects for input into monocle ###

library(Seurat)

filename = "path_to_seurat_object.Rdata"
output = "~/Documents/"

load(filename)

#rename seurat object
my.pancreas=seurat_object

## subcluster cells of interest, if desired

my.pancreas = SubsetData(my.pancreas, ident.use = c(1,2,3))

# regress on nUMI and batch, but do not scale or center
my.pancreas <- RegressOut(my.pancreas, latent.vars = c("nUMI","Labels"), do.scale=F,do.center=F)
my.pancreas <- MeanVarPlot(my.pancreas, x.low.cutoff = 0.1)
my.pancreas <- PCAFast(my.pancreas, pcs.compute=40)

save(my.pancreas, paste(output, "mon_input_seur_ob.Rdata"))