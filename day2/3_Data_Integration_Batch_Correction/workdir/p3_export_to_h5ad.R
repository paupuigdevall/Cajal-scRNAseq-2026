
library(Seurat)
library(SeuratDisk)

alldata_to_h5ad <- readRDS("../data/covid/results/seurat_covid_to_be_exported_h5ad.rds")
# ## By default, "h5ad" file already exists in the repository. We force to overwrite it with the conversion.
SaveH5Seurat(alldata_to_h5ad, filename = "../data/covid/results/alldata.h5Seurat", overwrite = T)
Convert("../data/covid/results/alldata.h5Seurat", dest = "h5ad", overwrite = T)
