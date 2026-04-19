
cran_pkgs <- c('ggplot2','Seurat', 'remotes', 'devtools','patchwork','reticulate','harmony','BiocManager')
install.packages(cran_pkgs)
BiocManager::install("batchelor")
devtools::install_github('satijalab/seurat-data')
remotes::install_github('satijalab/seurat-wrappers')
remotes::install_github("mojaveazure/seurat-disk")
devtools::install_version('uwot', version = '0.1.10', repos = 'http://cran.us.r-project.org')





