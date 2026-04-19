bioc_pkgs <- c('Seurat', 'scran', 'ggtree', 'MultiAssayExperiment', 'Biostrings', 'XVector',
               'SingleCellExperiment', 'SummarizedExperiment', 'GenomicRanges', 'GenomeInfoDb', 'IRanges', 'S4Vectors',
               'BiocGenerics', 'phyloseq', 'scran', 'tidySummarizedExperiment', 'HGNChelper', 'scDblFinder',
               'clusterProfiler', 'glmGamPoi')
BiocManager::install(bioc_pkgs)
BiocManager::install('org.Hs.eg.db', character.only = TRUE)
devtools::install_version('uwot', version = '0.1.10', repos = 'http://cran.us.r-project.org')
devtools::install_github('chris-mcginnis-ucsf/DoubletFinder@3b420df6', upgrade = FALSE)
remotes::install_github('https://github.com/chris-mcginnis-ucsf/DoubletFinder@3b420df68b8e2a0cc6ebd4c5c1c7ea170464c97f', upgrade = FALSE, dependencies = FALSE)
remotes::install_github("mjemons/spatialFDA")
remotes::install_github("sgunz/sosta")
