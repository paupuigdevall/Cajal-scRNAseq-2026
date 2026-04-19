
cran_pkgs <- c('ggplot2','enrichR','Seurat','dplyr','Matrix', 'remotes',
	       'devtools','patchwork','BiocManager','tidyverse','RColorBrewer',
	       'gghighlight','pheatmap','lme4')
install.packages(cran_pkgs)
bioconductor_pkgs <- c('edgeR','MAST','clusterProfiler','org.Hs.eg.db','enrichplot','msigdbr','fgsea')
BiocManager::install(bioconductor_pkgs)
devtools::install_version('uwot', version = '0.1.10', repos = 'http://cran.us.r-project.org')

