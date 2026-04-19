
cran_pkgs <- c('ggplot2','Seurat', 'remotes', 'devtools','patchwork','tidyverse','RColorBrewer','clustree','pheatmap')
install.packages(cran_pkgs)
devtools::install_version('uwot', version = '0.1.10', repos = 'http://cran.us.r-project.org')
devtools::install_github("davidsjoberg/ggsankey")

