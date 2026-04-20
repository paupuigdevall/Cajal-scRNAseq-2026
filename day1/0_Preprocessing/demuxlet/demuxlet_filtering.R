

## incorporation of demuxlet into data processing

.libPaths(c("/shared/projects/tp_2616_fnom_183960/conda/envs/PP_r_env_all/lib/R/library",
            "/shared/home/tp184323/R/x86_64-conda-linux-gnu-library/4.5"))

library(Seurat)
library(readxl)
library(ggplot2)
library(dplyr)
library(harmony)
library(readxl)


path_to_cellranger_output <- "/shared/projects/tp_2616_fnom_183960/Puigdevall_P/cr_out/"

qcTab <- sapply(paste0(path_to_cellranger_output,dir(path_to_cellranger_output)), function(x){
  
  summfiles <- paste0(x,"/outs/metrics_summary.csv")
  tt <- read.csv(summfiles)
  tt$analysis <- sapply(strsplit(summfiles,"/"), function(x) x[7])
  return(tt)
  
}, simplify=F)

qcTab <- do.call("rbind", qcTab); rownames(qcTab) <- NULL


metadata_sra <- read.csv("/shared/projects/tp_2616_fnom_183960/Puigdevall_P/midbrainDataset/rawData/SraRunTable.csv")
link_ids <- read.table("/shared/projects/tp_2616_fnom_183960/Puigdevall_P/midbrainDataset/rawData/README")
metadata_paper <- read_excel("/shared/projects/tp_2616_fnom_183960/Puigdevall_P/midbrainDataset/rawData/41467_2025_67779_MOESM3_ESM.xlsx",
                             sheet="Data S2", skip=2)

qcTab$SRA_identifier <- gsub("cr_output_","",qcTab$analysis)
qcTab$Sample.Name <- metadata_sra[match(qcTab$SRA_identifier, metadata_sra$Run),]$Sample.Name
qcTab$library_10x <- link_ids[match(qcTab$Sample.Name, link_ids$V1),]$V2
qcTab$source_name <- metadata_sra[match(qcTab$SRA_identifier, metadata_sra$Run),]$source_name
qcTab$cell_line <- metadata_sra[match(qcTab$SRA_identifier, metadata_sra$Run),]$cell_line

donor_metadata <- read_excel("/shared/projects/tp_2616_fnom_183960/Puigdevall_P/midbrainDataset/rawData/41467_2025_67779_MOESM3_ESM.xlsx",
                             sheet="Data S14", skip=2)

qcTab$vcfId <- sapply(unique(qcTab$library_10x), function(y){
  paste0(unique(subset(donor_metadata, midBrainId==y)$vcfId),collapse=",")
}, simplify=T)

qcTab$condition <- sapply(unique(qcTab$library_10x), function(y){
  paste0(subset(donor_metadata, midBrainId==y)$condition,collapse=",")
}, simplify=T)

qcTab$timePoint <- sapply(unique(qcTab$library_10x), function(y){
  unique(subset(donor_metadata, midBrainId==y)$timePoint)
}, simplify=T)



list_10xfiles <- paste0(path_to_cellranger_output,
                        dir(path_to_cellranger_output),
                        "/outs/filtered_feature_bc_matrix")


## The sapply function below computes the percentage of singletons per 10x library. It assumes that demuxlet has been run previously.

perc <- sapply(list_10xfiles, function(x){
  
  demuxlet_file <- dir(gsub("/filtered_feature_bc_matrix","", x), full.names=T)[grepl("noFilter.best",dir(gsub("/filtered_feature_bc_matrix","", x)))]
  
  ## Load demuxlet file
  demuxlet <- read.table(demuxlet_file, header=T, sep="\t")
  singletons_perc <- unname(round(table(grepl("SNG-", demuxlet$BEST))*100/sum(table(grepl("SNG-", demuxlet$BEST))),2)["TRUE"])
  tmp <- data.frame(sample=sapply(strsplit(demuxlet_file, "/"), function(x) x[7]),
                    percSingletons=singletons_perc)
  tmp$sample <- gsub("cr_output_","",tmp$sample)
  return(tmp)
  
}, simplify=F)

perc <- do.call("rbind", perc)
rownames(perc) <- NULL

qcTab <- merge(qcTab, perc, by.x="SRA_identifier", by.y="sample")
qcTab$sampleIndex <- 1:length(qcTab$SRA_identifier)

###
###


## It saves a list of Seurat objects corresponding to each 10x library.

allobj <- sapply(list_10xfiles, function(x){
  
  tenXrun <- Read10X(data.dir = x)
  tenXrun <- CreateSeuratObject(counts = tenXrun, project = "MLO")
  sample <- gsub("cr_output_","",sapply(strsplit(x,"/"), function(x) x[7]))
  
  demuxlet_file <- dir(gsub("/filtered_feature_bc_matrix","", x), full.names=T)[grepl("noFilter.best",dir(gsub("/filtered_feature_bc_matrix","", x)))]
  tenXrun@meta.data$origin <-  "iPSC"
  
  demuxlet <- read.table(demuxlet_file, header=T, sep="\t")
  
  ## Filter-out all cells that have not been demultiplexed as an unambiguous singleton
  
  tenXrun <- tenXrun[,which(!is.na(match(rownames(tenXrun@meta.data), demuxlet$BARCODE)))]
  tenXrun <- tenXrun[,which(grepl("SNG-", demuxlet$BEST))]
  demuxlet <- demuxlet[grepl("SNG-", demuxlet$BEST),]
  
  stopifnot(dim(demuxlet)[1]==dim(tenXrun)[2])
  
  tenXrun@meta.data$vcfId <- gsub("SNG-","",demuxlet[match(rownames(tenXrun@meta.data), demuxlet$BARCODE),]$BEST)
  tenXrun@meta.data$SRA_identifier <- sample
  tenXrun@meta.data$midBrainId <- qcTab[match(tenXrun@meta.data$SRA_identifier, qcTab$SRA_identifier),]$library_10x
  
  donor_metadata <- read_excel("/shared/projects/tp_2616_fnom_183960/Puigdevall_P/midbrainDataset/rawData/41467_2025_67779_MOESM3_ESM.xlsx",
                               sheet="Data S14", skip=2)
  donor_metadata <- as.data.frame(donor_metadata)
  tenXrun$donorId <-donor_metadata[match(tenXrun$vcfId, donor_metadata$vcfId),]$donorId
  tenXrun$timePoint <- donor_metadata[match(tenXrun$midBrainId, donor_metadata$midBrainId),]$timePoint
  tenXrun$originDimensions <- donor_metadata[match(tenXrun$midBrainId, donor_metadata$midBrainId),]$originDimensions
  tenXrun@meta.data$donorDimensions <- paste0(tenXrun$donorId,"-",tenXrun$originDimensions)

  limsNFeatures <- c(1500,6000)
  limsMito <- 10
    
  
  idDef <- gsub("cr_output_","",unname(sapply(strsplit(x,"/"), function(y) y[7])))
  tenXrun[["percent.mt"]] <- PercentageFeatureSet(tenXrun, pattern = "^MT-")
  
  plot1 <- FeatureScatter(tenXrun, feature1 = "nCount_RNA", feature2 = "percent.mt")+
    ggtitle(idDef)+theme(plot.title=element_text(hjust=0.5, size=11),
                         legend.position="none")
  plot2 <- FeatureScatter(tenXrun, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
    ggtitle(idDef)+
    theme(plot.title=element_text(hjust=0.5, size=11),
          legend.position="none")
  qcplot1 <- plot1 + plot2
  plot(qcplot1)
  
  tmpTabQC <- data.frame(nFeature_RNA=sort(tenXrun$nFeature_RNA),
                         index=1:length(sort(tenXrun$nFeature_RNA)))
  
  qcplot2 <- ggplot(tmpTabQC, aes(x=index, y=nFeature_RNA))+
    geom_point(alpha=0.5, colour="black", size=3, shape=1)+
    theme_bw()+
    ggtitle(paste0("nFeature_RNA distribution: ", idDef))+
    theme(plot.title=element_text(hjust=0.5, size=12, face="bold"))+
    geom_hline(yintercept=limsNFeatures, linetype="dashed", col="red")+
    ylab("Number of expressed genes per cell")+
    xlab("Cells sorted by nFeature_RNA")
  
  plot(qcplot2)
  
  
  tmpTabQC <- data.frame(percent.mt=sort(tenXrun$percent.mt),
                         index=1:length(sort(tenXrun$percent.mt)))
  
  qcplot3<- ggplot(tmpTabQC, aes(x=index, y=percent.mt))+
    geom_point(alpha=0.5, colour="black", size=3, shape=1)+
    theme_bw()+
    ggtitle(paste0("Mitochondrial content: ", idDef))+
    theme(plot.title=element_text(hjust=0.5, size=14, face="bold"))+
    geom_hline(yintercept=limsMito, linetype="dashed", col="red")+
    ylab("% of Mitochondrial counts")+
    xlab("Cells sorted by mitochondrial content")
  
  plot(qcplot3)

  
  tenXrun <- subset(tenXrun, subset = nFeature_RNA > limsNFeatures[1] & nFeature_RNA < limsNFeatures[2] & percent.mt < limsMito)
  
  ## List of filtered-cells per 10x library ready to be merged
  tenXrun
  
})



names(allobj) <- gsub("cr_output_","",sapply(strsplit(names(allobj), "/"), function(x) x[7]))

mlo.merged.subset120 <- merge(allobj[[1]], y=do.call("c",allobj[2:length(allobj)]),
                              add.cell.ids = qcTab[match(names(allobj), qcTab$SRA_identifier),]$sampleIndex,
                              project="MLO")

saveRDS(mlo.merged.subset120, file="/shared/projects/tp_2616_fnom_183960/Puigdevall_P/midbrainDataset/merged_subset_mlo_day120.rds")


## In the next processing steps, we will take a closer look to the step of removing doublets & multiplets, as part of our post-Cell Ranger quality control. 
## We can revisit these demuxlet results and do a results cross-comparison with those tools later, so keep in mind this script for pulling demuxlet results and the 
## matrices with pre-QA cells.

## Further QC-tests (load Seurat objects)

test1 <- readRDS("/shared/projects/tp_2616_fnom_183960/Puigdevall_P/midbrainDataset/DataS3_mlo_resolution075_Annot.RDS")
dim(test1)
# [1] 22032 49284

test2 <- readRDS("/shared/projects/tp_2616_fnom_183960/Puigdevall_P/midbrainDataset/DataS4_midBrainIntegration.RDS")
dim(test2)
# [1]  34703 194521







