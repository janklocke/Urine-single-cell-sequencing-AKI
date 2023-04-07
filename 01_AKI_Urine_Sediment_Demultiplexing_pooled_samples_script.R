# description ----

# This is the R script for demultiplexing previously pooled and HTO-tagged urine samples from raw data found in
# "Urinary single-cell sequencing captures intrarenal injury and repair processes 
# in human acute kidney injury" by Klocke et al. 
# We will generate seurat objects with patient metadata which will be further analyzed in the following scripts.
# The script will mostly follow the workflow introduced in this seurat vignette: 
# https://satijalab.org/seurat/articles/hashing_vignette.html#8-hto-dataset-from-human-pbmcs

# loading packages and data ----

## loading packages ----
library(Seurat)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(paletteer) 
library(patchwork)

## loading data ----
Pool1neg.data <- Read10X(data.dir = "~/UR_AKI_Pool1_cd45-/filtered_feature_bc_matrix")
Pool1pos.data <- Read10X(data.dir = "~/UR_AKI_Pool1_cd45+/filtered_feature_bc_matrix")
Pool2neg.data <- Read10X(data.dir = "~/UR_AKI_Pool2_cd45-/filtered_feature_bc_matrix")
Pool2pos.data <- Read10X(data.dir = "~/UR_AKI_Pool2_cd45+/filtered_feature_bc_matrix")
Pool3.data <- Read10X(data.dir = "~/UR_AKI_Pool3/filtered_feature_bc_matrix")
Pool4.data <- Read10X(data.dir = "~/UR_AKI_Pool4/filtered_feature_bc_matrix")
Pool5.data <- Read10X(data.dir = "~/UR_AKI_Pool5/filtered_feature_bc_matrix")
Pool6.data <- Read10X(data.dir = "~/UR_AKI_Pool6/filtered_feature_bc_matrix")

# Preparing the data ----

## seperating UMI and HTO libraries ----
Pool1neg.htos <- Pool1neg.data[["Antibody Capture"]]
Pool1neg.umis <- Pool1neg.data[["Gene Expression"]]
Pool1pos.htos <- Pool1pos.data[["Antibody Capture"]]
Pool1pos.umis <- Pool1pos.data[["Gene Expression"]]
Pool2neg.htos <- Pool2neg.data[["Antibody Capture"]]
Pool2neg.umis <- Pool2neg.data[["Gene Expression"]]
Pool2pos.htos <- Pool2pos.data[["Antibody Capture"]]
Pool2pos.umis <- Pool2pos.data[["Gene Expression"]]
Pool3.htos <- Pool3.data[["Antibody Capture"]]
Pool3.umis <- Pool3.data[["Gene Expression"]]
Pool4.htos <- Pool4.data[["Antibody Capture"]]
Pool4.umis <- Pool4.data[["Gene Expression"]]
Pool5.htos <- Pool5.data[["Antibody Capture"]]
Pool5.umis <- Pool5.data[["Gene Expression"]]
Pool6.htos <- Pool6.data[["Antibody Capture"]]
Pool6.umis <- Pool6.data[["Gene Expression"]]

PoolHtoList <- c(Pool1neg.htos, Pool1pos.htos, Pool2neg.htos, Pool2pos.htos, Pool3.htos, Pool4.htos, Pool5.htos, Pool6.htos)
PoolUmiList <- c(Pool1neg.umis, Pool1pos.umis, Pool2neg.umis, Pool2pos.umis, Pool3.umis, Pool4.umis, Pool5.umis, Pool6.umis)

##  finding cells that are tagged with at least one AB and subset to those cells.
for (i in 1:length(PoolHtoList)) {
  Pool.htos2 <- as.matrix(PoolHtoList[[i]])
  n <- c(rownames(Pool.htos2))
  Pool.htos2 <- as.data.frame(t(Pool.htos2[,-1]))
  colnames(Pool.htos2) <- n
  Pool.htos2$cells <- factor(row.names(Pool.htos2))
  # Check the column types
  str(Pool.htos2) 
  # filter for cells with >0 AB
  Pool.htos2 <- Pool.htos2 %>%
    mutate(TotalAB = CD298B0251 + CD298B0252 + CD298B0253) %>%
    filter(TotalAB > 0) 
  cell.list <- Pool.htos2$cells
  # Subset RNA and HTO counts by cell barcodes that were AB tagged
  PoolUmiList[[i]] <- PoolUmiList[[i]][, cell.list]
  PoolHtoList[[i]] <- as.matrix(PoolHtoList[[i]][, cell.list])
  # Confirm that the HTO have the correct names
  rownames(PoolHtoList[[i]])
}

## create seurat objects from UMI libraries ----
Pool1neg <- CreateSeuratObject(counts = PoolUmiList[[1]], project = "Pool1")
Pool1pos <- CreateSeuratObject(counts = PoolUmiList[[2]], project = "Pool1")
Pool2neg <- CreateSeuratObject(counts = PoolUmiList[[3]], project = "Pool2")
Pool2pos <- CreateSeuratObject(counts = PoolUmiList[[4]], project = "Pool2")
Pool3 <- CreateSeuratObject(counts = PoolUmiList[[5]], project = "Pool3")
Pool4 <- CreateSeuratObject(counts = PoolUmiList[[6]], project = "Pool4")
Pool5 <- CreateSeuratObject(counts = PoolUmiList[[7]], project = "Pool5")
Pool6 <- CreateSeuratObject(counts = PoolUmiList[[8]], project = "Pool6")


## add HTO libraries as a new assay to Seurat objects ----
PoolList <- c(Pool1neg, Pool1pos, Pool2neg, Pool2pos, Pool3, Pool4, Pool5, Pool6)

for (i in 1:length(PoolList)) {
  # Add HTO data as a new assay independent from RNA
  PoolList[[i]][["HTO"]] <- CreateAssayObject(counts = PoolHtoList[[i]])
  # Normalize HTO data, here we use centered log-ratio (CLR) transformation
  PoolList[[i]] <- NormalizeData(PoolList[[i]], assay = "HTO", normalization.method = "CLR")
}

# Demultiplexing ----

## demultiplex using HTODemux function ----
## using low quantile because of noisy data in some cases

Pool1neg <- MULTIseqDemux(PooList[[1]], quantile = 0.20) # HTODemux does not function in this case, MULTIseqDemux was applied. 
Pool1pos <- HTODemux(PoolList[[2]], positive.quantile = 0.95)
Pool2neg <- HTODemux(PoolList[[3]], positive.quantile = 0.90)
Pool2pos <- HTODemux(PoolList[[4]], positive.quantile = 0.99)
Pool3 <- HTODemux(PoolList[[5]], positive.quantile = 0.99)
Pool4 <- HTODemux(PoolList[[6]], positive.quantile = 0.90)
Pool5 <- HTODemux(PoolList[[7]], positive.quantile = 0.90)
Pool6 <- HTODemux(PoolList[[8]], positive.quantile = 0.99)

Pool1neg$hash.ID <-Pool1neg$MULTI_ID

## show the Demux results in table ---- 
table((Pool3$hash.ID))
table((Pool4$hash.ID))
table((Pool5$hash.ID))
table((Pool6$hash.ID))

Pool1 <- merge(Pool1neg, Pool1pos)
Pool2 <- merge(Pool2neg, Pool2pos)

PoolList <- c(Pool1, Pool2, Pool3, Pool4, Pool5, Pool6)

# Suppl. Fig. 3, Visualization of Demux results ----

## preparation of SOs ----
# First, we will remove negative cells from the object
plot.list <- list(1,2,3,4,5,6)
for (i in 1:length(PoolList)) {
Pool.subset <- subset(PoolList[[i]], idents = "Negative", invert = TRUE)
# Calculate a distance matrix using HTO
hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = Pool.subset, assay = "HTO"))))
# Calculate tSNE embeddings with a distance matrix
Pool.subset <- RunTSNE(Pool.subset, distance.matrix = hto.dist.mtx, perplexity = 100)
plot.list[[i]] <- DimPlot(Pool.subset, pt.size= 2) 
}

## color schemes ----
Demuxcolors <- c(paletteer_d("ggsci::red_material")[5], paletteer_d("palettesForR::Pastels")[c(3,5,6,9,13)])
names(Demuxcolors) <- levels(Pool3)

## Suppl. Fig 3A, TSNE ----
for (i in 1:4) { plot.list[[i]] <- plot.list[[i]]+ 
  theme_minimal()+
  scale_color_manual(values = Demuxcolors)+
  ggtitle(paste("sample pool", i))}
for (i in 2:4) { plot.list[[i]] <- plot.list[[i]]+ 
  NoLegend()}

plot.list[[3]] <- plot.list[[3]]+ggtitle("sample pool 1")
plot.list[[5]] <- plot.list[[5]]+ggtitle("sample pool 2")
plot.tsne <- plot.list[[3]] + plot.list[[5]]+ plot_layout(guides= "collect") & theme(title = element_text(size=24, face="bold", hjust=0.5),
                                                                                      text = element_text(size=24),
                                                                                      legend.text = element_text(size=24))

png("plottsne.png", width=1600, height=700, type="cairo-png")
plot.tsne
dev.off()

## Suppl. Fig 3C, heatmaps ----
plot.heatmap <- HTOHeatmap(Pool3, assay = "HTO")+xlab("sample pool 1")+theme(axis.text.y=element_text(size=24, face="bold"))+
  HTOHeatmap(Pool5, assay = "HTO")+xlab("sample pool 2")+theme(axis.text.y=element_blank())+
 plot_layout(guides="collect") & theme(axis.title.x = element_text(size=24, face="bold"), legend.text = element_text(size=24),legend.title = element_text(size=24))

png("plotheatmap.png", width=1600, height=400, type="cairo-png")
plot.heatmap
dev.off()

## Suppl. Fig 3B, ridge plots ----
a <-   RidgePlot(Pool3, assay = "HTO", features = "CD298B0251")&  
  scale_fill_manual(values = Demuxcolors)&
  scale_y_discrete(limits= c("Negative", "Doublet", "CD298B0254", "CD298B0253", "CD298B0252", "CD298B0251"))&
  theme_minimal()&
  theme(axis.title = element_blank(),
        text = element_text(size=24))&
  NoLegend()
b <-   RidgePlot(Pool3, assay = "HTO", features = "CD298B0252")&  
  scale_fill_manual(values = Demuxcolors)&
  scale_y_discrete(limits= c("Negative", "Doublet", "CD298B0254", "CD298B0253", "CD298B0252", "CD298B0251"))&
  theme_minimal()&
  theme(axis.title = element_blank(),
        text = element_text(size=24),
        axis.text.y = element_blank())&
  NoLegend()
c <-   RidgePlot(Pool3, assay = "HTO", features = "CD298B0253")&  
  scale_fill_manual(values = Demuxcolors)&
  scale_y_discrete(limits= c("Negative", "Doublet", "CD298B0254", "CD298B0253", "CD298B0252", "CD298B0251"))&
  theme_minimal()&
  theme(axis.title = element_blank(),
        text = element_text(size=24),
        axis.text.y = element_blank())&
  NoLegend()
d <-   RidgePlot(Pool3, assay = "HTO", features = "CD298B0254")&  
  scale_fill_manual(values = Demuxcolors)&
  scale_y_discrete(limits= c("Negative", "Doublet", "CD298B0254", "CD298B0253", "CD298B0252", "CD298B0251"))&
  theme_minimal()&
  theme(axis.title = element_blank(),
        text = element_text(size=24),
        axis.text.y = element_blank())&
  NoLegend()

layout <- "ABCD"
plot.ridge <- a+b+c+d+plot_layout(design=layout)

png("plotridge.png", width=1600, height=400, type="cairo-png")
plot.ridge
dev.off()

# Add sample names, remove doublets, save ----

## remove doublets ----
for (i in 1:length(PoolList)) {
  PoolList[[i]] <- SetIdent(PoolList[[i]], "hash.ID")
  PoolList[[i]] <- subset(PoolList[[i]], ident = "Doublet", invert = T)
  }

## load barcode patient list ----
barcodes <- read.delim("~/Urine_AKI_barcode_list.csv", sep=";")
rownames(barcodes) <- barcodes[,1]
barcodes <- barcodes[,2:5]


## assign sample names to cells via barcode ----
for (i in 1:6) {
  PoolList[[i]]$orig.ident <- NA
  PoolList[[i]]@meta.data[PoolList[[i]]$hash.ID == "CD298B0251", ]$orig.ident <- barcodes[rownames(barcodes)[i], "CD298B0251"]
  PoolList[[i]]@meta.data[PoolList[[i]]$hash.ID == "CD298B0252", ]$orig.ident <- barcodes[rownames(barcodes)[i], "CD298B0252"]
  PoolList[[i]]@meta.data[PoolList[[i]]$hash.ID == "CD298B0253", ]$orig.ident <- barcodes[rownames(barcodes)[i], "CD298B0253"]
}
for (i in c(1,2,3,5,6)) {
  PoolList[[i]]@meta.data[PoolList[[i]]$hash.ID == "CD298B0254", ]$orig.ident <- barcodes[rownames(barcodes)[i], "CD298B0254"]
}

PoolList[[1]]@meta.data[is.na(PoolList[[1]]$orig.ident),]$orig.ident <- "pool1"
PoolList[[2]]@meta.data[is.na(PoolList[[2]]$orig.ident),]$orig.ident <- "pool2"
PoolList[[3]]@meta.data[is.na(PoolList[[3]]$orig.ident),]$orig.ident <- "pool3"
PoolList[[4]]@meta.data[is.na(PoolList[[4]]$orig.ident),]$orig.ident <- "pool4"
PoolList[[5]]@meta.data[is.na(PoolList[[5]]$orig.ident),]$orig.ident <- "pool5"
PoolList[[6]]@meta.data[is.na(PoolList[[6]]$orig.ident),]$orig.ident <- "pool6"

## save pool Seurat objects ----
PoolNames <- c("AKI_Urine_sediment_Pool1.rds", "AKI_Urine_sediment_Pool3.rds","AKI_Urine_sediment_Pool2.rds", "AKI_Urine_sediment_Pool4.rds", "AKI_Urine_sediment_Pool5.rds", "AKI_Urine_sediment_Pool6.rds")

for (i in 1:length(PoolList)) {
  saveRDS(PoolList[[i]], file = PoolNames[[i]])
}


