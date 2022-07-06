# description ----

# This is the R script for generating Supplemental Figure 15 in 
# "Urinary single-cell sequencing captures intrarenal injury and repair processes 
# in human acute kidney injury" by Klocke et al. 

# loading packages and data ----

## loading packages ----
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(scDblFinder)
library(harmony)
library(celldex)
library(patchwork)
library(RColorBrewer)
library(viridis)
library(paletteer) 
library(hrbrthemes)

## loading own data ----
URINE <- readRDS("~/SO_all_urine_cells.rds")

## loading repository data ----
raw.dir <- ""

# data from GSE157640, which included a pooled urine sample from 10 healthy volunteers (otherwise diabetic nephropathy patients)
H1.data <- Read10X(data.dir = file.path(raw.dir, "GSE157640_urine_healthy/SRR12607414/raw_feature_bc_matrix"))
H2.data <- Read10X(data.dir = file.path(raw.dir, "GSE157640_urine_healthy/SRR12607415/raw_feature_bc_matrix"))
H3.data <- Read10X(data.dir = file.path(raw.dir, "GSE157640_urine_healthy/SRR12607416/raw_feature_bc_matrix"))
H4.data <- Read10X(data.dir = file.path(raw.dir, "GSE157640_urine_healthy/SRR12607417/raw_feature_bc_matrix"))

# data from GSE165396, which included a pooled urine sample from 10 healthy volunteers (otherwise diabetic nephropathy patients)
H5.data <- Read10X(data.dir = file.path(raw.dir, "GSE165396_urine_healthy"))

                           
## loading color schemes ----
cellgroup_colors <- brewer.pal(5, "Paired")
names(cellgroup_colors) <- levels(factor(H_AKI$cellgroup))[c(4,3,2,1,5)]

## functions ----
number_ticks <- function(n) {function(limits) pretty(limits, n)}

# create and merge Seurat Objects ----

## create separate objects ----
H1 <- CreateSeuratObject(counts = H1.data, project = "H1", min.cells = 3, min.features = 200)
H2 <- CreateSeuratObject(counts = H2.data, project = "H2", min.cells = 3, min.features = 200)
H3 <- CreateSeuratObject(counts = H3.data, project = "H3", min.cells = 3, min.features = 200)
H4 <- CreateSeuratObject(counts = H4.data, project = "H4", min.cells = 3, min.features = 200)
H5 <- CreateSeuratObject(counts = H5.data, project = "H5", min.cells = 3, min.features = 200)


H_List <- list(H1, H2, H3, H4)

## determine percent mitochondrial RNA, doublets and merging of samples. ----
compileSO <- function(URINEList) { for (i in 1:length(URINEList)) {
  URINEList[[i]][["percent.mt"]] <- PercentageFeatureSet(URINEList[[i]], pattern = "^MT-")
  doublets <- scDblFinder(GetAssayData(URINEList[[i]], assay = "RNA", slot = "data"))
  doublets <- as.vector(doublets@colData@listData[["scDblFinder.class"]])
  URINEList[[i]]@meta.data$multiplet_class <- doublets
  URINEList[[i]]@project.name <- levels(URINEList[[i]]@active.ident)
  URINEList[[i]] <- RenameCells(URINEList[[i]], add.cell.id = paste0(URINEList[[i]]$orig.ident, "_"))
}
  
  # merge all the objects in the list
  URINE <- purrr::reduce(URINEList, merge)
  URINE <- SetIdent(URINE, value = "orig.ident")
  return(URINE)
}

H <- compileSO(H_List)

# QC filtering, normalization, dim. red. and clustering ----
H0 <- H

processSO <- function(URINE) {
  ## subset by percent.mt, singlets and features
  URINE <- subset(URINE, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 500 & nCount_RNA < 50000 & percent.mt < 20 & multiplet_class == "singlet")
  
  URINE <- SCTransform(object = URINE, vars.to.regress = "percent.mt", method = "glmGamPoi", verbose = T)
  URINE <- RunPCA(object = URINE, dims = 1:45, verbose = T)
  
  ElbowPlot(URINE, ndims = 45)
  
  #perform batch correction with harmony
  URINE <- RunHarmony(URINE, group.by.vars = c("orig.ident"), lambda = 2, tau = 1000, theta = 1, assay.use = "SCT")
  #run dimensionality reduction and find clusters
  URINE <- RunUMAP(object = URINE, reduction = "harmony", dims = 1:45, verbose = T)
  URINE <- RunTSNE(object = URINE, reduction = "harmony", dims = 1:45, verbose = T)
  URINE <- FindNeighbors(object = URINE, reduction = "harmony", dims = 1:45, verbose = T)
  URINE <- FindClusters(object = URINE, resolution = 0.3, verbose = T)
  #cell cycle scoring
  URINE <- CellCycleScoring(URINE, assay = 'SCT', s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
  
  return(URINE)
}

H <- processSO(H)

DimPlot(H)

# merging with AKI samples SO "URINE" ----

URINE$disease <- "AKI"
H$disease <- "healthy"

## downsampling URINE to cell count of healthy samples ----
URINE <- SetIdent(URINE, value="disease")
URINEsub <- subset(URINE, downsample = ncol(H))

DimPlot(URINEsub)

## merge AKI ("URINE") and healthy ("H") SOs ----
H_AKI <- merge(URINEsub, H)


H_AKI <- SCTransform(object = H_AKI, vars.to.regress = "percent.mt", method = "glmGamPoi", verbose = T)
H_AKI <- RunPCA(object = H_AKI, dims = 1:45, verbose = T)
#perform batch correction for disease type with harmony
H_AKI <- RunHarmony(H_AKI, group.by.vars = "disease", assay.use = "SCT")
#perform dimensionality reduction and find clusters
H_AKI <- RunUMAP(object = H_AKI, reduction = "harmony", dims = 1:45, verbose = T)
H_AKI <- RunTSNE(object = H_AKI, reduction = "harmony", dims = 1:45, verbose = T)
H_AKI <- FindNeighbors(object = H_AKI, reduction = "harmony", dims = 1:45, verbose = T)
H_AKI <- FindClusters(object = H_AKI, resolution = 0.2, verbose = T)

DimPlot(H_AKI, label=T)


H_AKI.markers <- FindAllMarkers(H_AKI)
markers <- H_AKI.markers %>% 
  group_by(cluster) %>% 
  top_n(2, avg_log2FC)
DotPlot(H_AKI, features = c(unique(markers$gene), "EPCAM", "CRYAB", "PSCA"))+ coord_flip()
FeaturePlot(H_AKI, )

H_AKIClusternames <- c("UGEC", "UGEC", "MYEL", "TEC", "TEC", "UGEC", "LYMP", "PDC")
# URINE$celltype includes names of clusters
names(H_AKIClusternames) <- levels(H_AKI)
H_AKI <- RenameIdents(H_AKI, H_AKIClusternames)
H_AKI$cellgroup <- H_AKI@active.ident
H_AKI <- SetIdent(H_AKI, value ="seurat_clusters")
H_AKI$cellgroup <- factor(H_AKI$cellgroup, levels(H_AKI$cellgroup)[c(3,5,2,4,1)])
levels(H_AKI$cellgroup)

# Suppl. Fig. 15 ----

a <-DimPlot(H_AKI, group.by = "cellgroup", split.by = "disease", label=F) + scale_color_manual(values=cellgroup_colors) + ggtitle("urinary cells")

b<- H_AKI@meta.data%>% group_by(disease, cellgroup) %>%
  summarise(count = n()) %>%
  spread(cellgroup, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('disease', 'total_cell_count', everything())) %>%
  arrange(factor(disease, levels = levels(H_AKI$disease))) %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'disease') %>%
  ggplot(aes(disease, value)) +
  geom_bar(aes(fill = variable), position = 'fill', stat = 'identity') +
  scale_fill_manual(values=cellgroup_colors) +
  scale_y_continuous(name = 'Celltype \n', breaks=number_ticks(4)) +
  coord_cartesian(clip = 'off') +
  theme_classic() +
  theme(
    legend.position = 'right',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(vjust = 0.5, hjust=0.5, angle=90, face= "bold", size=16),
    axis.title.x = element_blank(),
    axis.text.y = element_text(vjust = 0.5, hjust=0, face= "bold", size=16), 
    axis.title.y = element_blank(), 
    plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
  ) +NoLegend()

## compile and save ----
layout <- 
  "AAAAB"
pw <- a+b+plot_layout(design=layout, guides = "collect")

ggsave(plot = pw, "~/healthy.png", w = 4000, h = 2000, units = "px", scale = 1, type = "cairo-png")

