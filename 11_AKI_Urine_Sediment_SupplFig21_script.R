# description ----

# This is the R script for generating Supplemental Figure 21 in 
# "Urinary single-cell sequencing captures intrarenal injury and repair processes 
# in human acute kidney injury" by Klocke et al. 

# loading packages and data ----

## loading packages ----
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(celldex)
library(harmony)
library(patchwork)
library(celldex)
library(RColorBrewer)
library(viridis)
library(edgeR)
library(corrplot)
library("Hmisc")
library(lattice)
library("PerformanceAnalytics")
library(ggcorrplot2)

## loading data ----
URINE <- readRDS("~/SO_all_urine_cells.rds")
RENAL <-readRDS("~/SO_kidney_urine_cells.rds")

## Healthy tumor adjacent kidney tissue data was downloaded from the Kidney Precision Medicine Project (KPMP) database 
## (participants 17-1606, 18-139, 18-162, 18-342, 18142-5).

k1.data <- Read10X("~/Forschung/Rohdaten/NiereScRNAseq/17-1606-2-0")
k2.data <- Read10X("~/Forschung/Rohdaten/NiereScRNAseq/17-1606-2-1")
k3.data <- Read10X("~/Forschung/Rohdaten/NiereScRNAseq/18-139")
k4.data <- Read10X("~/Forschung/Rohdaten/NiereScRNAseq/18-162")
k5.data <- Read10X("~/Forschung/Rohdaten/NiereScRNAseq/18-342")
k6.data <- Read10X("~/Forschung/Rohdaten/NiereScRNAseq/18142-5-1")
k7.data <- Read10X("~/Forschung/Rohdaten/NiereScRNAseq/18142-5-2")

## scRNAseq data of 3 human bladder samples (GSE129845) was downloaded 
## from the National Center for the Biotechnology Information Database (NCBI) / Gene expression omnibus (GEO) repository.

u1.data <- Read10X("~/Forschung/Rohdaten/HarnblaseScRNAseq/1")
u2.data <- Read10X("~/Forschung/Rohdaten/HarnblaseScRNAseq/2")
u3.data <- Read10X("~/Forschung/Rohdaten/HarnblaseScRNAseq/3")


## load color schemes ----
renalcelltype_colors <- 
  c(paletteer_d("ggsci::purple_material")[c(3,5)],                  #PDC colors
    paletteer_d("ggsci::teal_material")[c(8,6,4,3)],                   #healthy TEC PT and LOH colors   
    paletteer_d("ggsci::light_blue_material")[c(3,5,7,9)],          #healthy TEC DT and CD colors
    paletteer_d("ggsci::amber_material")[c(3,4,6,8,10)],            #EMT/infl TEC
    paletteer_d("ggsci::red_material")[c(2,5,8)],                   #OxyStress TEC
    paletteer_d("ggsci::pink_material")[4],                         #prolif TEC
    paletteer_d("ggsci::light_green_material")[c(9,7,5,3)])         #progenitors 
names(renalcelltype_colors) <- levels(RENAL$renalcelltype_nr)
celltype_colors <- colorRampPalette(paletteer_d("palettesForR::Pastels"))(19)
HKcolors <- paletteer_d("palettesForR::Pastels")[c(2:8,10:15)]
JOINED_Ucolors <- c(renalcelltype_colors[c(1,3:10,12,16,19,22)], celltype_colors[c(9,15,16, 18)]) 
UROcolors <- paletteer_d("dichromat::Categorical_12")[c(1:3,5:9,11:12)]

names(JOINED_Ucolors) <- levels(JOINED$rcg_nr)[1:17]
names(UROcolors) <- levels(JOINED$UROcelltype)[1:10]
names(HKcolors) <- levels(JOINED$HKcelltype)[1:13]

names(UROcolors) <- levels(JOINED$URO_nr)[1:10]
names(HKcolors) <- levels(JOINED$HK_nr)[1:13]

## functions ----
makeTransparent <- function(someColor, alpha=100){
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}


# create and merge Seurat Objects ----

## Initialize Seurat objects ----
### healthy kidney ----
k1 <- CreateSeuratObject(counts = k1.data, project = "17-1606-2-0", min.cells = 3, min.features = 200)
k2 <- CreateSeuratObject(counts = k2.data, project = "17-1606-2-1", min.cells = 3, min.features = 200)
k3 <- CreateSeuratObject(counts = k3.data, project = "18-139", min.cells = 3, min.features = 200)
k4 <- CreateSeuratObject(counts = k4.data, project = "18-162", min.cells = 3, min.features = 200)
k5 <- CreateSeuratObject(counts = k5.data, project = "18-342", min.cells = 3, min.features = 200)
k6 <- CreateSeuratObject(counts = k6.data, project = "18142-5-1", min.cells = 3, min.features = 200)
k7 <- CreateSeuratObject(counts = k7.data, project = "18142-5-2", min.cells = 3, min.features = 200)
### healthy bladder ----
u1 <- CreateSeuratObject(counts = u1.data, project = "uro-s1", min.cells = 3, min.features = 200)
u2 <- CreateSeuratObject(counts = u2.data, project = "uro-s2", min.cells = 3, min.features = 200)
u3 <- CreateSeuratObject(counts = u3.data, project = "uro-s3", min.cells = 3, min.features = 200)

## make a list of all objects, determine percentage of mitochondrial RNA and doublets per sample. ----
HKList <- list(k1, k2, k3, k4, k5, k6, k7)
for (i in 1:length(HKList)) {
  HKList[[i]][["percent.mt"]] <- PercentageFeatureSet(HKList[[i]], pattern = "^MT-")
  HKList[[i]] <- RenameCells(HKList[[i]], add.cell.id = paste0(HKList[[i]]$orig.ident, "_"))
}

UROList <- list(u1, u2, u3)
for (i in 1:length(UROList)) {
  UROList[[i]][["percent.mt"]] <- PercentageFeatureSet(UROList[[i]], pattern = "^MT-")
  UROList[[i]] <- RenameCells(UROList[[i]], add.cell.id = paste0(UROList[[i]]$orig.ident, "_"))
}
## merge all the objects in the list ----
HK <- purrr::reduce(HKList, merge)
URO <- purrr::reduce(UROList, merge)


# QC analysis, subsetting and SCTransform, PCA ----

a <- VlnPlot(HK, c("nFeature_RNA"), group.by = "orig.ident", pt.size = 0.2) + 
  NoLegend() 
b <- VlnPlot(HK, c("nCount_RNA"), group.by = "orig.ident", pt.size = 0.2) + 
  NoLegend() 
c <- VlnPlot(HK, c("percent.mt"), group.by = "orig.ident", pt.size = 0.2) + 
  NoLegend()   
a / b / c

HK0 <- HK
HK <- subset(HK, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
HK <- SCTransform(object = HK, method = "glmGamPoi", verbose = T)
HK <- RunPCA(object = HK, dims = 1:15, verbose = T)
DimHeatmap(HK, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(HK, ndims = 15)

URO0 <- URO
URO <- subset(URO, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
URO <- SCTransform(object = URO, method = "glmGamPoi", verbose = T)
URO <- RunPCA(object = URO, dims = 1:15, verbose = T)
DimHeatmap(URO, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(URO, ndims = 15)

# batch correction with harmony ----
HK <- RunHarmony(HK, group.by.vars = c("orig.ident"), assay.use = "SCT")
URO <- RunHarmony(URO, group.by.vars = c("orig.ident"), assay.use = "SCT")

# dimensionality reduction and find clusters ----
HK <- RunUMAP(object = HK, reduction = "harmony", dims = 1:15, verbose = T)
HK <- RunTSNE(object = HK, reduction = "harmony", dims = 1:15, verbose = T)
HK <- FindNeighbors(object = HK, reduction = "harmony", dims = 1:15, verbose = T)
HK <- FindClusters(object = HK, resolution = 0.7, verbose = T)

URO <- RunUMAP(object = URO, reduction = "harmony", dims = 1:15, verbose = T)
URO <- RunTSNE(object = URO, reduction = "harmony", dims = 1:15, verbose = T)
URO <- FindNeighbors(object = URO, reduction = "harmony", dims = 1:15, verbose = T)
URO <- FindClusters(object = URO, resolution = 0.3, verbose = T)


DimPlot(HK, reduction = "umap", label = T, cols = custom_colors[[1]]) + NoLegend()
DimPlot(HK, reduction = "umap", label = T, group.by= "HKcelltype") 
HK.markers <- FindAllMarkers(HK, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
HK.top20 <- HK.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
HK.top3 <- HK.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)

DimPlot(URO, reduction = "umap", label = T) + NoLegend()
DimPlot(URO, reduction = "umap", label = F, group.by= "orig.ident") 

# cluster annotation ----

## determine marker genes in HK ----
HK.markers <- FindAllMarkers(HK, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
HK.top3 <- HK.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
## determine marker genes in URO ----
URO.markers <- FindAllMarkers(URO, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
URO.top3 <- URO.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)


## annotate HK clusters ----
HKClusterNames <- c("k_PT", "k_LYMP", "k_CD_PC", "k_DTL", 
                    "k_PT", "k_EC", "k_CD_IC", "k_PT", 
                    "k_EC", "k_TAL", "k_DCT", "k_CNT", 
                    "k_LYMP", "k_SMC", "k_MYEL", "k_PT", 
                    "k_ATL", "k_MYEL", "k_PDC", "k_EC")
names(HKClusterNames) <- levels(HK)
HK <- RenameIdents(HK, HKClusterNames)
HK$HKcelltype <- HK@active.ident
HK <- SetIdent(HK, value ="seurat_clusters")
HK$HKcelltype <- factor(HK$HKcelltype, levels(HK$HKcelltype)[c(13,1,4,12,7,8,9,3,6,5,10,11,2)])

HK <- SetIdent(HK, value ="HKcelltype")
levels(HK)
JOINEDClusterNrs  <- paste(1:13)
names(JOINEDClusterNrs) <- levels(HK)
HK <- RenameIdents(HK, JOINEDClusterNrs)
HK$HK_nr <- HK@active.ident
HK <- SetIdent(HK, value ="seurat_clusters")
levels(HK$HK_nr)

## annotate URO clusters ----
UROClusterNames <- c("b_intC", "b_basC", "b_FBR", "b_umbC", 
                     "b_MFB", "b_SMC", "b_basC2", "b_EC", 
                     "b_MYEL", "b_LYMP")
names(UROClusterNames) <- levels(URO)
URO <- RenameIdents(URO, UROClusterNames)
URO$UROcelltype <- URO@active.ident
URO <- SetIdent(URO, value ="seurat_clusters")
URO$UROcelltype <- factor(URO$UROcelltype, levels(URO$UROcelltype)[c(2,7,1,4,8,3,5,6,9,10)])

URO <- SetIdent(URO, value ="UROcelltype")
levels(URO)
JOINEDClusterNrs  <- paste(1:10)
names(JOINEDClusterNrs) <- levels(URO)
URO <- RenameIdents(URO, JOINEDClusterNrs)
URO$URO_nr <- URO@active.ident
URO <- SetIdent(URO, value ="seurat_clusters")
levels(URO$URO_nr)

# individual UMAP plots for healthy kdieny (HK) and bladder (URO) ----

d <- DimPlot(HK, reduction = "umap", group.by = "HK_nr", label = T, repel = T, pt.size = 1, label.size = 12) + 
  theme_void() + 
  scale_color_manual(labels = paste(levels(HK$HK_nr), levels(HK$HKcelltype)), values = HKcolors) + 
  theme() +
  ggtitle("kidney") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 32),
        plot.title = element_text(face= "bold", size = 48, hjust=0.5))

e <- DimPlot(URO, reduction = "umap", group.by= "URO_nr", label = T, repel = T, pt.size = 1, label.size = 12) + 
  theme_void() + 
  scale_color_manual(labels = paste(levels(URO$URO_nr), levels(URO$UROcelltype)), values = UROcolors) + 
  theme() +
  ggtitle("bladder") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 32),
        plot.title = element_text(face= "bold", size = 48, hjust=0.5))

png("supplfigindiv.png", width=2000, height= 1200, type="cairo-png") 
d+e
dev.off()


# merge kideny, bladder and urine datasets ----

## assign tissue information ----
URINE$tissue <- "urine"
HK$tissue <- "kidney"
URO$tissue <- "bladder"

## subset URINE data (exclude RBC) ----
URINE <- SetIdent(URINE, value = "renalcellgroup")
URINE_subset <- subset(URINE, idents = c("ERY"), invert = T)


## downsample URINE to same size as HK and URO ----
URINE_subset <- SetIdent(URINE_subset, value = "tissue")
URINE_subset <- subset(URINE_subset, downsample = (nrow(URO@meta.data) + nrow(HK@meta.data)))
table(URINE_subset$tissue)


## merge, sctransform and PCA ----
JOINED <- purrr::reduce(list(URINE_subset, URO, HK), merge)
JOINED <- SCTransform(object = JOINED, method = "glmGamPoi", verbose = T)
JOINED <- RunPCA(object = JOINED, dims = 1:30, verbose = T)
DimHeatmap(JOINED, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(JOINED, ndims = 40)

## perform batch correction for tisue type with harmony ----
JOINED <- RunHarmony(JOINED, group.by.vars = c("tissue"), assay.use = "SCT")
## run dimensionality reduction and find clusters ----
JOINED <- RunUMAP(object = JOINED, reduction = "harmony", dims = 1:30, verbose = T)
JOINED <- RunTSNE(object = JOINED, reduction = "harmony", dims = 1:30, verbose = T)
JOINED <- FindNeighbors(object = JOINED, reduction = "harmony", dims = 1:30, verbose = T)
JOINED <- FindClusters(object = JOINED, resolution = 0.7, verbose = T)

# visualization ----

JOINED <- SetIdent(JOINED, value = "tissue")

## assign other tissue identities to respective cell type metadata ----
JOINED@meta.data[JOINED$tissue == "kidney",]$renalcellgroup <- "kidney"
JOINED@meta.data[JOINED$tissue == "kidney",]$UROcelltype <- "kidney"
JOINED@meta.data[JOINED$tissue == "urine",]$UROcelltype <- "urine"
JOINED@meta.data[JOINED$tissue == "urine",]$HKcelltype <- "urine"
JOINED@meta.data[JOINED$tissue == "bladder",]$renalcellgroup <- "bladder"
JOINED@meta.data[JOINED$tissue == "bladder",]$HKcelltype <- "bladder"

## reorder clusters ----
JOINED$HKcelltype <- factor(JOINED$HKcelltype)
JOINED$UROcelltype <- factor(JOINED$UROcelltype)
JOINED$renalcellgroup <- factor(JOINED$renalcellgroup)
JOINED$HKcelltype <- factor(JOINED$HKcelltype, levels(JOINED$HKcelltype)[c(11,12,7,2,14,6,5,4,3,8,13,10,9,1,15)])
JOINED$UROcelltype <- factor(JOINED$UROcelltype, levels(JOINED$UROcelltype)[c(1,2,5,10,3,4,7,9,8,6,11,12)])
JOINED$renalcellgroup <- factor(JOINED$renalcellgroup, levels(factor(JOINED$renalcellgroup))[c(12,11,8,1,13,7,6,5,4,15,18,17,16,10,14,2,19,9,3)])

## add cluster numbering ----
JOINED <- SetIdent(JOINED, value ="renalcellgroup")
levels(JOINED)
JOINEDClusterNrs  <- paste(1:19)
names(JOINEDClusterNrs) <- levels(JOINED)
JOINED <- RenameIdents(JOINED, JOINEDClusterNrs)
JOINED$rcg_nr <- JOINED@active.ident
JOINED <- SetIdent(JOINED, value ="seurat_clusters")
levels(JOINED$rcg_nr)

JOINED <- SetIdent(JOINED, value ="UROcelltype")
levels(JOINED)
JOINEDClusterNrs  <- paste(1:12)
names(JOINEDClusterNrs) <- levels(JOINED)
JOINED <- RenameIdents(JOINED, JOINEDClusterNrs)
JOINED$URO_nr <- JOINED@active.ident
JOINED <- SetIdent(JOINED, value ="seurat_clusters")
levels(JOINED$URO_nr)

JOINED <- SetIdent(JOINED, value ="HKcelltype")
levels(JOINED)
JOINEDClusterNrs  <- paste(1:15)
names(JOINEDClusterNrs) <- levels(JOINED)
JOINED <- RenameIdents(JOINED, JOINEDClusterNrs)
JOINED$HK_nr <- JOINED@active.ident
JOINED <- SetIdent(JOINED, value ="seurat_clusters")
levels(JOINED$HK_nr)


HKcolors <- c(HKcolors, makeTransparent("gold", 8), makeTransparent("darkturquoise", 8))
names(HKcolors) <- levels(JOINED$HK_nr)
UROcolors <- c(UROcolors, makeTransparent("deeppink2", 8), makeTransparent("darkturquoise", 8))
names(UROcolors) <- levels(JOINED$URO_nr)
JOINED_Ucolors <- c(JOINED_Ucolors, makeTransparent("deeppink2", 8), makeTransparent("gold", 8))
names(JOINED_Ucolors) <- levels(JOINED$rcg_nr)

## Suppl. Fig 21 A,C,D, create plots ---- 
a <- DimPlot(JOINED, reduction = "umap", group.by = "HK_nr", label = T, repel = T, pt.size = 1, label.size = 12) + 
  theme_void() + 
  scale_color_manual(labels = paste(levels(JOINED$HK_nr), levels(JOINED$HKcelltype)), values = HKcolors) + 
  theme() +
  ggtitle("kidney") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 32),
        plot.title = element_text(face= "bold", size = 48, hjust=0.5))

b <- DimPlot(JOINED, reduction = "umap", group.by= "URO_nr", label = T, repel = T, pt.size = 1, label.size = 12) + 
  theme_void() + 
  scale_color_manual(labels = paste(levels(JOINED$URO_nr), levels(JOINED$UROcelltype)), values = UROcolors) + 
  theme() +
  ggtitle("bladder") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 32),
        plot.title = element_text(face= "bold", size = 48, hjust=0.5))

rcg_names <- str_replace(levels(JOINED$renalcellgroup), "DCT", "TAL/DCT")
rcg_names <- str_replace(rcg_names, "DTL", "TEC_VCAM1")
rcg_names <- str_replace(rcg_names, "TEC_emt", "TEC_inj")
c <- DimPlot(JOINED, reduction = "umap", group.by= "rcg_nr", label = T, repel = T, pt.size = 1, label.size = 12) + 
  theme_void() + 
  scale_color_manual(labels = paste(levels(JOINED$rcg_nr), rcg_names), values = JOINED_Ucolors) + 
  theme() +
  ggtitle("AKI urine sediment") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 32),
        plot.title = element_text(face= "bold", size = 48, hjust=0.5))

#### compile and save ----
png("supplfig21acd.png", width=3000, height= 1200, type="cairo-png") 
a+b+c
dev.off()

# Suppl. Fig. 21C, gene expression correlation plot ----

## function for Suppl Fig. 13C ----
# corr_df() generates a dataframe containing normalized mean expressions of all genes across clusters/cell types
corr_df <- function(SO, datacolumn) {
  KbyClusters.df <- data.frame(stringsAsFactors = FALSE, row.names = SO@assays$RNA@data@Dimnames[[1]])
  
  for(i in levels(factor(datacolumn))){
    print(i)
    ref_I <- subset(SO, cells = rownames(SO@meta.data[datacolumn %in% i, ]))
    nrow(SO@meta.data)
    ref_I.df <- as.data.frame(Matrix::rowMeans(ref_I@assays$RNA@counts))
    cpm <- apply(ref_I.df, 2, function(x) (x/sum(x))*1000000)
    KbyClusters.df <- cbind(KbyClusters.df, cpm)
  }
  
  # 1. remove if non-expressed genes
  KbyClusters.df <- KbyClusters.df[apply(KbyClusters.df == 0, 1, sum) != 9, ]
  
  # 2. rename column name
  colnames(KbyClusters.df) <- levels(factor(datacolumn))
  
  # 3. add pseudocount 1
  KbyClusters.df <- KbyClusters.df + 1
  
  # 4. log2 
  KbyClusters_log2.df <- as.data.frame(lapply(KbyClusters.df, log2), row.names = rownames(KbyClusters.df))
  return(KbyClusters_log2.df)
}
## calculating correlation matrix ----

# calculate mean gene expression per cell type
HK_corr <- corr_df(HK, HK$HKcelltype)
URO_corr <- corr_df(URO, URO$UROcelltype)
URINE_corr <- corr_df(URINE_subset, URINE_subset$renalcellgroup)


#filter for only the highly variable genes in the datasets
HK_corr2 <- HK_corr[rownames(HK_corr) %in% c(HK@assays[["SCT"]]@var.features[1:200], URO@assays[["SCT"]]@var.features[1:200], URINE@assays[["SCT"]]@var.features[1:200]),]
URO_corr2 <- URO_corr[rownames(URO_corr) %in% c(HK@assays[["SCT"]]@var.features[1:200], URO@assays[["SCT"]]@var.features[1:200], URINE@assays[["SCT"]]@var.features[1:200]),]
URINE_corr2 <- URINE_corr[rownames(URINE_corr) %in% c(HK@assays[["SCT"]]@var.features[1:200], URO@assays[["SCT"]]@var.features[1:200], URINE@assays[["SCT"]]@var.features[1:200]),]


# merge dataframes
MERGED_corr <- merge(HK_corr2, URO_corr2, by = 0, all = F)
rownames(MERGED_corr) <- MERGED_corr$Row.names
MERGED_corr$Row.names <- NULL
MERGED_corr <- merge(MERGED_corr, URINE_corr2, by = 0, all = F)
rownames(MERGED_corr) <- MERGED_corr$Row.names
MERGED_corr$Row.names <- NULL

# use rcorr() to calculate correlation matrix
MERGED.mtx <- rcorr(as.matrix(MERGED_corr))


# subset matrix to wanted comparisons (urine cell types in columns, kidney and bladder cell types in rows)
MERGED_subset <- MERGED.mtx$r[1:(ncol(HK_corr)+ncol(URO_corr)), ((ncol(HK_corr)+ncol(URO_corr))+1):ncol(MERGED_corr)]
colnames(MERGED_subset) <- str_replace(colnames(MERGED_subset), "DCT", "TAL/DCT")
colnames(MERGED_subset) <- str_replace(colnames(MERGED_subset), "DTL", "TEC_VCAM1")
colnames(MERGED_subset) <- str_replace(colnames(MERGED_subset), "TEC_emt", "TEC_inj")

## visualize with corrplot, save ----
png("Supplfig21b.png", width= 800, height =850, type="cairo-png")
corrplot(MERGED_subset, 
         method = "square", type = "full", addCoef.col = "white",
         cl.pos = "r", cl.ratio = 0.1, cl.lim = c(-0.1,1), cl.cex = 1.5,
         tl.col = "black", tl.srt = 90, tl.cex = 2,
         addCoefasPercent = T, 
         order= "original", addgrid.col = F)
dev.off()

