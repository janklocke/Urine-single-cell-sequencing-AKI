# description ----

# This is the R script for generating the whole data Seurat object ("URINE") 
# and subsetted renal cell Seurat object ("RENAL") from raw data shared in 
# "Urinary single-cell sequencing captures intrarenal injury and repair processes 
# in human acute kidney injury" by Klocke et al. 

# loading packages and data ----

## loading packages ----
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(celldex)
library(scDblFinder)
library(harmony)
library(patchwork)
library(celldex)
library(RColorBrewer)
library(viridis)
library(paletteer) 
library(hrbrthemes)

## loading data ----
P001.data <- Read10X(data.dir = "~/UR_AKI_P001/filtered_feature_bc_matrix")
P002.1.data <- Read10X(data.dir = "~/UR_AKI_P002.1/filtered_feature_bc_matrix")
P002.2.data <- Read10X(data.dir = "~/UR_AKI_P002.2/filtered_feature_bc_matrix")
P003.data <- Read10X(data.dir = "~/UR_AKI_P003/filtered_feature_bc_matrix")
P005.data <- Read10X(data.dir = "~/UR_AKI_P005/filtered_feature_bc_matrix")
P006.data <- Read10X(data.dir = "~/UR_AKI_P006/filtered_feature_bc_matrix")
P007.data <- Read10X(data.dir = "~/UR_AKI_P007/filtered_feature_bc_matrix")
P017.1.data <- Read10X(data.dir = "~/UR_AKI_P017.1/filtered_feature_bc_matrix")
P017.2.data <- Read10X(data.dir = "~/UR_AKI_P017.2/filtered_feature_bc_matrix")
P018.1.data <- Read10X(data.dir = "~/UR_AKI_P018.1/filtered_feature_bc_matrix")
P018.2.data <- Read10X(data.dir = "~/UR_AKI_P018.2/filtered_feature_bc_matrix")
P019.1.data <- Read10X(data.dir = "~/UR_AKI_P019.1/filtered_feature_bc_matrix")
P019.2.data <- Read10X(data.dir = "~/UR_AKI_P019.2/filtered_feature_bc_matrix")
P021.data <- Read10X(data.dir = "~/UR_AKI_P021/filtered_feature_bc_matrix")
P022.data <- Read10X(data.dir = "~/UR_AKI_P022/filtered_feature_bc_matrix")
P023.1.data <- Read10X(data.dir = "~/UR_AKI_P023.1/filtered_feature_bc_matrix")
P023.2.data <- Read10X(data.dir = "~/UR_AKI_P023.2/filtered_feature_bc_matrix")
P023.3.data <- Read10X(data.dir = "~/UR_AKI_P023.3/filtered_feature_bc_matrix")
P023.4.data <- Read10X(data.dir = "~/UR_AKI_P023.4/filtered_feature_bc_matrix")
P024.1.data <- Read10X(data.dir = "~/UR_AKI_P024.1/filtered_feature_bc_matrix")
P024.2.data <- Read10X(data.dir = "~/UR_AKI_P024.2/filtered_feature_bc_matrix")


## load seurat objects from demultiplexed samples (see "AKI_Urine_Sediment_Demultiplexing_pooled_samples_script") ----
Pool1 <- readRDS("~/AKI_Urine_sediment_Pool1.rds")
Pool2 <- readRDS("~/AKI_Urine_sediment_Pool2.rds")
Pool3 <- readRDS("~/AKI_Urine_sediment_Pool3.rds")
Pool4 <- readRDS("~/AKI_Urine_sediment_Pool4.rds")
Pool5 <- readRDS("~/AKI_Urine_sediment_Pool5.rds")
Pool6 <- readRDS("~/AKI_Urine_sediment_Pool6.rds")

### remove samples from patients without AKI in patient pools ----
Pool4 <- subset(Pool4, idents = "P054", invert=T)
Pool6 <- subset(Pool6, idents = c( "P116", "P118", "P120"), invert=T)

# create and merge Seurat Objects ----

## Initialize Seurat objects ----
P005 <- CreateSeuratObject(counts = P005.data, project = "P005", min.cells = 3, min.features = 200)
P006 <- CreateSeuratObject(counts = P006.data, project = "P006", min.cells = 3, min.features = 200)
P007 <- CreateSeuratObject(counts = P007.data, project = "P007", min.cells = 3, min.features = 200)
P017.1 <- CreateSeuratObject(counts = P017.1.data, project = "P017.1", min.cells = 3, min.features = 200)
P017.2 <- CreateSeuratObject(counts = P017.2.data, project = "P017.2", min.cells = 3, min.features = 200)
P018.1 <- CreateSeuratObject(counts = P018.1.data, project = "P018.1", min.cells = 3, min.features = 200)
P018.2 <- CreateSeuratObject(counts = P018.2.data, project = "P018.2", min.cells = 3, min.features = 200)
P019.1 <- CreateSeuratObject(counts = P019.1.data, project = "P019.1", min.cells = 3, min.features = 200)
P019.2 <- CreateSeuratObject(counts = P019.2.data, project = "P019.2", min.cells = 3, min.features = 200)
P021 <- CreateSeuratObject(counts = P021.data, project = "P021", min.cells = 3, min.features = 200)
P022 <- CreateSeuratObject(counts = P022.data, project = "P022", min.cells = 3, min.features = 200)
P023.1 <- CreateSeuratObject(counts = P023.1.data, project = "P023.1", min.cells = 3, min.features = 200)
P023.2 <- CreateSeuratObject(counts = P023.2.data, project = "P023.2", min.cells = 3, min.features = 200)
P023.3 <- CreateSeuratObject(counts = P023.3.data, project = "P023.3", min.cells = 3, min.features = 200)
P023.4 <- CreateSeuratObject(counts = P023.4.data, project = "P023.4", min.cells = 3, min.features = 200)
P024.1 <- CreateSeuratObject(counts = P024.1.data, project = "P024.1", min.cells = 3, min.features = 200)
P024.2 <- CreateSeuratObject(counts = P024.2.data, project = "P024.2", min.cells = 3, min.features = 200)
P001 <- CreateSeuratObject(counts = P001.data, project = "P001", min.cells = 3, min.features = 200)
P002.1 <- CreateSeuratObject(counts = P002.1.data, project = "P002.1", min.cells = 3, min.features = 200)
P002.2 <- CreateSeuratObject(counts = P002.2.data, project = "P002.2", min.cells = 3, min.features = 200)
P003 <- CreateSeuratObject(counts = P003.data, project = "P003", min.cells = 3, min.features = 200)

## make a list of all objects, determine percentage of mitochondrial RNA and doublets per sample. ----
URINEList <- list(P005, P006, P007, P017.1, P017.2, P018.1, P018.2, 
                  P019.1, P019.2, P021, P022, P023.1, P023.2, P023.3, 
                  P023.4, P024.1, P024.2, P001, P002.1, P002.2, P003, 
                  Pool1, Pool2, Pool3, Pool4, Pool5, Pool6)

for (i in 1:length(URINEList)) {
  URINEList[[i]][["percent.mt"]] <- PercentageFeatureSet(URINEList[[i]], pattern = "^MT-")
  doublets <- scDblFinder(GetAssayData(URINEList[[i]], assay = "RNA", slot = "data"))
  doublets <- as.vector(doublets@colData@listData[["scDblFinder.class"]])
  URINEList[[i]]@meta.data$multiplet_class <- doublets
  URINEList[[i]]@project.name <- levels(URINEList[[i]]@active.ident)
  URINEList[[i]] <- RenameCells(URINEList[[i]], add.cell.id = paste0(URINEList[[i]]$orig.ident, "_"))
}

## merge all the objects in the list ----
URINE <- purrr::reduce(URINEList, merge)

# Add clinical meta data ----

## create empty metadata slots ----
URINE$patient <- NA
URINE$gender <- NA
URINE$age <- NA
URINE$gating <- NA
URINE$sampling <- NA
URINE$COVID19_infection <- NA
URINE$CKD <- NA
URINE$KDIGO_AKI_stage <- NA
URINE$outcome_after_ninety_days <- NA
URINE$RRT <- NA
URINE$AKI_type <- NA
URINE$AKI_timepoint.main <- NA
URINE$AKI_timepoint.fine <- NA

## load clinical data ----
clin.data <- read.delim("~/URINE_AKI_patient_data.csv", sep = ";")
clin.data$sample.pseudonym
colnames(clin.data)

## add clinical metadata to SO ----
for (i in clin.data$sample.pseudonym) {
  URINE@meta.data[URINE@meta.data$orig.ident %in% i, ]$patient <- clin.data$patient[clin.data$sample.pseudonym %in% i]
  URINE@meta.data[URINE@meta.data$orig.ident %in% i, ]$gender <- clin.data$gender[clin.data$sample.pseudonym %in% i]
  URINE@meta.data[URINE@meta.data$orig.ident %in% i, ]$age <- clin.data$age[clin.data$sample.pseudonym %in% i]
  URINE@meta.data[URINE@meta.data$orig.ident %in% i, ]$gating <- clin.data$gating[clin.data$sample.pseudonym %in% i]
  URINE@meta.data[URINE@meta.data$orig.ident %in% i, ]$sampling <- clin.data$sampling[clin.data$sample.pseudonym %in% i]
  URINE@meta.data[URINE@meta.data$orig.ident %in% i, ]$COVID19_infection <- clin.data$COVID19_infection[clin.data$sample.pseudonym %in% i]
  URINE@meta.data[URINE@meta.data$orig.ident %in% i, ]$CKD <- clin.data$CKDgrade[clin.data$sample.pseudonym %in% i]
  URINE@meta.data[URINE@meta.data$orig.ident %in% i, ]$KDIGO_AKI_stage <- clin.data$KDIGO_AKI_stage[clin.data$sample.pseudonym %in% i]
  URINE@meta.data[URINE@meta.data$orig.ident %in% i, ]$outcome_after_ninety_days <- clin.data$outcome.main[clin.data$sample.pseudonym %in% i]
  URINE@meta.data[URINE@meta.data$orig.ident %in% i, ]$RRT <- clin.data$RRT[clin.data$sample.pseudonym %in% i]
  URINE@meta.data[URINE@meta.data$orig.ident %in% i, ]$AKI_type <- clin.data$AKI.type[clin.data$sample.pseudonym %in% i]
  URINE@meta.data[URINE@meta.data$orig.ident %in% i, ]$AKI_timepoint.main <- clin.data$AKI.timepoint.main[clin.data$sample.pseudonym %in% i]
  URINE@meta.data[URINE@meta.data$orig.ident %in% i, ]$AKI_timepoint.fine <- clin.data$AKI.timepoint.fine[clin.data$sample.pseudonym %in% i]
}

# QC analysis and subsetting ----

a <- VlnPlot(URINE, c("nFeature_RNA"), group.by = "orig.ident", pt.size = 0.2) + 
  NoLegend() + 
  geom_hline(yintercept = 200, linetype = "solid", color = "red", size = 1) + 
  geom_hline(yintercept = 4000, linetype = "solid", color = "red", size = 1)
b <- VlnPlot(URINE, c("nCount_RNA"), group.by = "orig.ident", pt.size = 0.2) + 
  NoLegend() + 
  geom_hline(yintercept = 500, linetype = "solid", color = "red", size = 1) + 
  geom_hline(yintercept = 50000, linetype = "solid", color = "red", size = 1)
c <- VlnPlot(URINE, c("percent.mt"), group.by = "orig.ident", pt.size = 0.2) + 
  NoLegend() + 
  geom_hline(yintercept = 20, linetype = "solid", color = "red", size = 1) + 
  geom_hline(yintercept = 10, linetype = "longdash", color = "red", size = 1)

a / b / c

## subset by percent.mt, singlets, counts and features ----
URINE_noQC <- URINE

###################
# since the scDblFinder() function has slightly variable output, the list of quality filtered cells in URINE may vary based on which cells 
# have been assigned as singlets/doublets. To guarantee an output equal to the one shown in the paper, we provide the URINE Seurat object post QC via
# figshare: https://figshare.com/articles/dataset/KI_2022_Klocke_et_al_SO_all_urine_cells_rds/22567201
# otherwise you can continue with the code getting a similar albeit slightly differing object. 
###################

URINE <- subset(URINE, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 500 & nCount_RNA < 50000 & percent.mt < 20 & multiplet_class == "singlet")

## normalization, dim. reduction and clustering 1 ----
## in the first round of clustering, we will consider all cells with percent.mt < 20% and 
## determine clusters of leukocytes and urogenital cells, which should have a cut-off of 
## mitochondrial RNA of 10% for further analysis, as opposed to kidney cells (20%)

## perform normalization, scaling, etc. with SCTransform
URINE <- SCTransform(object = URINE, vars.to.regress = "percent.mt", method = "glmGamPoi", verbose = T)
URINE <- RunPCA(object = URINE, dims = 1:45, verbose = T)

## perform batch correction with harmony
## tau and theta are adjusted to prevent overcompenstation of batch-effects and overclustering of small vs. large samples 
URINE <- RunHarmony(URINE, group.by.vars = c("orig.ident"), lambda = 2, tau = 1000, theta = 1, assay.use = "SCT")

## run dimensionality reduction and find clusters
URINE <- RunUMAP(object = URINE, reduction = "harmony", dims = 1:45, verbose = T)
URINE <- FindNeighbors(object = URINE, reduction = "harmony", dims = 1:45, verbose = T)
URINE <- FindClusters(object = URINE, resolution = 0.3, verbose = T)

## Find differentially expressed genes per cluster
URINE.markers <- FindAllMarkers(URINE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
URINE.top10 <- URINE.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

## In clusters of leukocyte and urogenital origin, find cells with mitochondiral RNA > 10%
highmito <- WhichCells(URINE, idents = c(0, 1, 2, 6, 7, 8, 9, 11, 13, 14, 18), expression = percent.mt > 10, invert = F)
DimPlot(URINE, reduction = "umap", label = T, cells.highlight = highmito, cols = custom_colors[[1]]) + NoLegend()

## normalization, dim. reduction and clustering 2 ----
## repeat dim. red. and clustering, but remove leukocytes/urogenital cells with >10% mitochondrial RNA
URINE <- subset(URINE_noQC, cells = highmito, invert = T)
## subset by percent.mt, singlets and features
URINE <- subset(URINE, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 500 & nCount_RNA < 50000 & percent.mt < 20 & multiplet_class == "singlet")

# gather simple info on SO ----
URINE.info <- tibble(URINE_info = c("Total cells", 
                                    "Total genes", 
                                    "median cells/patient", 
                                    "median UMIs/cell", 
                                    "median genes/cell", 
                                    "median  %mtRNA/cell"), values = round(c(ncol(URINE),
                                                                             nrow(URINE),
                                                                             median(table(URINE$orig.ident)),
                                                                             median(URINE$nCount_RNA),
                                                                             median(URINE$nFeature_RNA),
                                                                             median(URINE$percent.mt)), digits = 1))       



# perform normalization, scaling, etc. with SCTransform
URINE <- SCTransform(object = URINE, vars.to.regress = "percent.mt", method = "glmGamPoi", verbose = T)
URINE <- RunPCA(object = URINE, dims = 1:40, verbose = T)
DimHeatmap(URINE, dims = 1:15, cells = 500, balanced = TRUE)

#perform batch correction with harmony
# tau and theta are adjusted to prevent overcompenstation of batch-effects and overclustering of small vs. large samples 
URINE <- RunHarmony(URINE, group.by.vars = c("orig.ident"), lambda = 2, tau = 1000, theta = 1, assay.use = "SCT")

#run dimensionality reduction and find clusters
URINE <- RunUMAP(object = URINE, reduction = "harmony", dims = 1:45, verbose = T)
URINE <- RunTSNE(object = URINE, reduction = "harmony", dims = 1:45, verbose = T)
URINE <- FindNeighbors(object = URINE, reduction = "harmony", dims = 1:45, verbose = T)
URINE <- FindClusters(object = URINE, resolution = 0.3, verbose = T)

# find differentially expressed genes in URINE ----
#add cell cycle scoring
URINE <- CellCycleScoring(URINE, assay = 'SCT', s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)

#find differentially expressed genes
## DEG in Supplemental Table 2 ----
URINE.markers <- FindAllMarkers(URINE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
URINE.top20 <- URINE.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
URINE.top5 <- URINE.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)


# annotation of URINE clusters ----
# assign names and numbers to clusters, rearrange cluster order
URINE <- SetIdent(URINE, value ="seurat_clusters")
URINEClusterNames <- c("UGEC1", "MO_infl", "UGEC2", "TEC_prg", "TEC_dmg1", 
                       "MO_kdnrs", "TEC_inj", "Tcells", "MO_MT+", "ERY", 
                       "MO_dmg", "MO_SPP1+", "GRAN", "PDC", "TEC_dmg2", 
                       "TEC_CD", "TEC_prlf", "TEC_PT", "Bcells")

# URINE$celltype includes names of clusters
names(URINEClusterNames) <- levels(URINE)
URINE <- RenameIdents(URINE, URINEClusterNames)
URINE$celltype <- URINE@active.ident
URINE <- SetIdent(URINE, value ="seurat_clusters")
URINE$celltype <- factor(URINE$celltype, levels(URINE$celltype)[c(14,18,16,7,5,15,17,4,2,6,9,11,12,13,8,19,10,1,3)])
levels(URINE$celltype)

# URINE$celltype_nr includes new numbers of clusters ordered by origin
URINE <- SetIdent(URINE, value ="celltype")
URINEClusterNrs  <- paste(1:19)
names(URINEClusterNrs) <- levels(URINE)
URINE <- RenameIdents(URINE, URINEClusterNrs)
URINE$celltype_nr <- URINE@active.ident
URINE <- SetIdent(URINE, value ="seurat_clusters")
levels(URINE$celltype_nr)

# assigning cellgroups for broad comparisons of differing urinary cell profiles, found in "URINE$cellgroup"
URINE$cellgroup <- as.character(URINE$celltype)
URINE@meta.data[URINE@meta.data$celltype %in% c("TEC_prg", "TEC_inj", 
                                                "TEC_dmg1", "TEC_dmg2", 
                                                "TEC_CD", "TEC_prlf", "TEC_PT"),]$cellgroup <- "TEC"
URINE@meta.data[URINE@meta.data$celltype %in% c("MO_dmg", "MO_infl", "MO_MT+", 
                                                "MO_SPP1+", "MO_kdnrs", 
                                                "Tcells", "Bcells", "GRAN"),]$cellgroup <- "LEUK"
URINE@meta.data[URINE@meta.data$celltype %in% c("UGEC1", "UGEC2"),]$cellgroup <- "UGEC"

DimPlot(URINE, group.by = "cellgroup")


# Subset renal cells to Seurat object RENAL ----

# adjust cluster resolution
URINE <- FindClusters(object = URINE, resolution = 0.5, verbose = T)

# find clusters with kidney markers (EPCAM, CRYAB, PAX8) but without urogenital or leukocyte markers
DotPlot(URINE, features = c("NPHS2", "EPCAM", "CRYAB", "PAX8",  "HBA1", "PSCA", "UPK2", "LYZ", "CD3E", "CD79A"), col.min = 0) + 
  theme(axis.text.y = element_text(size = 10), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=0.95, vjust=0.2),
        legend.title = element_text(size =10),
        legend.text = element_text(size = 10)) + 
  coord_flip()

# subset kidney parenchymal clusters
RENAL <- subset(URINE, idents = c(3, 4, 5, 11, 17, 18, 19, 20, 22, 23, 24))


# normalization, dim reduction and clustering for RENAL ----

RENAL@active.assay <- "RNA"

# perform normalization, scaling, etc. with SCTransform
RENAL <- SCTransform(object = RENAL, vars.to.regress = "percent.mt", method = "glmGamPoi", verbose = T)
RENAL <- RunPCA(object = RENAL, dims = 1:30, verbose = T)

#perform batch correction with harmony
RENAL <- RunHarmony(RENAL, group.by.vars = "orig.ident", assay.use = "SCT")

#run dimensionality reduction and find clusters
RENAL <- RunUMAP(object = RENAL, reduction = "harmony", dims = 1:30, verbose = T)
RENAL <- RunTSNE(object = RENAL, reduction = "harmony", dims = 1:30, verbose = T)
RENAL <- FindNeighbors(object = RENAL, reduction = "harmony", dims = 1:30, verbose = T)
RENAL <- FindClusters(object = RENAL, resolution = 0.5, verbose = T)

# Visualization RENAL ----

#visualization
DimPlot(RENAL, reduction = "umap", label = T, repel = F, pt.size = 1, label.size = 5) +ggtitle("urinary renal cells") + NoLegend()

# find remaining leukocytge and urogenital clusters
DotPlot(RENAL, features = c("NPHS2", "EPCAM", "CRYAB", "PAX8",  "HBA1", "PSCA", "UPK2", "LYZ", "CD3E", "CD79A"), col.min = 0) + 
  theme(axis.text.y = element_text(size = 10), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=0.95, vjust=0.2),
        legend.title = element_text(size =10),
        legend.text = element_text(size = 10)) + 
  coord_flip()

# remove remaining leukocyte and urogenital cluster
RENAL <- subset(RENAL, idents = c(9,10), invert=T)


# find differentially expressed genes in URINE ----

## DEG in Supplemental Table 3 ----
RENAL.markers <- FindAllMarkers(RENAL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
RENAL.top20 <- RENAL.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

# annotation of RENAL clusters ----
# assign names and numbers to clusters, rearrange cluster order
RENALClusterNames <- c("TEC_prg1", "TEC_MMP7", "TEC_prg2", "TEC_TXNRD1", 
                       "TEC_GCLM", "TEC_IL32", "TEC_prg4", "TAL/DCT", 
                       "PT", "TEC_prlf", "ATL", "DTL", "PDC_inj", 
                       "TEC_LIF", "TEC_dmg", "CD-PC", "PDC", "CNT/CD", 
                       "TEC_prg3", "TAL", "CD-IC", "TEC_MT", "TEC_IFIT")

# RENAL$renalcelltype includes names of clusters
RENAL <- SetIdent(RENAL, value ="seurat_clusters")
names(RENALClusterNames) <- levels(RENAL)
RENAL <- RenameIdents(RENAL, RENALClusterNames)
RENAL$renalcelltype <- RENAL@active.ident
RENAL <- SetIdent(RENAL, value ="seurat_clusters")
RENAL$renalcelltype <- factor(RENAL$renalcelltype, levels(RENAL$renalcelltype)[c(17,13,9,12,11,20,8,18,16,21,22,23,2,6,14,4,5,15,10,1,3,19,7)])

# RENAL$renalcelltype_nr includes new numbers of clusters ordered by origin
RENAL <- SetIdent(RENAL, value ="renalcelltype")
RENALClusterNrs  <- paste(1:23)
names(RENALClusterNrs) <- levels(RENAL)
RENAL <- RenameIdents(RENAL, RENALClusterNrs)
RENAL$renalcelltype_nr <- RENAL@active.ident
RENAL <- SetIdent(RENAL, value ="seurat_clusters")
levels(RENAL$renalcelltype_nr)

# add metadata combining several renalcelltypes into groups, found in "RENAL$renalcellgroup"
RENAL$renalcellgroup <- as.character(RENAL$renalcelltype)
RENAL@meta.data[RENAL@meta.data$renalcelltype %in% c("TEC_prg1", "TEC_prg2", "TEC_prg3", "TEC_prg4"),]$renalcellgroup <- "TEC_prg"
RENAL@meta.data[RENAL@meta.data$renalcelltype %in% c("TEC_MMP7", "TEC_IL32", "TEC_LIF", "TEC_MT", "TEC_IFIT"),]$renalcellgroup <- "TEC_emt"
RENAL@meta.data[RENAL@meta.data$renalcelltype %in% c("TEC_TXNRD1", "TEC_GCLM", "TEC_dmg"),]$renalcellgroup <- "TEC_str"
RENAL@meta.data[RENAL@meta.data$renalcelltype %in% c("PDC", "PDC_inj"),]$renalcellgroup <- "PDC"
RENAL@meta.data[RENAL@meta.data$renalcelltype %in% c("PT", "DTL", "ATL", "TAL", "TAL/DCT", "CNT/CD", "CD-PC", "CD-IC"),]$renalcellgroup <- "TEC_cnn"

RENAL$renalcellgroup <- factor(RENAL$renalcellgroup, levels(factor(RENAL$renalcellgroup))[c(1,2,3,6,5,4)])

DimPlot(RENAL, group.by = "renalcellgroup")

# save objects  ----

# Seurat object containing all urine cells ("URINE")
saveRDS(URINE, "~/SO_all_urine_cells.rds")
# table of DEG of all URINE clusters
saveRDS(URINE.markers, "~/DEG_all_urine_cells.rds")

# Seurat object containing all urine cells before quality control filtering("URINE_noQC")
saveRDS(URINE_noQC, "~/SO_all_urine_cells_noQC.rds")

# Seurat object containing kidney urine cells ("RENAL")
saveRDS(RENAL, "~/SO_kidney_urine_cells.rds")
# table of DEG of all RENAL clusters
saveRDS(RENAL.markers, "~/DEG_kidney_urine_cells.rds")
