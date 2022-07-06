# description ----

# This is the R script for generating Supplemental Figure 16,18,19,20 in 
# "Urinary single-cell sequencing captures intrarenal injury and repair processes 
# in human acute kidney injury" by Klocke et al. 

# loading packages and data ----

## loading packages ----
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(celldex)
library(patchwork)
library(RColorBrewer)
library(viridis)
library(paletteer) 
library(hrbrthemes)
library(scDblFinder)
library(harmony)
library(ggpubr)


## loading color schemes ----
renalcelltype_colors <- 
  c(paletteer_d("ggsci::purple_material")[c(3,5)],                  #PDC colors
    paletteer_d("ggsci::teal_material")[c(8,6,4,3)],                #healthy TEC PT and LOH colors   
    paletteer_d("ggsci::light_blue_material")[c(3,5,7,9)],          #healthy TEC DT and CD colors
    paletteer_d("ggsci::amber_material")[c(3,4,6,8,10)],            #EMT/infl TEC
    paletteer_d("ggsci::red_material")[c(2,5,8)],                   #OxyStress TEC
    paletteer_d("ggsci::pink_material")[4],                         #prolif TEC
    paletteer_d("ggsci::light_green_material")[c(9,7,5,3)])         #progenitors 

renalcellgroup_colors <- renalcelltype_colors[c(1,8,14,17,19,21)]
names(renalcelltype_colors) <- levels(factor(RENAL$renalcelltype))
names(renalcellgroup_colors) <- levels(factor(RENAL$renalcellgroup))
celltype_nr_colors <- c(colorRampPalette(paletteer_d("palettesForR::Pastels"))(19)[1], "#c98ed4", colorRampPalette(paletteer_d("palettesForR::Pastels"))(19)[2:19])
disease_celltype_colors <- c(renalcelltype_colors, renalcellgroup_colors, "#c98ed4")
names(disease_celltype_colors)[30] <- "MSC"
names(disease_celltype_colors)[3] <- "PT"
disease_celltype_colors <- disease_celltype_colors[names(disease_celltype_colors) %in% TEC$renalcelltype]
disease_celltype_colors <- disease_celltype_colors[c(1:10,13)]
names(disease_celltype_colors) <- paste(1:11)
## loading data ----

### load own raw data ----
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

### load seurat objects from demultiplexed samples (see "AKI_Urine_Sediment_Demultiplexing_pooled_samples_script") ----
Pool1 <- readRDS("~/AKI_Urine_sediment_Pool1.rds")
Pool2 <- readRDS("~/AKI_Urine_sediment_Pool2.rds")
Pool3 <- readRDS("~/AKI_Urine_sediment_Pool3.rds")
Pool4 <- readRDS("~/AKI_Urine_sediment_Pool4.rds")
Pool5 <- readRDS("~/AKI_Urine_sediment_Pool5.rds")
Pool6 <- readRDS("~/AKI_Urine_sediment_Pool6.rds")

#### remove samples from patients without AKI in patient pools ----
Pool4 <- subset(Pool4, idents = "P054", invert=T)
Pool6 <- subset(Pool6, idents = c( "FlUriCov116", "FlUriCov118", "FlUriCov120"), invert=T)

### load FSGS urine single-cell sequencing data (NCBI GEO, dataset GSE176465) ----
USC_DK_001.data <- Read10X(data.dir = "~/USC_DK_001/raw_gene_bc_matrices/GRCh38")
USC_DK_002.data <- Read10X(data.dir = "~/USC_DK_002/raw_gene_bc_matrices/GRCh38")
USC_DK_003.data <- Read10X(data.dir = "~/USC_DK_003/raw_gene_bc_matrices/GRCh38")
USC_DK_004.data <- Read10X(data.dir = "~/USC_DK_004/raw_gene_bc_matrices/GRCh38")
USC_DK_005.data <- Read10X(data.dir = "~/USC_DK_005_Fixed/raw_gene_bc_matrices/GRCh38")
USC_DK_006.data <- Read10X(data.dir = "~/USC_DK_006/raw_gene_bc_matrices/GRCh38")
USC_DK_007.data <- Read10X(data.dir = "~/USC_DK_007/raw_gene_bc_matrices/GRCh38")
USC_DK_008.data <- Read10X(data.dir = "~/USC_DK_008/raw_gene_bc_matrices/GRCh38")
USC_DK_009.data <- Read10X(data.dir = "~/USC_DK_009/raw_gene_bc_matrices/GRCh38")
USC_DK_010.data <- Read10X(data.dir = "~/USC_DK_010/raw_gene_bc_matrices/GRCh38")
USC_DK_011.data <- Read10X(data.dir = "~/USC_DK_011_Fixed/raw_gene_bc_matrices/GRCh38")
USC_DK_012.data <- Read10X(data.dir = "~/USC_DK_012/raw_gene_bc_matrices/GRCh38")
USC_DK_013.data <- Read10X(data.dir = "~/USC_DK_013/raw_gene_bc_matrices/GRCh38")
USC_DK_014.data <- Read10X(data.dir = "~/USC_DK_014/raw_gene_bc_matrices/GRCh38")
USC_DK_015.data <- Read10X(data.dir = "~/USC_DK_015/raw_gene_bc_matrices/GRCh38")
USC_DK_016.data <- Read10X(data.dir = "~/USC_DK_016/raw_gene_bc_matrices/GRCh38")
USC_DK_018.data <- Read10X(data.dir = "~/USC_DK_018/raw_gene_bc_matrices/GRCh38")
USC_DK_019.data <- Read10X(data.dir = "~/USC_DK_019/raw_gene_bc_matrices/GRCh38")
USC_DK_020.data <- Read10X(data.dir = "~/USC_DK_020/raw_gene_bc_matrices/GRCh38")
USC_DK_021.data <- Read10X(data.dir = "~/USC_DK_021_PL1/raw_gene_bc_matrices/GRCh38")
USC_DK_022.data <- Read10X(data.dir = "~/USC_DK_022_PL2/raw_gene_bc_matrices/GRCh38")
USC_DK_023.data <- Read10X(data.dir = "~/USC_DK_023_LA1/raw_gene_bc_matrices/GRCh38")
USC_DK_024.data <- Read10X(data.dir = "~/USC_DK_024_LA2/raw_gene_bc_matrices/GRCh38")

### load Diabetes urine single-cell dataset ----

# raw data matrices by sample were not available, downloaded merged matrix instead --> no by-patient data 
DM.data <-  read.delim2("~/Human_urine_sigle_cell_datamatrix.txt")

### load external AKI samples ----

COVAKI1.data <- Read10X(data.dir = "~/Forschung/Rohdaten/GSE180595_urine_AKI/p01_AKI_COVID")
COVAKI2.data <- Read10X(data.dir = "~/Forschung/Rohdaten/GSE180595_urine_AKI/p02_AKI_COVID")
COVAKI3.data <- Read10X(data.dir = "~/Forschung/Rohdaten/GSE180595_urine_AKI/p03_AKI_COVID")
COVAKI4.data <- Read10X(data.dir = "~/Forschung/Rohdaten/GSE180595_urine_AKI/p04_AKI_COVID")
COVAKI5.data <- Read10X(data.dir = "~/Forschung/Rohdaten/GSE180595_urine_AKI/p05_AKI_COVID")
COVAKI9.data <- Read10X(data.dir = "~/Forschung/Rohdaten/GSE180595_urine_AKI/p09_AKI")
COVAKI10.data <- Read10X(data.dir = "~/Forschung/Rohdaten/GSE180595_urine_AKI/p10_AKI")
COVAKI11.data <- Read10X(data.dir = "~/Forschung/Rohdaten/GSE180595_urine_AKI/p11_AKI")
COVAKI12.data <- Read10X(data.dir = "~/Forschung/Rohdaten/GSE180595_urine_AKI/p12_AKI")


# create and merge Seurat Objects ----
## Initialize Seurat objects ----

#AKI samples
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

#FSGS samples
USC_DK_001 <- CreateSeuratObject(counts = USC_DK_001.data, project = "USC_DK_001", min.cells = 3, min.features = 200)
USC_DK_002 <- CreateSeuratObject(counts = USC_DK_002.data, project = "USC_DK_002", min.cells = 3, min.features = 200)
USC_DK_003 <- CreateSeuratObject(counts = USC_DK_003.data, project = "USC_DK_003", min.cells = 3, min.features = 200)
USC_DK_004<- CreateSeuratObject(counts = USC_DK_004.data, project = "USC_DK_004", min.cells = 3, min.features = 200)
USC_DK_005 <- CreateSeuratObject(counts = USC_DK_005.data, project = "USC_DK_005", min.cells = 3, min.features = 200)
USC_DK_006 <- CreateSeuratObject(counts = USC_DK_006.data, project = "USC_DK_006", min.cells = 3, min.features = 200)
USC_DK_007 <- CreateSeuratObject(counts = USC_DK_007.data, project = "USC_DK_007", min.cells = 3, min.features = 200)
USC_DK_008 <- CreateSeuratObject(counts = USC_DK_008.data, project = "USC_DK_008", min.cells = 3, min.features = 200)
USC_DK_009 <- CreateSeuratObject(counts = USC_DK_009.data, project = "USC_DK_009", min.cells = 3, min.features = 200)
USC_DK_010 <- CreateSeuratObject(counts = USC_DK_010.data, project = "USC_DK_010", min.cells = 3, min.features = 200)
USC_DK_011 <- CreateSeuratObject(counts = USC_DK_011.data, project = "USC_DK_011", min.cells = 3, min.features = 200)
USC_DK_012 <- CreateSeuratObject(counts = USC_DK_012.data, project = "USC_DK_012", min.cells = 3, min.features = 200)
USC_DK_013 <- CreateSeuratObject(counts = USC_DK_013.data, project = "USC_DK_013", min.cells = 3, min.features = 200)
USC_DK_014<- CreateSeuratObject(counts = USC_DK_014.data, project = "USC_DK_014", min.cells = 3, min.features = 200)
USC_DK_015 <- CreateSeuratObject(counts = USC_DK_015.data, project = "USC_DK_015", min.cells = 3, min.features = 200)
USC_DK_016 <- CreateSeuratObject(counts = USC_DK_016.data, project = "USC_DK_016", min.cells = 3, min.features = 200)
USC_DK_018 <- CreateSeuratObject(counts = USC_DK_018.data, project = "USC_DK_018", min.cells = 3, min.features = 200)
USC_DK_019 <- CreateSeuratObject(counts = USC_DK_019.data, project = "USC_DK_019", min.cells = 3, min.features = 200)
USC_DK_020 <- CreateSeuratObject(counts = USC_DK_020.data, project = "USC_DK_020", min.cells = 3, min.features = 200)
USC_DK_021 <- CreateSeuratObject(counts = USC_DK_021.data, project = "USC_DK_021", min.cells = 3, min.features = 200)
USC_DK_022 <- CreateSeuratObject(counts = USC_DK_022.data, project = "USC_DK_022", min.cells = 3, min.features = 200)
USC_DK_023 <- CreateSeuratObject(counts = USC_DK_023.data, project = "USC_DK_023", min.cells = 3, min.features = 200)
USC_DK_024 <- CreateSeuratObject(counts = USC_DK_024.data, project = "USC_DK_024", min.cells = 3, min.features = 200)

# DM sample
DM <- CreateSeuratObject(counts = DM.data, project = "DM", min.cells = 3, min.features = 200)


# external AKI samples 
COVAKI1 <- CreateSeuratObject(counts = COVAKI1.data, project = "COVAKI01", min.cells = 3, min.features = 200)
COVAKI2<- CreateSeuratObject(counts = COVAKI2.data, project = "COVAKI02", min.cells = 3, min.features = 200)
COVAKI3 <- CreateSeuratObject(counts = COVAKI3.data, project = "COVAKI03", min.cells = 3, min.features = 200)
COVAKI4 <- CreateSeuratObject(counts = COVAKI4.data, project = "COVAKI04", min.cells = 3, min.features = 200)
COVAKI5 <- CreateSeuratObject(counts = COVAKI5.data, project = "COVAKI05", min.cells = 3, min.features = 200)
COVAKI6 <- CreateSeuratObject(counts = COVAKI6.data, project = "COVAKI06", min.cells = 3, min.features = 200)
COVAKI7 <- CreateSeuratObject(counts = COVAKI7.data, project = "COVAKI07", min.cells = 3, min.features = 200)
COVAKI8 <- CreateSeuratObject(counts = COVAKI8.data, project = "COVAKI08", min.cells = 3, min.features = 200)
COVAKI9 <- CreateSeuratObject(counts = COVAKI9.data, project = "COVAKI09", min.cells = 3, min.features = 200)
COVAKI10 <- CreateSeuratObject(counts = COVAKI10.data, project = "COVAKI10", min.cells = 3, min.features = 200)
COVAKI11 <- CreateSeuratObject(counts = COVAKI11.data, project = "COVAKI11", min.cells = 3, min.features = 200)
COVAKI12 <- CreateSeuratObject(counts = COVAKI12.data, project = "COVAKI12", min.cells = 3, min.features = 200)

## make a list of objects based on their different AKI etiology or disease. ----
## CSA = cardiac surgery assoc. AKI
## SEP = septic penumonia assoc. AKI
## PRE = other prerenal causes of AKI
# determine percentage of mitochondrial RNA and doublets per sample. ----

CSAList <- list(ANV05, ANV06, ANV17.1, ANV17.2, ANV19.1, ANV19.2, ANV21, ANV22, ANV23.1, ANV23.2, ANV23.3, ANV23.4, ANV24.1, ANV24.2, ANV0812neg, ANV0812pos)
SEPList <- list(ANV18.1, ANV18.2, CoV01, CoV02.1, CoV02.2, CoV03, Pool3, Pool4, Pool5, Pool6)
PREList <- list(ANV07, ANV21, ANV22, Pool1, Pool2)
FSGSList <- c(USC_DK_001,USC_DK_002,USC_DK_003,USC_DK_004,USC_DK_005,USC_DK_006,
              USC_DK_007,USC_DK_008,USC_DK_009,USC_DK_010,USC_DK_011,USC_DK_012,
              USC_DK_013,USC_DK_014,USC_DK_015,USC_DK_016,USC_DK_018,USC_DK_019,
              USC_DK_020,USC_DK_021,USC_DK_022,USC_DK_023,USC_DK_024)
COVAKIList <- c(COVAKI1, COVAKI2, COVAKI3, COVAKI4, COVAKI5, COVAKI9, COVAKI10, COVAKI11, COVAKI12)

### determine percent mitochondrial RNA, doublets and merging of samples. ----
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

CSA <- compileSO(CSAList)
SEP <- compileSO(SEPList)
PRE <- compileSO(PREList)
FSGS <- compileSO(FSGSList)
DM <- compileSO(list(DM))
COVAKI <- compileSO(COVAKIList)

# QC filtering, normalization, dim. red. and clustering ----
CSA0 <- CSA
SEP0 <- SEP
PRE0 <- PRE
FSGS0 <- FSGS
DM0 <- DM
COVAKI0 <- COVAKI
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

CSA <- processSO(CSA)
SEP <- processSO(SEP)
PRE <- processSO(PRE)
FSGS <- processSO(FSGS)
DM <- processSO(DM)
COVAKI <- processSO(COVAKI)

## clustering ----
CSA <- FindClusters(object = CSA, resolution = 0.6, verbose = T)
SEP <- FindClusters(object = SEP, resolution = 0.7, verbose = T)
PRE <- FindClusters(object = PRE, resolution = 0.5, verbose = T)
DM <- FindClusters(object = DM, resolution = 0.5, verbose = T)
FSGS@active.assay <- "SCT"
FSGS <- FindClusters(object = FSGS, resolution = 0.5, verbose = T)
COVAKI <- FindClusters(object = COVAKI, resolution = 0.6, verbose = T)
## find marker genes ----
CSA.markers <- FindAllMarkers(CSA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SEP.markers <- FindAllMarkers(SEP, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
PRE.markers <- FindAllMarkers(PRE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
DM.markers <- FindAllMarkers(DM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
FSGS.markers <- FindAllMarkers(FSGS, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
COVAKI.markers <- FindAllMarkers(COVAKI, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

CSA.top10 <- CSA.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
SEP.top10 <- SEP.Markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
PRE.top10 <- PRE.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DM.top10 <- DM.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
FSGS.top10 <- FSGS.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

markers <- COVAKI.markers %>% 
  group_by(cluster) %>% 
  top_n(2, avg_log2FC)
DotPlot(COVAKI, features = c(unique(markers$gene)))+ coord_flip()
saveRDS(COVAKI, "~/20220906COVAKI.rds")

DM$patient <- "DM_Healthy"
DM@meta.data[grep("S1155843", rownames(DM@meta.data)),]$patient <- "DM_P01_1"
DM@meta.data[grep("S1155844", rownames(DM@meta.data)),]$patient <- "DM_P01_2"
DM@meta.data[grep("S1155834", rownames(DM@meta.data)),]$patient <- "DM_P01_3"
DM@meta.data[grep("S1155835", rownames(DM@meta.data)),]$patient <- "DM_P01_4"
DM@meta.data[grep("S1155841", rownames(DM@meta.data)),]$patient <- "DM_P02_1"
DM@meta.data[grep("S1155842", rownames(DM@meta.data)),]$patient <- "DM_P02_2"
DM@meta.data[grep("S1155832", rownames(DM@meta.data)),]$patient <- "DM_P02_3"
DM@meta.data[grep("S1155833", rownames(DM@meta.data)),]$patient <- "DM_P02_4"
DM@meta.data[grep("S1155839", rownames(DM@meta.data)),]$patient <- "DM_P03_1"
DM@meta.data[grep("S1155840", rownames(DM@meta.data)),]$patient <- "DM_P03_2"
DM@meta.data[grep("S1155830", rownames(DM@meta.data)),]$patient <- "DM_P03_3"
DM@meta.data[grep("S1155831", rownames(DM@meta.data)),]$patient <- "DM_P03_4"
DM@meta.data[grep("S1155837", rownames(DM@meta.data)),]$patient <- "DM_P04_1"
DM@meta.data[grep("S1155838", rownames(DM@meta.data)),]$patient <- "DM_P04_2"
DM@meta.data[grep("S1155828", rownames(DM@meta.data)),]$patient <- "DM_P04_3"
DM@meta.data[grep("S1155829", rownames(DM@meta.data)),]$patient <- "DM_P04_4"
DM@meta.data[grep("S1155836", rownames(DM@meta.data)),]$patient <- "DM_P05"
table(DM$patient)

DM$orig.ident <- DM$patient

DM <- SetIdent(DM , value = "patient")
DM <- subset(DM, ident = "DM_Healthy", inv=T)

# annotation of clusters ----
## cluster annotation was done manually by examination of DEG and review of previously published marker genes 
CSAClusterNames <- c("UGEC1", "UGEC2", "MO_infl", "TEC_prg", "TEC_inj", "TEC_str", "MO_infl", "MO_kdnrs", "MO_MT+", "TEC_inj", "UGEC1", "UGEC2", "GRAN", "Tcells", "MO_dmg", "TEC_dmg", "TEC_prg", "MO_dmg", "TEC_prlf", "TEC_PT", "TEC_CD", "TEC_prg", "UGEC2", "TEC_CD", "PDC", "ERY") 
SEPClusterNames <- c("UGEC1", "MO_infl", "TEC_str", "Tcells", "TEC_prg", "MO_infl", "TEC_inj", "UGEC2", "TEC_inj", "TEC_str", "TEC_str", "Tcells", "MO_MT+", "PDC", "MO_kdnrs", "MO_dmg", "Bcells", "GRAN", "TEC_CD", "TEC_prlf", "PDC", "GRAN", "Tcells", "TEC_PT")
PREClusterNames <- c("MO_kdnrs", "UGEC1", "MO_infl", "UGEC1", "MO_kdnrs", "MO_dmg", "MO_MT+", "TEC_inj", "GRAN", "MO_SPP1+", "UGEC2", "TEC_prg", "Tcells", "Tcells")
DMClusterNames <- c("UGEC1", "UGEC1", "UGEC1", "UGEC1", "UGEC1", 
                    "MO_kdnrs", "UGEC1", "TEC_inj", "UGEC2", "UGEC1", 
                    "UGEC1", "TEC_inj", "Bcells", "TEC_PT", 
                    "UGEC1", "UGEC1", "Tcells")
FSGSClusterNames <- c("UGEC1", "UGEC1", "MO_infl", "MO_kdnrs", 
                      "UGEC1", "UGEC1", "MSC", "TEC_prg", 
                      "TEC_inj", "Tcells", "MO_MT+", "TEC_inj",
                      "GRAN", "TEC_PT", "UGEC1", "TEC_CD", "UGEC2")
COVAKIClusterNames <- c("MO_infl", "UGEC1", "TEC_inj", "GRAN", 
                      "UGEC1", "UGEC1", "TEC_str", "TEC_dmg", 
                      "UGEC1", "TEC_PT", "GRAN", "MO_MT+",
                      "MO_kdnrs", "TEC_PT", "TEC_prg", "UGEC1", "TEC_CD",
                      "UGEC1", "Tcells", "PDC", "Bcells", "ERY")

## make a list of all cell types present in any of the datasets and number them ----
## (to make sure that numbering of clusters is consistent across datasets)
allclusters <- c("PDC", "MSC", "TEC_PT", "TEC_CD","TEC_inj",  "TEC_str",  "TEC_dmg", "TEC_prlf", "TEC_prg",  "MO_infl" ,
                 "MO_kdnrs", "MO_SPP1+", "MO_MT+",   "MO_dmg",  "GRAN", "Tcells",  "Bcells",      "ERY",      "UGEC1",    "UGEC2")   
names(allclusters) <- paste(1:20)

### name DM clusters ----
names(DMClusterNames) <- levels(DM)
DM <- RenameIdents(DM, DMClusterNames)
DM$celltype <- DM@active.ident
DM <- SetIdent(DM, value ="seurat_clusters")
DM$celltype <- factor(DM$celltype, levels(DM$celltype)[c(6,3,2,7,5,1,4)])
levels(DM$celltype)

DM <- SetIdent(DM, value ="celltype")
levels(DM)
DMClusterNrs  <- names(allclusters[allclusters %in% DM$celltype])
names(DMClusterNrs) <- levels(DM)
DM <- RenameIdents(DM, DMClusterNrs)
DM$celltype_nr <- DM@active.ident
DM <- SetIdent(DM, value ="seurat_clusters")
levels(DM$celltype_nr)

### name FSGS clusters ----
names(FSGSClusterNames) <- levels(FSGS)
FSGS <- RenameIdents(FSGS, FSGSClusterNames)
FSGS$celltype <- FSGS@active.ident
FSGS <- SetIdent(FSGS, value ="seurat_clusters")
FSGS$celltype <- factor(FSGS$celltype, levels(FSGS$celltype)[c(4,10,11,6,5,2,3,8,9,7,1,12)])
levels(FSGS$celltype)

FSGS <- SetIdent(FSGS, value ="celltype")
levels(FSGS)
FSGSClusterNrs  <- names(allclusters[allclusters %in% FSGS$celltype])
names(FSGSClusterNrs) <- levels(FSGS)
FSGS <- RenameIdents(FSGS, FSGSClusterNrs)
FSGS$celltype_nr <- FSGS@active.ident
FSGS <- SetIdent(FSGS, value ="seurat_clusters")
levels(FSGS$celltype_nr)

### name COVAKI clusters ----
names(COVAKIClusterNames) <- levels(COVAKI)
COVAKI <- RenameIdents(COVAKI, COVAKIClusterNames)
COVAKI$celltype <- COVAKI@active.ident
COVAKI <- SetIdent(COVAKI, value ="seurat_clusters")
COVAKI$celltype <- factor(COVAKI$celltype, levels(COVAKI$celltype)[c(13,7,11,3,5,6,10,1,9,8,4,12,14,15,2)])
levels(COVAKI$celltype)

COVAKI <- SetIdent(COVAKI, value ="celltype")
levels(COVAKI)
COVAKIClusterNrs  <- names(allclusters[allclusters %in% COVAKI$celltype])
names(COVAKIClusterNrs) <- levels(COVAKI)
COVAKI <- RenameIdents(COVAKI, COVAKIClusterNrs)
COVAKI$celltype_nr <- COVAKI@active.ident
COVAKI <- SetIdent(COVAKI, value ="seurat_clusters")
levels(COVAKI$celltype_nr)

### name CSA clusters ----
CSA <- SetIdent(CSA, value ="seurat_clusters")
names(CSAClusterNames) <- levels(CSA)
CSA <- RenameIdents(CSA, CSAClusterNames)
CSA$celltype <- CSA@active.ident
CSA <- SetIdent(CSA, value ="seurat_clusters")
levels(CSA$celltype)
CSA$celltype <- factor(CSA$celltype, levels(CSA$celltype)[c(16,14,15,5,6,12,13,4,3,7,8,11,9,10,17,1,2)])
levels(CSA$celltype)

CSA <- SetIdent(CSA, value ="celltype")
levels(CSA$celltype)
CSAClusterNrs  <- names(allclusters[allclusters %in% CSA$celltype])
names(CSAClusterNrs) <- levels(CSA)
CSA <- RenameIdents(CSA, CSAClusterNrs)
CSA$celltype_nr <- CSA@active.ident
CSA <- SetIdent(CSA, value ="seurat_clusters")
levels(CSA$celltype_nr)


### name SEP clusters ----
names(SEPClusterNames) <- levels(SEP)
SEP <- RenameIdents(SEP, SEPClusterNames)
SEP$celltype <- SEP@active.ident
SEP <- SetIdent(SEP, value ="seurat_clusters")
levels(SEP$celltype)
SEP$celltype <- factor(SEP$celltype, levels(SEP$celltype)[c(9,16,14,6,3,15,5,2,10,8,11,13,4,12,1,7)])
levels(SEP$celltype)

SEP <- SetIdent(SEP, value ="celltype")
levels(SEP)
SEPClusterNrs  <- names(allclusters[allclusters %in% SEP$celltype])
names(SEPClusterNrs) <- levels(SEP)
SEP <- RenameIdents(SEP, SEPClusterNrs)
SEP$celltype_nr <- SEP@active.ident
SEP <- SetIdent(SEP, value ="seurat_clusters")
levels(SEP$celltype_nr)


### name PRE clusters ----
names(PREClusterNames) <- levels(PRE)
PRE <- RenameIdents(PRE, PREClusterNames)
PRE$celltype <- PRE@active.ident
PRE <- SetIdent(PRE, value ="seurat_clusters")
levels(PRE$celltype)
PRE$celltype <- factor(PRE$celltype, levels(PRE$celltype)[c(6,10,3,1,8,5,4,7,11,2,9)])
levels(PRE$celltype)

PRE <- SetIdent(PRE, value ="celltype")
levels(PRE)
PREClusterNrs  <- names(allclusters[allclusters %in% PRE$celltype])
names(PREClusterNrs) <- levels(PRE)
PRE <- RenameIdents(PRE, PREClusterNrs)
PRE$celltype_nr <- PRE@active.ident
PRE <- SetIdent(PRE, value ="seurat_clusters")
levels(PRE$celltype_nr)

# integrate all kidney cells from each dataset to SO "TEC" ----

## label disease type in metadata ----
CSA$disease <- "CS AKI"
SEP$disease <- "septic AKI"
PRE$disease <- "prerenal AKI"
DM$disease <- "DN"
FSGS$disease <- "FSGS"
COVAKI$disease <- "AKI/COVID"

## subset to only kidney parenchymal clusters ----
DM <- SetIdent(DM, value = "celltype")
FSGS <- SetIdent(FSGS, value = "celltype")
CSA <- SetIdent(CSA, value = "celltype")
SEP <- SetIdent(SEP, value = "celltype")
PRE <- SetIdent(PRE, value = "celltype")
COVAKI <- SetIdent(COVAKI, value = "celltype")
DM_TEC <- subset(DM, idents= c("TEC_PT", "TEC_inj"))
FSGS_TEC <- subset(FSGS, idents= c("TEC_PT", "TEC_CD", "TEC_inj", "TEC_prg", "MSC"))
CSA_TEC <- subset(CSA, idents= c("TEC_prg", "TEC_inj", "TEC_prlf", "TEC_PT", "TEC_CD", "PDC"))
SEP_TEC <- subset(SEP, idents= c("TEC_str", "TEC_prg","TEC_inj", "PDC", "TEC_CD", "TEC_prlf", "TEC_PT"))
PRE_TEC <- subset(PRE, idents=c("TEC_inj", "TEC_prg")) 
COVAKI_TEC <- subset(COVAKI, idents= c("TEC_str", "TEC_prg","TEC_inj", "PDC", "TEC_CD", "TEC_dmg", "TEC_PT"))

## merge kidney parenchymal subsets of above datasets ----
TEC <- purrr::reduce(list(CSA_TEC, SEP_TEC, PRE_TEC, DM_TEC, FSGS_TEC, COVAKI_TEC), merge)

## batch correction, dim red. and clustering ----
TEC@active.assay <- "RNA"
TEC <- SCTransform(object = TEC, vars.to.regress = "percent.mt", method = "glmGamPoi", verbose = T)
TEC <- RunPCA(object = TEC, dims = 1:30, verbose = T)
#perform batch correction for disease type with harmony
TEC <- RunHarmony(TEC, group.by.vars = "disease", assay.use = "SCT")
#perform dimensionality reduction and find clusters
TEC <- RunUMAP(object = TEC, reduction = "harmony", dims = 1:30, verbose = T)
TEC <- RunTSNE(object = TEC, reduction = "harmony", dims = 1:30, verbose = T)
TEC <- FindNeighbors(object = TEC, reduction = "harmony", dims = 1:30, verbose = T)
TEC <- FindClusters(object = TEC, resolution = 0.25, verbose = T)

DimPlot(TEC, label =T, split.by = "disease")
## find marker genes for clusters ----
TEC.markers <- FindAllMarkers(TEC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

markers <- TEC.markers %>% 
  group_by(cluster) %>% 
  top_n(2, avg_log2FC)
DotPlot(TEC, features = c(unique(markers$gene), "EPCAM", "CRYAB", "PSCA"))+ coord_flip()
FeaturePlot(H_AKI, )
## annotate clusters ----
# annotation of clusters was done analogous to RENAL dataset / Fig. 2 
TECClusternames <- c("TEC_MMP7", "TEC_prg", "TEC_TXNRD1", "PT", "TAL", "TEC_MMP7", "UGEC", "TEC_dmg", "TEC_MT", "PT", "CD-PC", "PDC", "PT", "TEC_prlf", "LEUK", "CD-IC", "LEUK")

# name and number TEC clusters
names(TECClusternames) <- levels(TEC)
TEC <- RenameIdents(TEC, TECClusternames)
TEC$renalcelltype <- TEC@active.ident
TEC <- SetIdent(TEC, value ="seurat_clusters")
TEC$renalcelltype <- factor(TEC$renalcelltype, levels(TEC$renalcelltype)[c(10,4,5,9,13,1,8,3,7,11,2,12,6)])

TEC <- SetIdent(TEC, value ="renalcelltype")
levels(TEC)
TECClusterNrs  <- paste(1:13)
names(TECClusterNrs) <- levels(TEC)
TEC <- RenameIdents(TEC, TECClusterNrs)
TEC$celltype_nr <- TEC@active.ident
TEC <- SetIdent(TEC, value ="seurat_clusters")
levels(TEC$celltype_nr)


TEC <- SetIdent(TEC, value = "renalcelltype")
TEC <- subset(TEC, ident = c("UGEC", "LEUK"), invert =T)


# Suppl. Fig 16A-C, Suppl. Fig. 18A-F ----

names(celltype_nr_colors) <- names(allclusters)

a <- DimPlot(CSA, group.by = "celltype_nr", label=T, repel=T, label.size = 5) + ggtitle("cardiac surgery assoc. AKI\n(20249 cells)")+scale_color_manual(values= celltype_nr_colors)+
  scale_color_manual(values = celltype_nr_colors, labels = paste(names(allclusters), allclusters))
b <- DimPlot(SEP, group.by = "celltype_nr", label=T, repel=T, label.size = 5) + ggtitle("pneumonia/AKI\n(15582 cells)")+ scale_color_manual(values= celltype_nr_colors) +  NoLegend()
c <- DimPlot(PRE, group.by = "celltype_nr", label=T, repel=T, label.size = 5) + ggtitle("prerenal AKI\n(6777 cells)")+ scale_color_manual(values= celltype_nr_colors) + NoLegend()
d <- DimPlot(DM, group.by = "celltype_nr", label=T, repel=T, label.size = 5) + ggtitle("Diabetic nephropathy\n(22161 cells, Abedini et al.)")+ scale_color_manual(values= celltype_nr_colors) + NoLegend()
e <- DimPlot(FSGS, group.by = "celltype_nr", reduction= "umap", label=T, repel=T, label.size = 5) + ggtitle("FSGS\n(7606 cells, Latt et al.)")+ NoLegend()+ scale_color_manual(values= celltype_nr_colors) + NoLegend()
f <- DimPlot(COVAKI, group.by = "celltype_nr", label=T, repel=T, label.size = 5) + ggtitle("2. AKI cohort\n(20065 cells, Cheung et al.)")+ NoLegend()+ scale_color_manual(values= celltype_nr_colors) + NoLegend()


## compile and save ----
pw1 <- a+b+c+plot_layout(widths=c(1,1,1), heights=c(1), guides = "collect")
pw2 <- a+b+c+d+e+f+plot_layout(widths=c(1,1,1), heights=c(1,1), guides = "collect")

ggsave(plot = pw1, "supplfig16abc.png", h = 1300, w = 4200, units = "px", scale = 1, type = "cairo-png")
ggsave(plot = pw2, "supplfig18abcdef.png", h = 2600, w = 4200, units = "px", scale = 1, type = "cairo-png")


# Suppl. Fig. 16D ----

comparison_plots3 <- function(SO, xa, xn, ygroup, ya) {
  groups <- levels(factor(SO@meta.data[[xn]]))
  my_comparisons <- list(groups[1:2]  , groups[c(1,3)],groups[2:3])
  a <-SO@meta.data %>% group_by({{ygroup}}, patient, {{xa}}) %>%
    dplyr::summarize(count = n()) %>%
    spread({{ygroup}}, count, fill = 0) %>%
    ungroup() %>%
    filter({{xa}} != "pool") %>%
    ggplot(aes(x={{xa}}, y={{ya}}, fill={{xa}})) +
    geom_boxplot() +
    geom_jitter(size = 1, width=0.25, height=0.25)+
    scale_y_log10(name= "abs. cell count")+
    scale_x_discrete(limits=c("CS", "pneumonia", "prerenal"), labels=c("CS", "pneumonia", "prerenal"))+
    scale_fill_manual(name = 'AKI type', values = AKI_type_colors)+
    theme_classic()+
    theme(axis.title.x =element_blank(), 
          axis.text.x = element_text(angle=45, hjust=1, size=10, face="bold"),
          axis.text.y = element_text(size=10, face="bold"),
          title = element_text(size =12, face="bold", vjust=0.5))+
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired =F, label="p.signif", hide.ns=F)+
    ggtitle(deparse(substitute(ya)))
  
  return(a)
  
}

pw5 <- comparison_plots3(URINE, AKI_type, "AKI_type", cellgroup, TEC)+
  comparison_plots3(URINE, AKI_type, "AKI_type", cellgroup, PDC)+
  comparison_plots3(URINE, AKI_type, "AKI_type", cellgroup, LEUK)+
  comparison_plots3(URINE, AKI_type, "AKI_type", cellgroup, UGEC)+
  comparison_plots3(URINE, AKI_type, "AKI_type", celltype, MO_infl)+
  comparison_plots3(URINE, AKI_type, "AKI_type", celltype, MO_kdnrs)+
  comparison_plots3(URINE, AKI_type, "AKI_type", celltype, Tcells)+
  comparison_plots3(URINE, AKI_type, "AKI_type", celltype, Bcells)+
  comparison_plots3(URINE, AKI_type, "AKI_type", celltype, TEC_inj)+
  comparison_plots3(URINE, AKI_type, "AKI_type", celltype, TEC_dmg1)+
  comparison_plots3(URINE, AKI_type, "AKI_type", celltype, TEC_prlf)+
  comparison_plots3(URINE, AKI_type, "AKI_type", celltype, TEC_prg)+
  plot_layout(guides = "collect")

for (i in c(2:4,6:8,10:12)) { pw5[[i]] <- pw5[[i]] + 
  theme(axis.title.y = element_blank())}
for (i in c(1:8)) { pw5[[i]] <- pw5[[i]] + 
  theme(axis.text.x = element_blank())}
ggsave(plot = pw5, "supplfig16d.png", h = 4000, w = 4200, units = "px", scale = 1, type = "cairo-png")


# Suppl. Fig. 18 G ----
SOlist <- list(CSA, SEP, PRE, COVAKI, DM, FSGS)
names <- levels(factor(TEC$disease))[c(2,6,5,1,3,4)]

for (i in 1:6) {
  SOlist[[i]]$ct <- NA
  SOlist[[i]]@meta.data[SOlist[[i]]$celltype_nr %in% c(3:9),]$ct <- "TEC"
  SOlist[[i]]@meta.data[SOlist[[i]]$celltype_nr %in% c(16:17),]$ct <- "LYMPH"
  SOlist[[i]]@meta.data[SOlist[[i]]$celltype_nr %in% c(10:15),]$ct <- "MYEL"
  SOlist[[i]]@meta.data[SOlist[[i]]$celltype_nr %in% c(19:20),]$ct <- "UGEC"
}

tbl <- as.data.frame(table(tibble(SOlist[[1]]$orig.ident, SOlist[[1]]$ct)))
tbl$disease <- levels(factor(TEC$disease))[c(2)]
colnames(tbl) <- c("orig.ident", "ct", "Freq", "disease")

for (i in 1:length(SOlist)) {
x <- as.data.frame(table(tibble(SOlist[[i]]$orig.ident, SOlist[[i]]$ct)))
x$disease <- names[i]
colnames(x) <- c("orig.ident", "ct", "Freq", "disease")
tbl <- rbind(tbl,x)
}
tbl <- tbl[tbl$orig.ident != "POOL1" & tbl$orig.ident != "POOL2",]
tbl <- tbl[!duplicated(tbl),]

groups <- names
my_comparisons <- list(groups[c(1,2)], groups[c(1,3)], groups[c(1,4)], groups[c(1,5)], groups[c(1,6)], 
                       groups[c(2,3)], groups[c(2,4)], groups[c(2,5)], groups[c(2,6)],
                       groups[c(3,4)], groups[c(3,5)], groups[c(3,6)])

compPlot <- function(x) {tbl %>% filter(ct == x) %>%
ggplot(aes(x=disease, y=Freq, fill=disease)) +
  geom_boxplot() +
  geom_jitter(size = 1, width=0.25, height=0.25)+
  scale_x_discrete(limits = c("CS AKI", "septic AKI", "prerenal AKI", "AKI/COVID", "DN", "FSGS"),
                   labels = c("CS AKI", "pneumonia", "prerenal AKI", "2. AKI cohort", "DN", "FSGS"))+
  scale_fill_manual(values = diseasecolors, limits = c("CS AKI", "septic AKI", "prerenal AKI", "AKI/COVID", "DN", "FSGS"),
                    labels = c("CS AKI", "pneumonia", "prerenal AKI", "2. AKI cohort", "DN", "FSGS"))+
  scale_y_log10(name= "abs. cell count")+
  theme_classic()+
  theme(axis.title.x =element_blank(), 
        axis.text.x = element_text(angle=45, hjust=1, size=10, face="bold"),
        axis.text.y = element_text(size=10, face="bold"),
        title = element_text(size =12, face="bold", vjust=0.5))+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired =F, label="p.signif", hide.ns=F)+
  ggtitle(x)
}

pw <- compPlot("TEC")+
compPlot("MYEL")+
compPlot("LYMPH")+
compPlot("UGEC")+plot_layout(widths=c(1,1,1,1), height = c(1),guides="collect")
  
ggsave(plot = pw, "supplfig18g.png", h = 2000, w = 4200, units = "px", scale = 1, type = "cairo-png")


# Suppl. Fig 19 ----

## plotting ----

# split TEC dataset by disease 
TEC <- SetIdent(TEC, value = "disease")

TEC_sub_FSGS <- subset(TEC, idents=c("FSGS"))
TEC_sub_DN <- subset(TEC, idents=c("DN"))
TEC_sub_CSA <- subset(TEC, idents=c("CS AKI"))
TEC_sub_SEP <- subset(TEC, idents=c("septic AKI"))
TEC_sub_PRE <- subset(TEC, idents=c("prerenal AKI"))
TEC_sub_COVAKI <- subset(TEC, idents=c("AKI/COVID"))

# plots
h <- DimPlot(TEC_sub_CSA, pt.size=1, label=T, repel=T, label.size = 5, group.by = "celltype_nr") +
  ggtitle(paste("cardiac surgery assoc. AKI\n", ncol(TEC_sub_CSA), "renal cells (33%)"))
i <- DimPlot(TEC_sub_SEP, pt.size=1, label=T, repel=T, label.size = 5, group.by = "celltype_nr") +
  ggtitle(paste("pneumonia AKI\n", ncol(TEC_sub_SEP), "renal cells (36%)"))
j <- DimPlot(TEC_sub_PRE, pt.size=1, label=T, repel=T, label.size = 5, group.by = "celltype_nr") +
  ggtitle( paste("prerenal AKI\n",ncol(TEC_sub_PRE), "renal cells (7%)"))
k <- DimPlot(TEC_sub_DN, pt.size=1, label=T, repel=T, label.size = 5, group.by = "celltype_nr")  +
  ggtitle(paste("Diabetic nephropathy\n",ncol(TEC_sub_DN), "renal cells (9%)"))
l <- DimPlot(TEC_sub_FSGS, pt.size=1, label=T, repel=T, label.size = 5, group.by = "celltype_nr")+
  ggtitle(paste("FSGS\n",ncol(TEC_sub_FSGS), "renal cells (17%)"))
m <- DimPlot(TEC_sub_COVAKI, pt.size=1, label=T, repel=T, label.size = 5, group.by = "celltype_nr") +
  ggtitle(paste("2. AKI cohort\n", ncol(TEC_sub_COVAKI), "renal cells (31%)"))


### compile and save ----

pw3 <- h+i+j+k+l+m+plot_layout(widths=c(1,1,1), height = c(1,1),guides="collect")&
  scale_color_manual(values=disease_celltype_colors, labels =paste(levels(TEC_sub_CSA$celltype_nr), levels(TEC_sub_CSA$renalcelltype)), limits = levels(TEC_sub_CSA$celltype_nr)[1:11])

ggsave(plot = pw3, "supplfig19.png", h = 2600, w = 4200, units = "px", scale = 1, type = "cairo-png")



# Suppl. Fig. 20 ----
AKI_type_colors <- paletteer_d("unikn::pal_signal")
diseasecolors <- c(AKI_type_colors, "lightblue", "pink", "#c98ed4")
names(diseasecolors) <- c("CS AKI", "prerenal AKI", "septic AKI", "FSGS", "DN", "AKI/COVID")
TEC$disease <- factor(TEC$disease, levels(factor(TEC@meta.data$disease))[c(2,6,5,1,3,4)])

comparison_plots_disease <- function(SO, xa, xn, ygroup, ya) {
  groups <- levels(factor(SO@meta.data[[xn]]))
  my_comparisons <- list(groups[c(1,2)], groups[c(1,3)], groups[c(1,4)], groups[c(1,5)], groups[c(1,6)], 
                         groups[c(2,3)], groups[c(2,4)], groups[c(2,5)], groups[c(2,6)],
                         groups[c(3,4)], groups[c(3,5)], groups[c(3,6)])
  a <-SO@meta.data %>% group_by({{ygroup}}, orig.ident, {{xa}}) %>%
    dplyr::summarize(count = n()) %>%
    spread({{ygroup}}, count, fill = 0) %>%
    ungroup() %>%
    filter({{xa}} != "pool") %>%
    ggplot(aes(x={{xa}}, y={{ya}}, fill={{xa}})) +
    geom_boxplot() +
    geom_jitter(size = 1, width=0.25, height=0.25)+
    scale_y_log10(name= "abs. cell count")+
    #scale_x_discrete(limits=c("CS", "pneumonia", "prerenal"), labels=c("CS", "pneumonia", "prerenal"))+
    scale_fill_manual(name = 'disease type', values = diseasecolors)+
    theme_classic()+
    theme(axis.title.x =element_blank(), 
          axis.text.x = element_text(angle=45, hjust=1, size=10, face="bold"),
          axis.text.y = element_text(size=10, face="bold"),
          title = element_text(size =12, face="bold", vjust=0.5))+
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired =F, label="p.signif", hide.ns=F)+
    ggtitle(deparse(substitute(ya)))
  
  return(a)
}


pw4 <- comparison_plots_disease(TEC, disease, "disease", renalcelltype, PDC)+
  comparison_plots_disease(TEC, disease, "disease", renalcelltype, PT)+
  comparison_plots_disease(TEC, disease, "disease", renalcelltype, TAL)+
  comparison_plots_disease(TEC, disease, "disease", renalcelltype, `CD-PC`)+
  comparison_plots_disease(TEC, disease, "disease", renalcelltype, `CD-IC`)+
  comparison_plots_disease(TEC, disease, "disease", renalcelltype, TEC_MMP7)+
  comparison_plots_disease(TEC, disease, "disease", renalcelltype, TEC_MT)+
  comparison_plots_disease(TEC, disease, "disease", renalcelltype, TEC_TXNRD1)+
  comparison_plots_disease(TEC, disease, "disease", renalcelltype, TEC_dmg)+
  comparison_plots_disease(TEC, disease, "disease", renalcelltype, TEC_prlf)+
  comparison_plots_disease(TEC, disease, "disease", renalcelltype, TEC_prg)+
  plot_layout(guides = "collect")


for (i in c(2:4,6:8,10:11)) { pw4[[i]] <- pw4[[i]] + 
  theme(axis.title.y = element_blank())}
for (i in c(1:7)) { pw4[[i]] <- pw4[[i]] + 
  theme(axis.text.x = element_blank())}
ggsave(plot = pw4, "supplfig20.png", h = 4000, w = 4200, units = "px", scale = 1, type = "cairo-png")


