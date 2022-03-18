# description ----

# This is the R script for generating Figure 4 in 
# "Urinary single-cell sequencing captures intrarenal injury and repair processes 
# in human acute kidney injury" by Klocke et al. 

# in here, we download the proximal tubular cell ischemia reperfusion dataset by Kirita et al. and use it as 
# reference for our own urinary tubular epithelial cells. 
# Kirita Y, Wu H, Uchimura K, Wilson PC, Humphreys BD. Cell profiling of mouse acute kidney injury reveals conserved cellular responses to injury. Proc Natl Acad Sci. 2020;117(27):15874-15883. doi:10.1073/pnas.2005477117

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
library(biomaRt)

## loading own data ----
RENAL <-readRDS("~/SO_kidney_urine_cells.rds")

## loading color schemes ----
renalcelltype_colors <- 
  c(paletteer_d("ggsci::purple_material")[c(3,5)],                  #PDC colors
    paletteer_d("ggsci::teal_material")[c(8,6,4,3)],                #healthy TEC PT and LOH colors   
    paletteer_d("ggsci::light_blue_material")[c(3,5,7,9)],          #healthy TEC DT and CD colors
    paletteer_d("ggsci::amber_material")[c(3,4,6,8,10)],            #EMT/infl TEC
    paletteer_d("ggsci::red_material")[c(2,5,8)],                   #OxyStress TEC
    paletteer_d("ggsci::pink_material")[4],                         #prolif TEC
    paletteer_d("ggsci::light_green_material")[c(9,7,5,3)])         #progenitors 
PTgroupcolors <- paletteer_d("fishualize::Scarus_hoefleri")[c(4,3,5,1)]

## loading reference dataset GSE139107 ----
M_IRI_metadata <- read.delim("~/GSE139107_MouseIRI/GSE139107_MouseIRI.metadata.txt", stringsAsFactors = FALSE)
M_IRI_4h.data <- read.delim("~/GSE139107_MouseIRI/GSE139107_MouseIRI_4hours.dge.txt", stringsAsFactors = FALSE)
M_IRI_12h.data <- read.delim("~/GSE139107_MouseIRI/GSE139107_MouseIRI_12hours.dge.txt", stringsAsFactors = FALSE)
M_IRI_2d.data <- read.delim("~/GSE139107_MouseIRI/GSE139107_MouseIRI_2days.dge.txt", stringsAsFactors = FALSE)
M_IRI_14d.data <- read.delim("~/GSE139107_MouseIRI/GSE139107_MouseIRI_14days.dge.txt", stringsAsFactors = FALSE)
M_IRI_6w.data <- read.delim("~/GSE139107_MouseIRI/GSE139107_MouseIRI_6weeks.dge.txt", stringsAsFactors = FALSE)
M_IRI_c.data <- read.delim("~/GSE139107_MouseIRI/GSE139107_MouseIRI_control.dge.txt", stringsAsFactors = FALSE)

# create, subset and merge Seurat objects ----

## Initialize Seurat objects ----
M_IRI_4h <- CreateSeuratObject(counts = M_IRI_4h.data, project = "M_IRI_PT", min.cells = 3, min.features = 200)
M_IRI_12h <- CreateSeuratObject(counts = M_IRI_12h.data, project = "M_IRI_PT", min.cells = 3, min.features = 200)
M_IRI_2d <- CreateSeuratObject(counts = M_IRI_2d.data, project = "M_IRI_PT", min.cells = 3, min.features = 200)
M_IRI_14d <- CreateSeuratObject(counts = M_IRI_14d.data, project = "M_IRI_PT", min.cells = 3, min.features = 200)
M_IRI_6w <- CreateSeuratObject(counts = M_IRI_6w.data, project = "M_IRI_PT", min.cells = 3, min.features = 200)
M_IRI_c <- CreateSeuratObject(counts = M_IRI_c.data, project = "M_IRI_PT", min.cells = 3, min.features = 200)

## add metadata from M_IRI_metadata to each object
M_IRI.list <- list(M_IRI_4h, M_IRI_12h, M_IRI_2d, M_IRI_14d, M_IRI_6w, M_IRI_c)

groups <- c("4hours", "12hours", "2days", "14days", "6weeks", "Control")
for (i in 1:length(M_IRI.list)) {
  meta <- filter(M_IRI_metadata, Group == groups[i])
  celltype <- meta$celltype
  names(celltype) <- rownames(meta)
  group <- meta$Group
  names(group) <- rownames(meta)
  replicates <- meta$Replicates
  names(replicates) <- rownames(meta)
  M_IRI.list[[i]] <- AddMetaData(M_IRI.list[[i]], celltype, col.name = "celltype")
  M_IRI.list[[i]] <- AddMetaData(M_IRI.list[[i]], replicates, col.name = "replicates")
  M_IRI.list[[i]] <- AddMetaData(M_IRI.list[[i]], group, col.name = "group")
}

# normalization, dim. reduction and clustering ----

# proceed with preparation of PT cells like described in Kritia paper: "We extracted mouse proximal
# tubular cell clusters and then performed clustering without Harmony integration.
# The highly variable genes for principal component analysis were
# obtained by identifying the top 300 variable genes from each dataset with
# FindVariableFeatures and merging the list. We then performed principal
# component analysis ("RunPCA" function), clustering, and UMAP.

## normalize, subset each object to proximal tubular cells, find variable features ----
for (i in 1:length(M_IRI.list)) {
  M_IRI.list[[i]]$celltype <- factor(M_IRI.list[[i]]$celltype)
  M_IRI.list[[i]]@active.ident <- as.factor(M_IRI.list[[i]]$celltype)
  M_IRI.list[[i]] <- NormalizeData(M_IRI.list[[i]])
  M_IRI.list[[i]] <- subset(M_IRI.list[[i]], idents = c("PTS1", "PTS2", "PTS3", "NewPT1", "NewPT2"))
  M_IRI.list[[i]] <- FindVariableFeatures(M_IRI.list[[i]], nfeatures = 300)
}
var.features <- c()

for (i in 1:6) { 
  var.features <- c(var.features, M_IRI.list[[i]]@assays[["RNA"]]@var.features)
}

var.features <- unique(var.features)


## merge objects ----
M_IRI_PT <- purrr::reduce(M_IRI.list, merge)
### add variable features 
M_IRI_PT@assays[["RNA"]]@var.features <- var.features

## PCA and scaling ----
M_IRI_PT <- ScaleData(object = M_IRI_PT, verbose = T)
M_IRI_PT <- RunPCA(object = M_IRI_PT, dims = 1:15, verbose = T)

## dimensionality reduction and clustering ----
M_IRI_PT <- RunUMAP(object = M_IRI_PT, reduction = "pca", dims = 1:15, verbose = T)
M_IRI_PT <- RunTSNE(object = M_IRI_PT, reduction = "pca", dims = 1:15, verbose = T)
M_IRI_PT <- FindNeighbors(object = M_IRI_PT, reduction = "pca", dims = 1:15, verbose = T)
M_IRI_PT <- FindClusters(object = M_IRI_PT, resolution = 0.1, verbose = T)

# Annotation of PT cell clusters ----
## Marker gene (Kritia et al Fig 2a) expression ----
DotPlot(M_IRI_PT, features = rev(c("Slc7a13", "Slc22a30", "Myc", "Havcr1", "Hspa1a", "Krt20","Top2a", "Sema5a", "Dcdc2a",  "Vcam1")), col.min = 0, dot.scale = 5) + 
  coord_flip()

## annotation of clusters to "M_IRI_PT$PTgroup" ----
clusternames <- c("Healthy", "Injured", "Failed repair", "Healthy", "Repairing", "Healthy", "Injured", "Healthy")
names(clusternames) <- levels(M_IRI_PT)
M_IRI_PT <- RenameIdents(M_IRI_PT, clusternames)
M_IRI_PT$PTgroup <- M_IRI_PT@active.ident
M_IRI_PT <- SetIdent(M_IRI_PT, value ="seurat_clusters")
M_IRI_PT$PTgroup <- factor(M_IRI_PT$PTgroup, levels(M_IRI_PT$PTgroup)[c(1,2,4,3)])

# Fig. 4A, reference dataset ----
## Dimplot ----
a <- DimPlot(M_IRI_PT, reduction = "umap", group.by = "PTgroup", label = T, label.size = 6, cols =PTgroupcolors) + theme_void() + 
  NoLegend() + ggtitle("Murine prox. tubule ischemia\nreperfusion injury (Kirita et al.)") +
  theme(plot.title = element_text(face= "bold", size = 20, hjust=0.5, vjust = 1))

## Dimplots split by timepoint ----
b <- DimPlot(M_IRI_PT, reduction = "umap", group.by = "PTgroup", split.by = "group", label = F, cols = PTgroupcolors, ncol=3) &NoLegend() &
  theme(plot.title = element_blank(), 
        axis.line = element_blank(),
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(), text = element_text(size=12))

## Dotplot marker genes ----
c <- DotPlot(M_IRI_PT, features = rev(c("Slc7a13", "Slc22a30", "Havcr1", "Krt20","Top2a", "Sema5a", "Dcdc2a")), group.by = "PTgroup",  cols = c("lightgrey", "blue"), col.min = 0, dot.scale = 5) + 
  theme_ipsum() +
  theme(axis.text.y = element_text(size =12, color = PTgroupcolors, face = "bold", hjust = 1), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12, face = "bold", angle = 90, hjust= 1, vjust=0.5),
        legend.title = element_text(size = 12),
        
        legend.text = element_text(size = 12),
        legend.position='right',
        legend.direction = "vertical",
        legend.box = "vertical") +
  scale_colour_viridis(option = "viridis", direction = -1)

## compile and save ----

layout <- "
AAAACCC
AAAABBB
AAAABBB"

ggsave(plot = a+b+c + plot_layout(design = layout, guides="collect"), "fig4a.png", w = 4200, h = 2000, units = "px", scale = 1, type = "cairo-png")

# Fig. 4B automatic annotation of urine cells with mouse PT as reference ----

# Here, we use singleR to automatically annotate our urine cells and use the mosue PT cells as reference. 
# We use the workflow for single-cell references described here:
# https://bioconductor.org/packages/devel/bioc/vignettes/SingleR/inst/doc/SingleR.html#1_Introduction

## Create human orthologous gene list from mouse data ----
ensembl = useEnsembl(biomart="ensembl")
view(listDatasets(ensembl))
ensembl.human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl.mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
ortholog.genes <- getLDS(attributes = c("mgi_symbol", "chromosome_name"),
                         filters = "mgi_symbol", values = M_IRI_PT@assays[["RNA"]]@data@Dimnames[[1]], mart = ensembl.mouse,
                         attributesL = c("hgnc_symbol"), martL = ensembl.human)

for (i in 1:length(ortholog.genes$MGI.symbol)) {
  M_IRI_PT@assays[["RNA"]]@data@Dimnames[[1]][M_IRI_PT@assays[["RNA"]]@data@Dimnames[[1]] == ortholog.genes$MGI.symbol[i]] <- ortholog.genes$HGNC.symbol[i]
}

## retrieve data, lognormalize expression matrices for singleR
test <- GetAssayData(RENAL, assay = 'RNA', slot = 'data')
test <- LogNormalize(test)
ref <- GetAssayData(M_IRI_PT, assay = "RNA", slot = "data")
ref <- LogNormalize(ref)

## use singleR to compare mouse data with our REANL SO ----
singler_results_kidney_ref <- SingleR::SingleR(
  test = test,
  ref = ref,
  labels = M_IRI_PT$PTgroup, de.method="wilcox"
)

## annotate urine TEC in RENAL SO via new metadata $ref_labels ----
RENAL@meta.data$ref_labels <- singler_results_kidney_ref@listData$labels

## visualize in heatmap ----
names(renalcelltype_colors) <- levels(RENAL$renalcelltype)
names(PTgroupcolors) <- levels(M_IRI_PT$PTgroup)
annot.colors <- list(mouse_PT = PTgroupcolors, urine = renalcelltype_colors)
annot.col <- data.frame(
  mouse_PT = RENAL$ref_labels,
  urine = RENAL$renalcelltype)

p <- SingleR::plotScoreHeatmap(
  singler_results_kidney_ref,
  main = "Urine TEC annotation by mouse ischemia reference",
  show.labels = F,
  annotation_col = annot.col, 
  annotation_colors = annot.colors, 
  annotation_legend = F,
  fontsize = 16
)

## save plot ----
ggsave(plot = p, "fig4b.png", w = 3300, h = 1000, units = "px", scale = 1, type = "cairo-png")


# Fig. 4C Reference labels in DimPlot ----

## find cells for each label ----
failed_repair <- WhichCells(RENAL, expression = ref_labels == "Failed repair")
repairing <- WhichCells(RENAL, expression = ref_labels == "Repairing")
healthy <- WhichCells(RENAL, expression = ref_labels == "Healthy")
injury <- WhichCells(RENAL, expression = ref_labels == "Injured")

## plot in UMAP ----
d <- DimPlot(RENAL, cells.highlight = healthy, cols.highlight = PTgroupcolors[1], pt.size =1) + 
  ggtitle("Healthy") + 
  theme_void() + NoLegend() + 
  theme(plot.title = element_text(face= "bold", size = 20, hjust=0.5))

e <- DimPlot(RENAL, cells.highlight = repairing, cols.highlight = PTgroupcolors[3], pt.size =1) + 
  ggtitle("Repairing") + theme_void() + NoLegend() + 
  theme(plot.title = element_text(face= "bold", size = 20, hjust=0.5))

f <- DimPlot(RENAL, cells.highlight = failed_repair, cols.highlight = PTgroupcolors[4], pt.size =1) + 
  ggtitle("Failed repair") + theme_void() + NoLegend() + 
  theme(plot.title = element_text(face= "bold", size = 20, hjust=0.5))

h <- DimPlot(RENAL, cells.highlight = injury, cols.highlight = PTgroupcolors[1], sizes.highlight = 1, pt.size =1) + 
  ggtitle("Injured") + 
  theme_void() + NoLegend() + 
  theme(plot.title = element_text(face= "bold", size = 30, hjust=0.5))

## compile and save ----
ggsave(plot = d + e + f, "fig4c.png", w = 3500, h = 1200, units = "px", scale = 1, type = "cairo-png")

# Fig. 4D Barplot ----

## dataframe ----
StatePerType.df <- as.data.frame(table(tibble(celltype= RENAL$renalcellgroup, state= RENAL$ref_labels)))
StatePerType.df <- StatePerType.df[StatePerType.df$celltype != "PDC",]
StatePerType.df <- filter(StatePerType.df, state != "Injured")
## plot ----
g <- ggplot(StatePerType.df, aes(celltype, Freq)) +
  geom_bar(aes(fill = state), position = 'stack', stat = 'identity') +
  scale_fill_manual(values = PTgroupcolors) +
  coord_cartesian(clip = 'off') +
  theme_void() +
  theme(
    legend.position = 'top',
    legend.direction = "vertical",
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(vjust = 0.5, hjust=1, face= "bold", size=14, angle=90),
    axis.text.y = element_text(vjust = 0.5, hjust=1, face= "bold", size=14), 
    legend.title = element_text(vjust = 0.5, hjust=0, face= "bold", size=14),
    legend.text = element_text(vjust = 0.5, hjust=0, size=14)
  ) 
## save ----
ggsave(plot = g, "fig4d.png", w = 600, h = 2200, units = "px", scale = 1, type = "cairo-png")
