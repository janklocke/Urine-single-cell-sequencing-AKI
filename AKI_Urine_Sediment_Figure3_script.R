# description ----

# This is the R script for generating Figure 3 and Supplemental Figure 11-12 in 
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
library(rlist)
library(data.table)
library(rlang)
library(paletteer)
library(corrplot)
library("Hmisc")
library(lattice)
library("PerformanceAnalytics")
library(ggcorrplot2)
library(hrbrthemes)

## loading data ----
URINE <- readRDS("~/SO_all_urine_cells.rds")
RENAL <-readRDS("~/SO_kidney_urine_cells.rds")

## the following seurat object contains post-mortem acute kidney injury kidney biopsy (AKIKB) samples and 
## was provided by Christian Hinze from his recently published preprint "Transcriptomic 
## responses of the human kidney to acute injury at single cell resolution". 
## the data will become publicly available upon publication of the peer reviewed version. 
## https://doi.org/10.1101/2021.12.15.472619 
AKIKB <-  readRDS("~/AKIKB.rds")

## loading color schemes ----
renalcelltype_colors <- 
  c(paletteer_d("ggsci::purple_material")[c(3,5)],                  #PDC colors
    paletteer_d("ggsci::teal_material")[c(8,6,4,3)],                #healthy TEC PT and LOH colors   
    paletteer_d("ggsci::light_blue_material")[c(3,5,7,9)],          #healthy TEC DT and CD colors
    paletteer_d("ggsci::amber_material")[c(3,4,6,8,10)],            #EMT/infl TEC
    paletteer_d("ggsci::red_material")[c(2,5,8)],                   #OxyStress TEC
    paletteer_d("ggsci::pink_material")[4],                         #prolif TEC
    paletteer_d("ggsci::light_green_material")[c(9,7,5,3)])         #progenitors 
names(renalcelltype_colors) <- levels(RENAL$renalcelltype_nr)
celltype_colors <- colorRampPalette(paletteer_d("palettesForR::Pastels"))(19)
AKIKBcolors <- c(paletteer_d("ggsci::red_material")[c(4:10)], paletteer_d("palettesForR::Pastels")[c(2:8,10:14)], paletteer_d("lisa::JackBush_1")[c(2,5)])
JOINED_Ucolors <- c(renalcelltype_colors[c(1,3:10,12,16,19,22)], celltype_colors[c(9,10,14,15,16)]) 

locationcolors <- paletteer_d("lisa::JackBush_1")
statecolors <- paletteer_d("unikn::pal_signal")[c(3,1)]

names(locationcolors) <- levels(RENAL$location)
names(statecolors) <- c("healthy", "injury")

highlightcolors <- c(paletteer_d("palettesForR::Pastels")[9], paletteer_d("rcartocolor::SunsetDark")[6]) 

## functions ----
number_ticks <- function(n) {function(limits) pretty(limits, n)}
makeTransparent <- function(someColor, alpha=100){
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

# Fig.3AB, Integration of URINE and AKIKB datasets. ----

## preparation of URINE dataset ----

## annotating urinary origin for all cells 
URINE$tissue <- "urine"

## we exclude urogenital cells and erythrocytes from integration
URINE <- SetIdent(URINE, value = "renalcellgroup")
URINE_subset <- subset(URINE, idents = c("ERY", "UGEC"), invert = T)

## preparation of AKIKB dataset ----

## annotating tissue origin for all cells 
AKIKB$tissue <- "kidney"

## downsampling AKIKB to the same size of URINE_subset
AKIKB <- SetIdent(AKIKB, value = "tissue")
AKIKB_subset <- subset(AKIKB, downsample = nrow(URINE_subset@meta.data))

## merging urine (URINE) and kidney (AKIKB) data ----
JOINED <- merge(AKIKB_subset, URINE_subset)

## integration, dim. reduction and clustering ----
JOINED@active.assay <- "RNA"
JOINED <- SCTransform(object = JOINED, method = "glmGamPoi", verbose = T)
JOINED <- RunPCA(object = JOINED, dims = 1:45, verbose = T)
#perform batch correction for tissue with harmony
JOINED <- RunHarmony(JOINED, group.by.vars = "tissue", assay.use = "SCT")
#run dimensionality reduction and find clusters
JOINED <- RunUMAP(object = JOINED, reduction = "harmony", dims = 1:45, verbose = T)
JOINED <- RunTSNE(object = JOINED, reduction = "harmony", dims = 1:45, verbose = T)
JOINED <- FindNeighbors(object = JOINED, reduction = "harmony", dims = 1:45, verbose = T)
JOINED <- FindClusters(object = JOINED, resolution = 0.3, verbose = T)

## DimPlot visualization ----

# add names to unnamed cells in previous cluster annotations
JOINED@meta.data[JOINED$tissue == "kidney_biopsy",]$renalcellgroup <- "biopsy"
JOINED@meta.data[JOINED$tissue == "urine",]$exp.ct <- "urine"

# reorder clusternames and add cluster nrs for "$exp.ct" (kidney biopsy cell type annotation) and "$renalcellgroup" (urine cell type annotation). 
JOINED$exp.ct <- factor(JOINED$exp.ct)
JOINED$renalcellgroup <- factor(JOINED$renalcellgroup)
JOINED$exp.ct <- factor(JOINED$exp.ct, levels(JOINED$exp.ct)[c(16,21,19,9,7,5,3,14,15,10,20,18,8,6,4,2,11,12,13,17,1,22)])
JOINED$renalcellgroup <- factor(JOINED$renalcellgroup, levels(factor(JOINED$renalcellgroup))[c(12,13,8,1,14,7,6,5,4,16,19,18,17,10,11,9,15,2,3)])

JOINED <- SetIdent(JOINED, value ="exp.ct")
levels(JOINED)
JOINEDClusterNrs  <- paste(1:22)
names(JOINEDClusterNrs) <- levels(JOINED)
JOINED <- RenameIdents(JOINED, JOINEDClusterNrs)
JOINED$exp_nr <- JOINED@active.ident
JOINED <- SetIdent(JOINED, value ="seurat_clusters")
levels(JOINED$exp_nr)

JOINED <- SetIdent(JOINED, value ="renalcellgroup")
levels(JOINED)
JOINEDClusterNrs  <- paste(1:19)
names(JOINEDClusterNrs) <- levels(JOINED)
JOINED <- RenameIdents(JOINED, JOINEDClusterNrs)
JOINED$rcg_nr <- JOINED@active.ident
JOINED <- SetIdent(JOINED, value ="seurat_clusters")
levels(JOINED$rcg_nr)


# DimPlot and labels of kidney cells (urine cells transparent)
a <- DimPlot(JOINED, reduction = "umap", group.by = "exp_nr", cols = c(AKIKBcolors, makeTransparent("darkturquoise", 3)), label = T, repel = T, pt.size = 0.5, label.size = 6) + 
  theme_void() + 
  scale_color_manual(labels = paste(levels(JOINED$exp_nr), levels(JOINED$exp.ct)), values = c(AKIKBcolors, makeTransparent("darkturquoise", 3))) + 
  theme() +
  ggtitle("AKI post mortem biopsy") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 14),
        plot.title = element_text(face= "bold", size = 20, hjust=0.5))

# DimPlot and labels of urine cells (kidney cells transparent)
b <- DimPlot(JOINED, reduction = "umap", group.by= "rcg_nr", cols = c(JOINED_Ucolors, makeTransparent("deeppink2", 1)), label = T, repel = T, pt.size = 0.5, label.size = 6) + 
  theme_void() + 
  scale_color_manual(labels = paste(levels(JOINED$rcg_nr), str_replace(levels(JOINED$renalcellgroup), "DCT", "TAL/DCT")), values = c(JOINED_Ucolors, makeTransparent("deeppink2", 3))) + 
  theme() +
  ggtitle("AKI urine sediment") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 14),
        plot.title = element_text(face= "bold", size = 20, hjust=0.5))

## compile and save ----
ggsave(plot = b+a, "fig3integrated.png", h = 2400, w = 4200, units = "px", scale = 1, type = "cairo-png")

# Fig.3C, gene expression correlation plot ----

## function for Fig.3C ----
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

## preparation of URINE dataset ----
# exclude erythrocytes from analysis
URINE <- SetIdent(URINE, value = "renalcellgroup")
URINE_subset <- subset(URINE, idents = c("ERY"), invert = T)

## calculating correlation matrix ----

# calculate mean gene expression per cell type(AKIKB$exp.ct)
AKIKB_corr <- corr_df(AKIKB, AKIKB$exp.ct)
# calculate mean gene expression per cell type(URINE_subset$renalcellgroup)
URINE_corr <- corr_df(URINE_subset, URINE_subset$renalcellgroup)

# filter for only the highly variable genes in the datasets
AKIKB_corr2 <- AKIKB_corr[rownames(AKIKB_corr) %in% c(AKIKB@assays[["integrated"]]@var.features, URINE@assays[["SCT"]]@var.features),]
URINE_corr2 <- URINE_corr[rownames(URINE_corr) %in% c(AKIKB@assays[["integrated"]]@var.features, URINE@assays[["SCT"]]@var.features),]

# merge dataframes
MERGED_corr <- merge(AKIKB_corr2, URINE_corr2, by = 0, all = T)
rownames(MERGED_corr) <- MERGED_corr$Row.names
MERGED_corr$Row.names <- NULL

# use rcorr() to calculate correlation matrix
MERGED.mtx <- rcorr(as.matrix(MERGED_corr))

# subset matrix to wanted comparisons (urine cell types in columns, biopsy cell types in rows)
MERGED_subset <- MERGED.mtx$r[1:ncol(AKIKB_corr), (ncol(AKIKB_corr)+1):ncol(MERGED_corr)]

# rename rows/columns
rownames(MERGED_subset) <- paste(levels(factor(AKIKB$exp.ct)))
colnames(MERGED_subset) <- paste(levels(factor(URINE_subset$renalcellgroup)))

## visualize with corrplot, save ----
png("fig3corrplotCircle.png", width= 1800, height =1800, type="cairo-png")
corrplot(MERGED_subset, 
         method = "square", type = "full", addCoef.col = "white",
         cl.pos = "r", cl.ratio = 0.1, cl.lim = c(-0.1,1), cl.cex = 3.5,
         tl.col = "black", tl.srt = 90, tl.cex = 4.5, number.cex=3.5,
         addCoefasPercent = T, 
         order= "original", addgrid.col = F)
dev.off()

# Suppl Fig. 12A, automatic annotation of urine cells with AKI biopsy as reference ----
  
# Here, we use singleR to automatically annotate our urine cells and use the AKI biopsy tissue cells as reference. 
# We use the workflow for single-cell references described here:
# https://bioconductor.org/packages/devel/bioc/vignettes/SingleR/inst/doc/SingleR.html#1_Introduction
  
## retrieve data, lognormalize expression matrices for singleR
test <- GetAssayData(RENAL, assay = 'RNA', slot = 'data')
test <- LogNormalize(test)
ref <- GetAssayData(AKIKB, assay = "RNA", slot = "data")
ref <- LogNormalize(ref)

## use singleR to compare AKI kidney biopsy data with our REANL SO ----
singler_results_kidney_ref <- SingleR::SingleR(
  test = test,
  ref = ref,
  labels = AKIKB$exp.ct, de.method="wilcox"
)

## annotate urine TEC in RENAL SO via new metadata $ref_labels ----
RENAL@meta.data$ref_labels <- singler_results_kidney_ref@listData$labels

## visualize in heatmap ----
RENAL.tbl <- RENAL@meta.data
RENAL.tbl$cell <- rownames(RENAL.tbl)
RENAL.tbl <- RENAL.tbl %>% arrange(renalcelltype)
RENAL.tbl$renalcelltype

names(renalcelltype_colors) <- levels(RENAL.tbl$renalcelltype)
names(AKIKBcolors) <- levels(AKIKB$exp.ct)
annot.colors <- list(biopsy = AKIKBcolors, urine =  renalcelltype_colors)
annot.col <- data.frame(
  biopsy = RENAL.tbl$AKIKB_labels,
  urine = RENAL.tbl$renalcelltype)
rownames(annot.col) <- rownames(RENAL.tbl)

p <- SingleR::plotScoreHeatmap(
  singler_results_kidney_ref,
  show.labels = F,
  annotation_col = annot.col, 
  annotation_colors = annot.colors, 
  fontsize = 24,
  annotation_legend = F,
  color = turbo(100),
  main = "Biopsy referencing of urine cells")

## save ----
png("supplfig12a.png", width=1500, height= 600, type="cairo-png") 
p
dev.off()

# Suppl Fig. 12B, reference annotation scores in UMAP ----

## create dataframe with UMAP coordinates and singleR reference annotation score for each cell in RENAL object ----
scores <- as.data.frame(singler_results_kidney_ref@listData$scores)
umaps <- as.data.frame(RENAL@reductions[["umap"]]@cell.embeddings)
score.tbl <- bind_cols(umaps, scores)
colnames(score.tbl) <- str_replace_all(colnames(score.tbl), "-", "_")

## function for Suppl Fig.12B ----
# score_plot() plots the singleR annotation score of a chosen reference cell type for each cell
score_plot <- function(cluster) { score.tbl <- arrange(score.tbl, eval(parse_expr(cluster)))
ggplot(score.tbl, aes_string(x="UMAP_1", y="UMAP_2", color= cluster)) +
  geom_point(size=0.5, alpha=1) +
  scale_color_viridis(limits=c(0,0.5), option = "H", name = "score") +
  theme_void() + 
  ggtitle(cluster) + 
  theme(plot.title = element_text(face= "bold", size = 20, hjust=0.5),
        legend.text = element_text(face= "bold", size = 20, hjust=0.5),
        legend.title = element_text(face= "bold", size = 20, hjust=0.5)) 
}

## plot scores for several TEC reference cell types ----
d2 <- score_plot("PT_New") +
  score_plot("TL_New")+
  score_plot("TAL_New")+
  score_plot("DCT_New")+
  score_plot("CNT_New")+
  score_plot("CD_PC_New")+
  score_plot("CD_IC_New")+  
  score_plot("PT")+
  score_plot("DTL")+
  score_plot("TAL")+
  score_plot("DCT")+
  score_plot("CNT")+
  score_plot("CD_PC")+
  score_plot("CD_IC")

for (i in 2:14) {d2[[i]] <- d2[[i]] + NoLegend()}

## compile and save ----
layout = "
ABCDEFG
HIJKLMN"

png("fig3ScorePlots.png", width=1500, height= 500, type="cairo-png") 
d2 + plot_layout(design = layout, guides = "collect")
dev.off()

# Fig.3D, comparison of location and cell state between tissue and urine samples ----

## assign nephron location ----
## assign broad nephron location to urine tubular epithelial cells based on singleR reference annotation (see Suppl. Fig. 12A above) 
RENAL$location<- NA
RENAL@meta.data[RENAL@meta.data$AKIKB_labels %in% c("PT", "PT_New"),]$location <- "PROX"
RENAL@meta.data[RENAL@meta.data$AKIKB_labels %in% c("PDC"),]$location <- "PDC"
RENAL@meta.data[RENAL@meta.data$AKIKB_labels %in% c("TAL", "TAL_New", "TL", "TL_New", "DTL"),]$location <- "LoH"
RENAL@meta.data[RENAL@meta.data$AKIKB_labels %in% c("CD-PC", "CD-PC_New", "CD-IC_New", "CD-IC"),]$location <- "CD"
RENAL@meta.data[RENAL@meta.data$AKIKB_labels %in% c("DCT", "DCT_New", "CNT_New", "CNT"),]$location <- "DIST"

RENAL$location <- factor(RENAL$location)
RENAL$location <- factor(RENAL$location, levels(RENAL$location)[c(4,5,3,2,1)])

## do the same for AKI kidney biopsy cells based on $exp.ct cell type annotation.
AKIKB$location <- NA
AKIKB@meta.data[AKIKB@meta.data$exp.ct %in% c("PT", "PT_New"),]$location <- "PROX"
AKIKB@meta.data[AKIKB@meta.data$exp.ct %in% c("PDC"),]$location <- "PDC"
AKIKB@meta.data[AKIKB@meta.data$exp.ct %in% c("TAL", "TAL_New", "TL", "TL_New", "DTL"),]$location <- "LoH"
AKIKB@meta.data[AKIKB@meta.data$exp.ct %in% c("CD-PC", "CD-PC_New", "CD-IC_New", "CD-IC"),]$location <- "CD"
AKIKB@meta.data[AKIKB@meta.data$exp.ct %in% c("DCT", "DCT_New", "CNT_New", "CNT"),]$location <- "DIST"

AKIKB$location <- factor(AKIKB$location)
AKIKB$location <- factor(AKIKB$location, levels(AKIKB$location)[c(4,5,3,2,1)])

## assign cell state ----

## we assume that cells annotated as "_New" cell types by the biopsy reference are injured/ injury-reactive cell states while other cells are more healthy
RENAL$state <- NA
RENAL@meta.data[RENAL@meta.data$AKIKB_labels %in% c("PT", "TAL", "TL",  "DTL", "CD-PC",  "CD-IC", "DCT", "CNT"),]$state <- "healthy"
RENAL@meta.data[RENAL@meta.data$AKIKB_labels %in% c("PT_New","TL_New", "TAL_New", "DCT_New", "CNT_New", "CD-PC_New", "CD-IC_New"),]$state <- "injury"

## do the same for AKI kidney biopsy cells based on $exp.ct cell type annotation.
AKIKB$state <- NA
AKIKB@meta.data[AKIKB@meta.data$exp.ct %in% c("PT", "TAL", "TL",  "DTL", "CD-PC",  "CD-IC", "DCT", "CNT"),]$state <- "healthy"
AKIKB@meta.data[AKIKB@meta.data$exp.ct %in% c("PT_New","TL_New", "TAL_New", "DCT_New", "CNT_New", "CD-PC_New", "CD-IC_New"),]$state <- "injury"

## DimPlots nephron location and cell state ----

a <- DimPlot(AKIKB, reduction = "umap", group.by = "location", label = F, repel = T, cols = locationcolors, shuffle = T) + theme_void() +theme(title = element_text(size = 20, face = "bold", vjust =0.5)) + NoLegend()
c <- DimPlot(RENAL, reduction = "umap", group.by = "location", label = F, repel = T, pt.size=0.25, cols = locationcolors) + theme_void() +theme(title = element_blank())+ NoLegend()
d <- DimPlot(AKIKB, reduction = "umap", group.by = "state", label = F, repel = T, cols = statecolors)+ theme_void()+ theme(title = element_text(size=20, face = "bold", vjust =0.5)) + NoLegend()
f <- DimPlot(RENAL, reduction = "umap", group.by = "state",  label = F, repel = T, pt.size=0.25, cols = statecolors) + theme_void()+ theme(title = element_blank()) + NoLegend()

## barplots of relative proportions for nephron location and cell state ----
location.df <- as.data.frame(tibble(
  Location = levels(factor(AKIKB$location)),
  Biopsy = table(AKIKB$location),
  Urine = table(RENAL$location)))
rownames(location.df) <- location.df$Location
location.df <- as.data.frame(t(location.df[,2:3]))
location.df$material <- rownames(location.df)
location.df <- location.df[,c(1,2,3,5,4,6)]
location.df <- pivot_longer(location.df, colnames(location.df)[1:5], names_to = "location")
b <- ggplot(location.df, aes(material, value)) +
  geom_bar(aes(fill = location), stat = 'identity', position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values = locationcolors) +
  scale_y_continuous(name = '% of single cells') +
  scale_x_discrete(limits = c("Urine", "Biopsy")) +
  coord_cartesian(clip = 'off') +
  theme_void() +
  theme(
    legend.position = 'left',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank() ) + 
  coord_flip() +
  
  guides(fill = guide_legend(reverse = T))

state.df <- as.data.frame(tibble(
  state = levels(factor(AKIKB$state)),
  Biopsy = table(AKIKB$state),
  Urine = table(RENAL$state)))
rownames(state.df) <- state.df$state
state.df <- as.data.frame(t(state.df[,2:3]))
state.df$material <- rownames(state.df)
state.df <- pivot_longer(state.df, colnames(state.df)[1:2], names_to = "state")
e <- ggplot(state.df, aes(material, value)) +
  geom_bar(aes(fill = state), stat = 'identity', position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values = statecolors) +
  scale_y_continuous(name = '% of single cells') +
  scale_x_discrete(limits = c("Urine", "Biopsy")) +
  coord_cartesian(clip = 'off') +
  theme_void() +
  theme(
    legend.position = 'left',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()) + 
  coord_flip()

## compile and save ----

layout2 <- "
AAAAAAAAAADDDDDDDDDD
AAAAAAAAAADDDDDDDDDD
AAAAAAAAAADDDDDDDDDD
AAAAAAAAAADDDDDDDDDD
AAAAAAAAAADDDDDDDDDD
AAAAAAAAAADDDDDDDDDD
AAAAAAAAAADDDDDDDDDD
AAAAAAAAAADDDDDDDDDD
AAAAAAAAAADDDDDDDDDD
AAAAAAAAAADDDDDDDDDD
BBBBBBBBBBEEEEEEEEEE
BBBBBBBBBBEEEEEEEEEE
CCCCCCCCCCFFFFFFFFFF
CCCCCCCCCCFFFFFFFFFF
CCCCCCCCCCFFFFFFFFFF
CCCCCCCCCCFFFFFFFFFF
CCCCCCCCCCFFFFFFFFFF
CCCCCCCCCCFFFFFFFFFF
CCCCCCCCCCFFFFFFFFFF
CCCCCCCCCCFFFFFFFFFF
CCCCCCCCCCFFFFFFFFFF
CCCCCCCCCCFFFFFFFFFF"

ggsave(plot = a+b+c+d+e+f+plot_layout(design =layout2, guides = "collect"), "fig3d.png", h = 2000, w = 2200, units = "px", scale = 1, type = "cairo-png")

# Suppl Fig. 11 ----

## Suppl. Fig 11 A-C ----
### Find cells with co-expression of SOX4 and SOX9 in AKI kidney biopsy (AKIKB) dataset ----
AKIKB_sox <- WhichCells(AKIKB, expression = SOX4 >= 1 & SOX9 >= 1)

AKIKB$sox <- "neg"
AKIKB@meta.data[AKIKB@assays[["RNA"]]@data@Dimnames[[2]] %in% AKIKB_sox, ]$sox <- "pos"

### visualize
a <- DimPlot(AKIKB, reduction = "umap", group.by = "exp.ct_nr", label = T, repel = T, label.size = 8) + 
  theme_void() + 
  scale_color_manual(labels = paste(levels(AKIKB_subset$exp.ct_nr), levels(AKIKB_subset$exp.ct2)), values = AKIKBcolors) + 
  theme() +
  ggtitle("AKI post mortem biopsy") +
  theme(legend.position = "top",
        legend.text = element_text(size = 20),
        plot.title = element_text(face= "bold", size = 30, hjust=0.5))
b <- DimPlot(AKIKB, cells.highlight = AKIKB_sox,  order= T, sizes.highlight = 0.25, raster = F) +
  theme_void() + 
  scale_color_manual(values = highlightcolors, labels=c("neg.", "pos."))+
  ggtitle("SOX4+ SOX9+ cells") +
  theme(legend.position = "top",
        legend.text = element_text(size = 20),
        plot.title = element_text(face= "bold", size = 30, hjust=0.5))
c <- DimPlot(AKIKB, cells.highlight = AKIKB_sox, cols = highlightcolors,  order= T, sizes.highlight = 0.25, split.by = "group", raster = F)+
  theme_void() +
  scale_color_manual(values = highlightcolors, labels=c("neg.", "pos."))+
  theme(text = element_text(face="bold", size = 30, hjust=0.5))+ NoLegend()

layout <- "AABB
AABB
CCCC
"
png("AKIKBProgDim.png", width = 1600, height = 1300, type="cairo-png")
c1+b+c+plot_layout(design=layout)
dev.off()

## Suppl. Fig 11 D ----

### find fractions of SOX4+SOX9+ cells in TEC subsets, compare between control, non-COVID AKI, and COVID-AKI ----
comparisons <- list(c("Control", "COVID"), c("COVID", "Non-COVID"), c("Control", "Non-COVID"))
d <- as.data.frame(table(tibble(AKIKB$orig.ident, AKIKB$sox, AKIKB$group, AKIKB$ct))) %>%
  pivot_wider(names_from = AKIKB.stem, values_from = Freq) %>%
  mutate(frac=pos/neg) %>%
  filter(neg !=0) %>%
  filter(AKIKB.finer.ct == "PT") %>%
  ggplot(aes(x=AKIKB.group, y=frac))+
  geom_boxplot(fill=paletteer_d("palettesForR::Pastels")[2])+
  geom_jitter()+
  theme_minimal()+
  ggtitle("PT")+
  scale_y_continuous(name = "SOX4+SOX9+\nfraction")+
  theme(legend.position = "bottom", 
        plot.title = element_text(size=30, vjust=0.5, hjust=0.5, face="bold"),
        axis.title.y = element_text(size=30,hjust=0.5), 
        axis.title.x =  element_text(size=30,hjust=10, vjust=-1), 
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18), 
        axis.text = element_text(vjust = 0.5, hjust=0.5, face= "bold", size=24))+
  stat_compare_means(comparisons = comparisons,  label = "p.signif", size= 10, stroke=2, hide.ns=F)
e <- as.data.frame(table(tibble(AKIKB$orig.ident, AKIKB$stem, AKIKB$group, AKIKB$finer.ct))) %>%
  pivot_wider(names_from = AKIKB.stem, values_from = Freq) %>%
  mutate(frac=pos/neg) %>%
  filter(neg !=0) %>%
  filter(AKIKB.finer.ct == "TAL") %>%
  ggplot(aes(x=AKIKB.group, y=frac))+
  geom_boxplot(fill=paletteer_d("palettesForR::Pastels")[9])+
  geom_jitter()+
  theme_minimal()+
  ggtitle("TAL")+
  scale_y_continuous(name = "SOX4+SOX9+\nfraction")+
  theme(legend.position = "bottom", 
        plot.title = element_text(size=30, vjust=0.5, hjust=0.5, face="bold"),
        axis.title.y = element_blank(), 
        axis.title.x =  element_text(size=30,hjust=10, vjust=-1), 
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18), 
        axis.text = element_text(vjust = 0.5, hjust=0.5, face= "bold", size=24))+
  stat_compare_means(comparisons = comparisons,  label = "p.signif", size= 10, stroke=2, hide.ns=F)
f <- as.data.frame(table(tibble(AKIKB$orig.ident, AKIKB$stem, AKIKB$group, AKIKB$finer.ct))) %>%
  pivot_wider(names_from = AKIKB.stem, values_from = Freq) %>%
  mutate(frac=pos/neg) %>%
  filter(neg !=0) %>%
  filter(AKIKB.finer.ct == "CD-PC") %>%
  ggplot(aes(x=AKIKB.group, y=frac))+
  geom_boxplot(fill=paletteer_d("palettesForR::Pastels")[12])+
  geom_jitter()+
  theme_minimal()+
  ggtitle("CD-PC")+
  scale_y_continuous(name = "SOX4+SOX9+\nfraction")+
  theme(legend.position = "bottom", 
        plot.title = element_text(size=30, vjust=0.5, hjust=0.5, face="bold"),
        axis.title.y = element_blank(), 
        axis.title.x =  element_text(size=30,hjust=10, vjust=-1), 
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18), 
        axis.text = element_text(vjust = 0.5, hjust=0.5, face= "bold", size=24))+
  stat_compare_means(comparisons = comparisons,  label = "p.signif", size= 10, stroke=2, hide.ns=F)

png("AKIKBProgBox.png", width = 1500, height = 600, type="cairo-png")
d+e+f
dev.off()