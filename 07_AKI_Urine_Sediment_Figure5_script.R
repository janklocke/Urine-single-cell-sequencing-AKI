# description ----

# This is the R script for generating Figure 5 and Supplemental Figures 11, 12 in 
# "Urinary single-cell sequencing captures intrarenal injury and repair processes 
# in human acute kidney injury" by Klocke et al. 

# loading packages and data ----

## loading packages ----
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(celldex)
library(patchwork)
library(RColorBrewer)
library(viridis)
library(paletteer) 
library(rrvgo)
library(topGO)
library(org.Hs.eg.db)
library(progeny)
library(ComplexHeatmap)
library(scales)

## loading data ----
RENAL <-readRDS("~/SO_kidney_urine_cells.rds")

## loading color schemes ----
highlightcolors <- c(paletteer_d("palettesForR::Pastels")[9], paletteer_d("rcartocolor::SunsetDark")[4]) 
celltype_colors <- colorRampPalette(paletteer_d("palettesForR::Pastels"))(19)
cellgroup_colors <- brewer.pal(5, "Paired")

renalcelltype_colors <- 
  c(paletteer_d("ggsci::purple_material")[c(3,5)],                  #PDC colors
    paletteer_d("ggsci::teal_material")[c(8,6,4,3)],                   #healthy TEC PT and LOH colors   
    paletteer_d("ggsci::light_blue_material")[c(3,5,7,9)],          #healthy TEC DT and CD colors
    paletteer_d("ggsci::amber_material")[c(3,4,6,8,10)],            #EMT/infl TEC
    paletteer_d("ggsci::red_material")[c(2,5,8)],                   #OxyStress TEC
    paletteer_d("ggsci::pink_material")[4],                         #prolif TEC
    paletteer_d("ggsci::light_green_material")[c(9,7,5,3)])         #progenitors 
names(renalcelltype_colors) <- str_replace(levels(factor(RENAL$renalcelltype)), "DCT", "TAL/DCT")
names(renalcelltype_colors) <- str_replace(names(renalcelltype_colors), "PCT", "PT")
names(renalcelltype_colors) <- str_replace(names(renalcelltype_colors), "DTL", "TEC_VCAM1")
renalcelltype_colors2 <- 
  c(paletteer_d("ggsci::purple_material")[c(3,5)],                  #PDC colors
    paletteer_d("ggsci::teal_material")[c(8,6,4,3)],                   #healthy TEC PT and LOH colors   
    paletteer_d("ggsci::light_blue_material")[c(3,5,7,9)],          #healthy TEC DT and CD colors
    paletteer_d("ggsci::amber_material")[c(3,4,6,8,10)],            #EMT/infl TEC
    paletteer_d("ggsci::red_material")[c(2,5,8)],                   #OxyStress TEC
    paletteer_d("ggsci::pink_material")[4],                         #prolif TEC
    paletteer_d("ggsci::light_green_material")[c(9,7,5,3)])         #progenitors 
names(renalcelltype_colors2) <- levels(factor(RENAL$renalcelltype))

## functions ----
number_ticks <- function(n) {function(limits) pretty(limits, n)}

# preparing the SO ----
RENAL$AKI_timepoint.fine <- as.numeric(RENAL$AKI_timepoint.fine)

## subsetting to TEC of interest ----
# excluding other cell types (podocytes) for Fig. 5A-C and subsetting to prox. tubule cells for Suppl. Fig. 11
table(RENAL$renalcelltype)
RENAL <- SetIdent(RENAL, value = "renalcelltype")
RENALsub <- subset(RENAL, idents= c("PDC", "PDC_inj"), invert=T)
RENALsub_PT <- subset(RENAL, idents= c("PT"), invert=F)

# Fig. 5A-C ----

## transfer data to monocle object ----
MO <- as.cell_data_set(RENALsub)
MO <- cluster_cells(MO)
MO <- learn_graph(MO)
## assign nodes in high UMOD/SLC12A1 (TAL marker) and low LCN2 (injury marker) expression as root.
MO <- order_cells(MO, reduction_method = "UMAP")

## Fig 5AB, visualize in Dimplot ----
MO@colData$renalcelltype <- str_replace(str_replace(MO@colData$renalcelltype, "PCT", "PT"), "DTL", "TEC_VCAM1")
a <- plot_cells(MO, alpha = 0.5, label_cell_groups=FALSE, color_cells_by = "renalcelltype",
                show_trajectory_graph=T, cell_size = 2, cell_stroke=0, graph_label_size =4, trajectory_graph_segment_size=1.2, label_branch_points =F, label_roots=T, label_leaves=F) +
  theme_void()  +
  scale_color_manual(values=renalcelltype_colors2, labels = str_replace(str_replace(names(renalcelltype_colors), "PCT", "PT"), "DTL", "TEC_VCAM1")[3:23], limits = names(renalcelltype_colors2)[3:23])

b <- plot_cells(MO, color_cells_by = "pseudotime", alpha = 0.3, label_cell_groups=FALSE,
                show_trajectory_graph=T, cell_size = 2, cell_stroke=0, graph_label_size =4,  trajectory_graph_segment_size=1.2, label_branch_points =F, label_roots=T, label_leaves=F) +
  theme_void()

### save plot ----
ggsave(plot = a+b&
         theme(text=element_text(size=16)), "fig5ab.png", w = 4000, h = 1500, units = "px", scale = 1, type = "cairo-png")

# Fig. 5C ----

## get marker gene expression from SO ----
RENALsub@active.assay <- "RNA"
assayRNA <- GetAssayData(RENALsub)
assayRNA <- NormalizeData(assayRNA)

## add marker gene expression to MO ----
MO@colData$UMOD <- assayRNA["UMOD",]
MO@colData$SLC12A1 <- assayRNA["SLC12A1",]
MO@colData$LCN2 <- assayRNA["LCN2",]
MO@colData$PAX2 <- assayRNA["PAX2",]
MO@colData$VCAM1 <- assayRNA["VCAM1",]

## Visualize expression along trajectory ----
d <- plot_cells(MO, alpha = 0.2, label_cell_groups=FALSE, color_cells_by ="UMOD",
                show_trajectory_graph=T, cell_size = 1.5, cell_stroke=0, graph_label_size =4, 
                trajectory_graph_segment_size=0.5, label_branch_points =F, label_roots=F, label_leaves=F)+
  ggtitle("UMOD")
e <- plot_cells(MO, alpha = 0.2, label_cell_groups=FALSE, color_cells_by ="LCN2",
                show_trajectory_graph=T, cell_size = 1.5, cell_stroke=0, graph_label_size =4, 
                trajectory_graph_segment_size=0.5, label_branch_points =F, label_roots=F, label_leaves=F)+
  ggtitle("LCN2")
f <- plot_cells(MO, alpha = 0.2, label_cell_groups=FALSE, color_cells_by ="VCAM1",
                show_trajectory_graph=T, cell_size = 1.5, cell_stroke=0, graph_label_size =4, 
                trajectory_graph_segment_size=0.5, label_branch_points =F, label_roots=F, label_leaves=F)+
  ggtitle("VCAM1")
g <- plot_cells(MO, alpha = 0.2, label_cell_groups=FALSE, color_cells_by ="PAX2",
                show_trajectory_graph=T, cell_size = 1.5, cell_stroke=0, graph_label_size =4, 
                trajectory_graph_segment_size=0.5, label_branch_points =F, label_roots=F, label_leaves=F)+
  ggtitle("PAX2")

### compile and save ----
pw <- d+e+f+g+plot_layout(widths=c(1,1,1,1), heights= c(1), guides = "collect") & 
  theme_void() & 
  scale_color_gradient(name = "expression", low = "#ebebeb", high = "#130a45", limits = c(0,1.5), oob= squish) &
  theme(text=element_text(size=20))

ggsave(plot = pw, "fig5c.png", w = 4000, h = 900, units = "px", scale = 1, type = "cairo-png")

# Figure 5D, pathway activity by cluster via PROGENy ----
## Here, we apply PROGENy to infer the activity of 14 relevant signaling pathways based on the suggested workflow
## https://bioconductor.org/packages/release/bioc/vignettes/progeny/inst/doc/ProgenySingleCell.html

## subset TEC cluster from RENAL dataset ----
RENAL  <- SetIdent(RENAL, value="renalcelltype")
RENAL_TEC <- subset(RENAL, idents = c("TEC_prg1", "TEC_MMP7", "TEC_prg2", "TEC_TXNRD1", 
                                      "TEC_GCLM", "TEC_IL32", "TEC_prg4", "TAL/DCT", 
                                      "PT", "TEC_prlf", "ATL", "DTL", "TEC_LIF", "TEC_dmg", 
                                      "CD-PC",  "CNT/CD","TEC_prg3", "TAL", "CD-IC", "TEC_MT"), invert=F)

## assign cells and identities ----
RENAL_TEC <- SetIdent(RENAL_TEC, value="renalcelltype")
CellsClusters <- data.frame(Cell = names(Idents(RENAL_TEC)), 
                            CellType = as.character(Idents(RENAL_TEC)),
                            stringsAsFactors = FALSE)
## compute the Progeny activity scores and add them to the Seurat object as a new assay called Progeny ----
RENAL_TEC <- progeny(RENAL_TEC, scale=FALSE, organism="Human", top=500, perm=1, 
                     return_assay = TRUE)

## scale the pathway activity scores ----
RENAL_TEC <- Seurat::ScaleData(RENAL_TEC, assay = "progeny") 

## transform Progeny scores into a data frame ----
progeny_scores_df <- 
  as.data.frame(t(GetAssayData(RENAL_TEC, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell")%>%
  gather(Pathway, Activity, -Cell) 

## match Progeny scores with the cell clusters ----
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

## summarize the Progeny scores by cellpopulation ----
summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

## prepare the data for the plot ----
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 
rownames(summarized_progeny_scores_df) <- str_replace(rownames(summarized_progeny_scores_df), "DCT", "TAL/DCT")
paletteLength <- 100
myColor <- viridis(101)

progenyBreaks <- c(seq(min(summarized_progeny_scores_df), 0, 
                       length.out=ceiling(paletteLength/2) + 1),
                   seq(max(summarized_progeny_scores_df)/paletteLength, 
                       max(summarized_progeny_scores_df), 
                       length.out=floor(paletteLength/2)))

## create heatmap and save ----
progeny_hmap <- pheatmap(t(data.matrix(summarized_progeny_scores_df)),fontsize=60, 
                         fontsize_row = 60, legend=T,
                         color=myColor, breaks = progenyBreaks, 
                         cellwidth = 80, cellheight = 80,
                         main = "pathway activity (PROGENy)", 
                         treeheight_col = 0,  border_color = NA)

png("fig5d.png", width=2000, height=2000, type="cairo-png")
progeny_hmap
dev.off()

# Figure 5E, treemap plot of enriched Gene ontology (GO) terms in PAX2 marker positive cells in RENAL dataset ----

## create a new metadata column with info on PAX2 marker expression from Figure 5B ----
RENAL$pax2 <- "neg"
RENAL@meta.data[RENAL@assays[["SCT"]]@data@Dimnames[[2]] %in% RENAL_pax2, ]$pax2 <- "pos"

## determine DEG of pax2 marker positive cells with FindMarkers() ----
markers.pax2 <- FindMarkers(RENAL, ident.1 =  "pos", ident.2 = "neg", group.by = "pax2", test.use = "wilcox", logfc.threshold = 0.5, min.pct = 0.01, assay = "SCT", slot="data")

## determine enriched gene ontology (GO) terms ----
## here, we use the topGP package to determine enriched GO terms 
## https://bioconductor.org/packages/release/bioc/html/topGO.html

### topGO function ----
## the function below determines the 30 top "Biological Process" (BP) ontology and tests significance with both with the "elim" and 
## "classic" methods. 
## "elim" =  processes the GO terms by traversing the GO hierarchy from bottom to top, (more conservative, less false positives)
## "classic" = each GO term is tested independently, not taking the GO hierarchy into account.
topGO_BP <- function(genes) {
  geneUniverse <- rownames(RENAL@assays$RNA)
  genesOfInterest <- rownames(genes[genes$avg_log2FC > 0.5, ])
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse
  onts <- "BP"
  tab <- as.list(onts)
  names(tab) <- onts
  GOdata <- new("topGOdata",
                description = "GOanalysis",
                ontology = "BP",
                allGenes = geneList,
                annot = annFUN.org,
                mapping = "org.Hs.eg.db",
                ID = "SYMBOL",
                nodeSize = 20)  
  res.result1 <- runTest(GOdata, statistic = "fisher", algorithm = "elim")
  res.result2 <- runTest(GOdata, statistic = "fisher", algorithm = "classic")
  tab <- data.frame(GenTable(GOdata, Fisher.elim = res.result1,
                             Fisher.classic = res.result2,
                             orderBy = "Fisher.elim" , topNodes = 30))
  topGOResults <- plyr::rbind.fill(tab)
  topGOResults.df <- as.data.frame(topGOResults)
  topGOResults.df$gene_ratio <- topGOResults.df$Significant / topGOResults.df$Annotated
  # modification appropriate for plot
  topGOResults.df$Fisher.elim <- as.numeric(topGOResults.df$Fisher.elim)
  topGOResults.df$Fisher.classic <- as.numeric(topGOResults.df$Fisher.classic)
  topGOResults.df$Term <- factor(topGOResults.df$Term, levels = rev(unique(topGOResults.df$Term)))  
  return(topGOResults.df)
}

### apply function on DEG ----
topGO_pax2 <- topGO_BP(markers.pax2)

## reduce terms with rrvgo ----
## the rrvgo package groups GO terms by similarity, therefore reducing redundant terms
## https://ssayols.github.io/rrvgo/

### calculate a similarity matrix for terms found above ----
l <- topGO_pax2
simMatrix <- calculateSimMatrix(l$GO.ID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")

### set -log10(p value) calculated by elim method of each GO term as score to determine most significant term among groups)
scores <- setNames(-log10(l$Fisher.elim), l$GO.ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.8,
                                orgdb="org.Hs.eg.db")

### create treemap plot, save ----
png("fig5e.png", width = 2000, height = 1600, type="cairo-png")
treemapPlot(reducedTerms, border.lwds=c(4,2), fontsize.label=c(30,10), fontface.labels="bold", inflate.labels =T)
dev.off()

# Supplements ######

# Suppl. Fig. 11 ----

## construct new monocle object of RENALsub_PT ----
RENALsub_PT
data <- as(as.matrix(GetAssayData(RENALsub_PT, assay = "RNA", slot = "data")), 'sparseMatrix')
pd <- data.frame(RENALsub_PT@meta.data)
#keep only the columns that are relevant
pData <- pd %>% dplyr::select(orig.ident, nCount_RNA, nFeature_RNA, percent.mt, renalcelltype)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

#Construct monocle cds
MO <- new_cell_data_set(expression_data = data, cell_metadata = pData, gene_metadata = fData)

#preprocess
MO = preprocess_cds(MO, num_dim = 10, norm_method = "size_only", pseudo_count = 0)
MO <- align_cds(MO)
MO = reduce_dimension(MO, reduction_method = "UMAP")
MO = cluster_cells(MO)
MO = learn_graph(MO)

# assign nodes with high LRP2/GPX3 expression (PT markers) as roots
MO = order_cells(MO, reduction_method = "UMAP")

## Suppl. Fig. 11B-C, visualize trajectory and pseudotime ----
f <- plot_cells(MO, alpha = 0.5, label_cell_groups=FALSE, genes = "GPX3",
                show_trajectory_graph=T, cell_size = 4, cell_stroke=0, graph_label_size =4, trajectory_graph_segment_size=1.2, label_branch_points =F, label_roots=T, label_leaves=T) +
  theme_void()  + scale_color_viridis(name = "expression", option="D", limits = c(0,3), oob= squish)
g <- plot_cells(MO, color_cells_by = "pseudotime", alpha = 0.3, label_cell_groups=FALSE,
                show_trajectory_graph=T, cell_size = 4, cell_stroke=0, graph_label_size =4,  trajectory_graph_segment_size=1.2, label_branch_points =F, label_roots=T, label_leaves=T) +
  theme_void()

### compile and save ----
ggsave(plot = f+g&theme(text=element_text(size=20)), "supplfig11bc.png", w = 2500, h = 1500, units = "px", scale = 1, type = "cairo-png")


## Suppl. Fig. 11D, visualize marker gene expression ----
a <- plot_cells(MO, genes= c("GPX3"),
                alpha = 0.3, label_cell_groups=FALSE,
                show_trajectory_graph=F, cell_size = 2.5, trajectory_graph_segment_size=1.5, 
                label_branch_points =F, label_roots=F, label_leaves=F, graph_label_size=4)
b <- plot_cells(MO, genes= c("LRP2"),
                alpha = 0.3, label_cell_groups=FALSE,
                show_trajectory_graph=F, cell_size = 2.5, trajectory_graph_segment_size=1.5, 
                label_branch_points =F, label_roots=F, label_leaves=F, graph_label_size=4)
c <- plot_cells(MO, genes= c("HAVCR1"),
                alpha = 0.3, label_cell_groups=FALSE,
                show_trajectory_graph=F, cell_size = 2.5, trajectory_graph_segment_size=1.5, 
                label_branch_points =F, label_roots=F, label_leaves=F, graph_label_size=4)
d <- plot_cells(MO, genes= c("VCAM1"),
                alpha = 0.3, label_cell_groups=FALSE,
                show_trajectory_graph=F, cell_size = 2.5, trajectory_graph_segment_size=1.5, 
                label_branch_points =F, label_roots=F, label_leaves=F, graph_label_size=4)
e <- plot_cells(MO, genes= c("PROM1"),
                alpha = 0.3, label_cell_groups=FALSE,
                show_trajectory_graph=F, cell_size = 2.5, trajectory_graph_segment_size=1.5, 
                label_branch_points =F, label_roots=F, label_leaves=F, graph_label_size=4)

### compile and save ----
pw <- b+c+d+e+plot_layout(widths=c(1,1,1,1), guides = "collect") &
  scale_color_viridis(name = "expression", option="D", limits = c(0,3), oob= squish)&theme_void()&theme(text=element_text(size=20))

ggsave(plot = pw, "supplfig11d.png", w = 2500, h = 700, units = "px", scale = 1, type = "cairo-png")


# Suppl. Fig. 12A, stem marker positive cells in RENAL dataset ----

## Find cells expressing high levels of select marker genes ---- 
RENAL_pax2 <- WhichCells(RENAL, expression = PAX2 >= 1)
RENAL_sox9 <- WhichCells(RENAL, expression = SOX9 >= 1)
RENAL_cd24 <- WhichCells(RENAL, expression = CD24 >= 2)
RENAL_prom1 <- WhichCells(RENAL, expression = PROM1 >= 1)

## Highlight marker+ cells in UMAP ----

a <- DimPlot(RENAL, cells.highlight = RENAL_sox9, cols = highlightcolors,  pt.size = 0.5, sizes.highlight = 0.5,order= T) + 
  theme_void() + 
  ggtitle("SOX9+") + 
  theme(title = element_text(size = 16, face = "bold", hjust = 0.5)) +
  NoLegend()
b <- DimPlot(RENAL, cells.highlight = RENAL_pax2, cols = highlightcolors, pt.size = 0.5, sizes.highlight = 0.5,order= T) + 
  theme_void() + 
  ggtitle("PAX2+") + 
  theme(title = element_text(size = 16, face = "bold", hjust = 0.5)) +
  NoLegend()
c <- DimPlot(RENAL, cells.highlight = RENAL_cd24, cols = highlightcolors,  pt.size = 0.5, sizes.highlight = 0.5,order= T) + 
  theme_void() + 
  ggtitle("CD24+") + 
  theme(title = element_text(size = 16, face = "bold", hjust = 0.5))+
  NoLegend()
d <- DimPlot(RENAL, cells.highlight = RENAL_prom1, cols = highlightcolors, pt.size = 0.5, sizes.highlight = 0.5, order= T) + 
  theme_void() + 
  ggtitle("PROM1+") + 
  theme(title = element_text(size = 16, face = "bold", hjust = 0.5)) +
  NoLegend()

ggsave(plot = a+b+c+d+plot_layout(heights=c(1,1), widths=c(1,1), guides="collect"), "fig5b.png", w = 2000, h = 2000, units = "px", scale = 1, type = "cairo-png")

# Suppl. Fig. 12C, violin plot of cell differentiation markers in RENAL dataset ----

p <- 
  VlnPlot(RENAL, "PLCG2", group.by = "renalcelltype", cols = renalcelltype_colors, pt.size = 0) + NoLegend() + 
  VlnPlot(RENAL, "SOX4", group.by = "renalcelltype", cols = renalcelltype_colors, pt.size = 0) + NoLegend() + 
  VlnPlot(RENAL, "FOS", group.by = "renalcelltype", cols = renalcelltype_colors, pt.size = 0) + NoLegend() +  
  VlnPlot(RENAL, "JUN", group.by = "renalcelltype", cols = renalcelltype_colors, pt.size = 0) + NoLegend() +
  plot_layout(heights = c(1,1,1,1))

for (i in 1:4) { p[[i]] <- p[[i]] + scale_y_continuous(breaks=number_ticks(2)) +
  scale_x_discrete(limits=levels(RENAL$renalcelltype)[c(1:13,18,16,19,14,15,23,22,21,20)], 
                   labels=str_replace(str_replace(str_replace(levels(RENAL$renalcelltype)[c(1:13,18,16,19,14,15,23,22,21,20)], "DCT", "TAL/DCT"), "PCT", "PT"), "DTL", "TEC_VCAM1")) +
  theme(title = element_text(size = 16), 
        axis.text.x = element_text(size=16, angle=90, vjust = 0.5),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=12), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  NoLegend() 
}
for (i in 1:3) { p[[i]] <- p[[i]] + scale_y_continuous(breaks=number_ticks(2)) +
  theme(title = element_text(size = 16), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=12), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  NoLegend() }
for (i in 2) { p[[i]] <- p[[i]] + scale_y_continuous(breaks=number_ticks(2)) +
  theme(title = element_text(size = 16), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=12), 
        axis.title.y = element_text(size=12, hjust=2, angle=90),
        axis.title.x = element_blank()) +
  NoLegend() }

ggsave(plot = p, "fig5a.png", w = 2000, h = 2000, units = "px", scale = 1, type = "cairo-png")

