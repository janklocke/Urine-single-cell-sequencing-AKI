# description ----

# This is the R script for generating Figure 2 and Supplemental Figure 6-8, 17 in 
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
library(EnhancedVolcano)

## loading data ----
URINE <- readRDS("~/SO_all_urine_cells.rds")
RENAL <-readRDS("~/SO_kidney_urine_cells.rds")

### renaming some clusters post-hoc 
x <- RENAL$renalcelltype
x <- str_replace(x, "PCT", "PT")
x <- str_replace(x, "TEC_emt", "TEC_inj")
x <- str_replace(x, "DTL", "TEC_VCAM1")
RENAL$renalcelltype2 <- x
RENAL$renalcelltype <- factor(RENAL$renalcelltype, levels(factor(RENAL$renalcelltype))[c(5,6,7,23,1,8,9,4,3,2,15,13,14,12,16,22,11,10,21,17,18,19,20)])

## loading color schemes ----
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
names(renalcelltype_colors) <- levels(RENAL$renalcelltype_nr)

renalcelltype_cnn_colors <- 
  c(paletteer_d("ggsci::purple_material")[c(3,5)],                  #PDC colors
    paletteer_d("ggsci::teal_material")[c(8,6,4,3)],                #healthy TEC PT and LOH colors   
    paletteer_d("ggsci::light_blue_material")[c(3,5,7,9)],          #healthy TEC DT and CD colors
    rep("lightgrey",13))      

renalcelltype_emt_colors <- 
  c(rep("lightgrey",10),
    paletteer_d("ggsci::amber_material")[c(3,4,6,8,10)],         
    rep("lightgrey",8))

renalcelltype_oxy_colors <- 
  c(rep("lightgrey",15),
    paletteer_d("ggsci::red_material")[c(2,5,8)],          
    rep("lightgrey",5))

renalcelltype_prlf_colors <- 
  c(rep("lightgrey",18),
    paletteer_d("ggsci::pink_material")[4],        
    rep("lightgrey",4))

renalcelltype_prg_colors <- 
  c(rep("lightgrey",19),
    paletteer_d("ggsci::light_green_material")[c(9,7,5,3)])

## functions ----
number_ticks <- function(n) {function(limits) pretty(limits, n)}

# Figure 2A, Dimensional plots (UMAP) ----

## Figure 2A Dimplot urinary renal epithelial cells ----
p1 <- DimPlot(RENAL, group.by = "renalcelltype_nr", cols = renalcelltype_colors, label = T, repel = T, label.size = 10, pt.size = 0.5) + 
  labs(title = "Urinary renal epithelial cells") +
  theme_void() + 
  scale_color_manual(labels = paste(levels(RENAL$renalcelltype_nr), levels(RENAL$renalcelltype)), values = renalcelltype_colors) + 
  theme(legend.position = "right",
        legend.text = element_text(size = 20),
        plot.title = element_text(face= "bold", size = 20, hjust=0.5))+
  NoLegend()

### save plot
ggsave(plot = p1, "fig2a.png", h = 2400, w = 2400, units = "px", scale = 1, type = "cairo-png")

## Figure 2A,  Dimplot of subsetted cells in URINE dataset ----
p2 <- DimPlot(URINE, group.by = "cellgroup", cols = c("lightgrey", "lightgrey", cellgroup_colors[3:4], "lightgrey")) +
  theme_void() + 
  theme(plot.title = element_blank()) + 
  NoLegend()

### save plot
ggsave(plot = p2, "fig2urine.png", h = 1000, w = 1000, units = "px", scale = 1, type = "cairo-png")


# Fig. 2B, nephorn segment specific clusters ----

## Fig. 2B, Dimplot ----
c1 <- DimPlot(RENAL, group.by = "renalcelltype_nr", cols = renalcelltype_cnn_colors) +
  theme_void() +
  theme(plot.title = element_blank()) + 
  NoLegend()

### save plot
c <- c1 + plot_spacer() + plot_annotation(title = "canonical TEC (TEC_cnn)", 
                                          theme = theme(plot.title = element_text(size = 18, face="bold")))

ggsave(plot = c, "fig2cnn1.png", w = 2400, h = 1200, units = "px", scale = 1, type = "cairo-png")


## Fig. 2B, Violinplots ----

### combining all non-segment specific clusters for violin plot
RENAL$cnngroup <- as.character(RENAL$renalcelltype)
RENAL@meta.data[RENAL@meta.data$renalcelltype %in% c("TEC_MMP7","TEC_TXNRD1","TEC_GCLM","TEC_IL32","TEC_LIF", "TEC_prg1", "TEC_prg2", "TEC_prg3", "TEC_prg4", "TEC_dmg", "TEC_prlf", "TEC_IFIT", "TEC_MT"),]$cnngroup <- "other"
RENAL$cnngroup <- factor(RENAL$cnngroup, levels(factor(RENAL$cnngroup))[c(5,6,7,8,11,1,9,10,4,3,2)])

### violin plots for segment specific markers
c2 <- 
  VlnPlot(RENAL, "PODXL", group.by = "cnngroup", cols = c("lightgrey", renalcelltype_cnn_colors), pt.size = 0) + NoLegend() + 
  VlnPlot(RENAL, "GPX3", group.by = "cnngroup", cols = c("lightgrey", renalcelltype_cnn_colors), pt.size = 0) + NoLegend() + 
  VlnPlot(RENAL, "VCAM1", group.by = "cnngroup", cols = c("lightgrey", renalcelltype_cnn_colors), pt.size = 0) + NoLegend() +  
  VlnPlot(RENAL, "CLDN10", group.by = "cnngroup", cols = c("lightgrey", renalcelltype_cnn_colors), pt.size = 0) + NoLegend() +
  VlnPlot(RENAL, "SLC12A1", group.by = "cnngroup", cols = c("lightgrey", renalcelltype_cnn_colors), pt.size = 0) + NoLegend() + 
  VlnPlot(RENAL, "UMOD", group.by = "cnngroup", cols = c("lightgrey", renalcelltype_cnn_colors), pt.size = 0) + NoLegend() + 
  VlnPlot(RENAL, "KCNJ1", group.by = "cnngroup", cols = c("lightgrey", renalcelltype_cnn_colors), pt.size = 0) + NoLegend() + 
  VlnPlot(RENAL, "AQP2", group.by = "cnngroup", cols = c("lightgrey", renalcelltype_cnn_colors), pt.size = 0) + NoLegend() +  
  VlnPlot(RENAL, "ATP6V0D2", group.by = "cnngroup", cols = c("lightgrey", renalcelltype_cnn_colors), pt.size = 0)+ NoLegend() + 
  plot_layout(heights = c(1,1,1,1,1,1,1,1,1))

for (i in 1:8) { c2[[i]] <- c2[[i]] + scale_y_continuous(breaks=number_ticks(2)) +
  theme(title = element_text(size = 16), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=16), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  NoLegend() }
c2[[9]] <- c2[[9]] + scale_y_continuous(breaks=number_ticks(2)) +
  scale_x_discrete(limits= c("other", "PDC", "PDC_inj", "PT", "TEC_VCAM1", "ATL", "TAL", "TAL/DCT", "CNT/CD", "CD-PC", "CD-IC")) +
  theme(title = element_text(size = 16), 
        axis.text.x = element_text(size=16, angle=90, hjust=1, vjust=0.5),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=16), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  NoLegend() 

### save plot
ggsave(plot = c2, "fig2cnn2.png", w = 800, h = 3000, units = "px", scale = 1, type = "cairo-png")


# Fig. 2C, INJ clusters ----

## Fig. 2C, Dimplot ----
d1 <- DimPlot(RENAL, group.by = "renalcelltype_nr", cols = renalcelltype_emt_colors) +
  theme_void() +
  theme(plot.title = element_blank()) + 
  NoLegend()

## Fig. 2C, Featureplots ----
d2 <- FeaturePlot(RENAL, c("MMP7", "VIM", "IGFBP7", "SPP1", "CCL2", "IL32", "IFITM3", "CD40", "LCN2", "HAVCR1", "TIMP1", "KRT8"), order=T, pt.size = 0.25, ncol =4) & 
  scale_colour_viridis(option = "turbo") & theme_void() & 
  theme(title= element_text(size=12)) & NoLegend()

## Fig. 2C, assembly + save ----
d <- (d1 + d2) + 
  plot_annotation(title = "Injury, Inflammation and tissue rearrangement (TEC_inj)", 
                  theme = theme(plot.title = element_text(size = 18, face="bold"))) + 
  plot_layout(widths = c(3,4))

ggsave(plot = d, "fig2c.png", w = 2600, h = 1200, units = "px", scale = 1, type = "cairo-png")


# Fig. 2D, oxidative stress clusters ----

## Fig. 2D, Dimplot ----
f1 <- DimPlot(RENAL, group.by = "renalcelltype_nr", cols = renalcelltype_oxy_colors) +
  theme_void() +
  theme(plot.title = element_blank()) + 
  NoLegend()

## Fig. 2D, Featureplots ----
f2 <- FeaturePlot(RENAL, c("TXNRD1", "GCLM", "SLC7A11"), order=T, pt.size = 1, ncol =1) & 
  scale_colour_viridis(option = "turbo") & 
  theme_void() & theme(title= element_text(size=12)) & NoLegend()

## Fig. 2D, assembly + save ----
f <- (f1 + f2) + 
  plot_annotation(title = "Oxidative stress / damage (TEC_str)", 
                  theme = theme(plot.title = element_text(size = 18, face="bold"))) + 
  plot_layout(widths = c(3,1))

ggsave(plot = f, "fig2d.png", w = 1400, h = 1200, units = "px", scale = 1, type = "cairo-png")


# Fig. 2E, progenitor clusters ----

## Fig. 2E, Dimplot ----
e1 <- DimPlot(RENAL, group.by = "renalcelltype_nr", cols = renalcelltype_prg_colors) +
  theme_void() +
  theme(plot.title = element_blank()) + 
  NoLegend()

## Fig. 2E, Featureplots ----
e2 <- FeaturePlot(RENAL, c("SOX4", "SOX9", "PAX2", "PLCG2", "UNCX", "HES1", "KLF4", "EMX2", "ELF3", "FGF9", "CXCL3", "CXCR4"), order=T, pt.size = 1, ncol =4) & 
  scale_colour_viridis(option = "turbo") & 
  theme_void()& theme(title= element_text(size=12)) & NoLegend()

## Fig. 2E, assembly + save ----
e <- (e1 + e2) + 
  plot_annotation(title = "Progenitor-like cells (TEC_prg)", 
                  theme = theme(plot.title = element_text(size = 18, face="bold"))) + 
  plot_layout(widths = c(3,4))

ggsave(plot = e, "fig2e.png", w = 2600, h = 1200, units = "px", scale = 1, type = "cairo-png")


# Fig. 2F, proliferating cluster ----

## Fig. 2F, Dimplot ----
g1 <- DimPlot(RENAL, group.by = "renalcelltype_nr", cols = renalcelltype_prlf_colors) +
  theme_void() +
  theme(plot.title = element_blank()) + 
  NoLegend()

## Fig. 2F, Featureplots ----
g2 <- FeaturePlot(RENAL, c("MKI67", "CENPF", "CCNB1"), order=T, pt.size = 1, ncol =1) & 
  scale_colour_viridis(option = "turbo") & theme_void() & theme(title= element_text(size=12)) & NoLegend()

## Fig. 2F, assembly + save ----
g <- (g1 + g2) + 
  plot_annotation(title = "Proliferation (TEC_prlf)", 
                  theme = theme(plot.title = element_text(size = 18, face="bold"))) + 
  plot_layout(widths = c(3,1))
ggsave(plot = g, "fig2f.png", w = 1400, h = 1200, units = "px", scale = 1, type = "cairo-png")


# Suppl. Fig. 6 ----

## Suppl. Fig. 6BC, Dimplots of genes/cell and percentage of mitochondrial RNA ----
p3 <- FeaturePlot(RENAL, c("nFeature_SCT"), order = T) + theme_void() + ggtitle("genes/cell") + theme(legend.position = "left") 
p4 <- FeaturePlot(RENAL, c("percent.mt"), order = T, cols = c("lightgrey", "red")) + theme_void() + ggtitle("% mitochondrial\nfeatures") + theme(legend.position = "left") 

## Suppl. Fig. 6D, Dotplots of marker gene expression ----

### assigning marker genes ----
RENALMarkerGenes <- c("PODXL", "NPHS2", "MME", "APOD",
                      "GPX3", "ANPEP", "CPE", "VCAM1", "AKR1B1", "CLDN10", "UMOD", "SLC12A1", 
                      "SMS", "KCNJ1", "STC1", "FXYD4", "AQP2", "ATP6V0D2",  "FOXI1", "MT1G", "MT2A", "IFIT1", "IFIT3", "MMP7", "LCN2", 
                      "IL32", "CD40", "ITGA2", "LIF", "TXNRD1", "GCLM", "SLC7A11", "HSPA6",
                      "MKI67", "UBE2S", "PAX2", "PAX8" , "SOX9", "SOX4", "PLCG2", "HES1", "KLF4", "SYNE2", "JUN",  "EPCAM", "CRYAB")

### creating dotplot ----
b <- DotPlot(RENAL, features= RENALMarkerGenes, col.min = 0, dot.scale = 5, group.by = "renalcelltype") + 
  scale_y_discrete(limits = rev(levels(RENAL$renalcelltype)), labels = paste(rev(levels(RENAL$renalcelltype_nr)), rev(levels(RENAL$renalcelltype)))) +
  scale_x_discrete(position = "top")  +
  theme_ipsum() +
  theme(axis.text.y = element_text(size = 16, color = rev(renalcelltype_colors), face = "bold", hjust = 0), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 16, face = "bold", angle = 90, hjust= 0.05, vjust=0.2),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position='right') +
  scale_color_viridis(name="avg. expr.", option = "viridis", direction = -1) +
  annotate("rect", xmin = 0.5, xmax = 4.5, ymin = 21.5, ymax = 23.5,
           alpha = .0, color = paletteer_d("ggsci::purple_material")[c(5)], linetype = "dashed", size = 1.2) +
  annotate("rect", xmin = 4.5, xmax = 19.5, ymin = 13.5, ymax = 21.5,
           alpha = .0, color = paletteer_d("ggsci::light_blue_material")[c(7)], linetype = "dashed", size = 1.2) +
  annotate("rect", xmin = 19.5, xmax = 29.5, ymin = 8.5, ymax = 13.5,
           alpha = .0, color = paletteer_d("ggsci::amber_material")[c(8)], linetype = "dashed", size = 1.2) +
  annotate("rect", xmin = 29.5, xmax = 33.5, ymin = 5.5, ymax = 8.5,
           alpha = .0, color = paletteer_d("ggsci::red_material")[c(5)], linetype = "dashed", size = 1.2) +  
  annotate("rect", xmin = 33.5, xmax = 35.5, ymin = 4.5, ymax = 5.5,
           alpha = .0, color =  paletteer_d("ggsci::pink_material")[4], linetype = "dashed", size = 1.2) +
  annotate("rect", xmin = 35.5, xmax = 44.5, ymin = 0.5, ymax = 4.5,
           alpha = .0, color = paletteer_d("ggsci::light_green_material")[c(7)], linetype = "dashed", size = 1.2) 

## compile plots and save ----
layout <- "
BBBB#CC
BBBB#CC
BBBB#DD
BBBB#DD
FFFFFFF
FFFFFFF
FFFFFFF
FFFFFFF
FFFFFFF
"

pw <-  p1  + p3 + p4 + b + plot_layout(design = layout, guides = "collect") 
for (i in 2:3) { pw[[i]] <- pw[[i]] + theme(title = element_text(size=16))}

ggsave(plot = pw, "supplfig6.png", h = 4000, w = 4400, units = "px", scale = 1, type = "cairo-png")

# Suppl. Fig. 7 ----

## Suppl. Fig. 7A ----
p <- FeaturePlot(RENAL, c("EPCAM", "CRYAB", "NPHS2", "PODXL", "LRP2", "CPE", "AQP1", "AKR1B1", "TNFSF10", "VCAM1", "SLC12A1", "UMOD", "DEFB1", "AQP2", "FXYD4", "KCNJ1", "SLC4A1", "ATP6V0D2"), pt.size =1, order =T, ncol=6) & 
  scale_color_viridis() &
  theme_void()


ggsave(plot = p, "SupplFig7a.png", w = 4000, h = 2000, units = "px", scale = 1, type = "cairo-png")

## Suppl. Fig. 7B ----
a <- FeaturePlot(RENAL, c("UMOD", "FXYD4"), blend=T, blend.threshold = 0.5, order=T, pt.size=1)
b <- FeaturePlot(RENAL, c("UMOD", "GPX3"), blend=T, blend.threshold = 0.5, order=T, pt.size=1)
c <- FeaturePlot(RENAL, c("GPX3", "FXYD4"), blend=T, blend.threshold = 0.5, order=T, pt.size=1)
d <- FeaturePlot(RENAL, c("PROM1", "VCAM1"), blend=T, blend.threshold = 0.5, order=T, pt.size=1)


ggsave(plot = a/b/c/d, "SupplFig7b.png", w = 4000, h = 4000, units = "px", scale = 1, type = "cairo-png")


# Suppl. Fig 8A, singleR automatic annotation with HPCA reference

# Here, we use singleR to automatically annotate our urine cells and use HumanPrimaryCellAtlas bulk data as reference. 
# We use the workflow for single-cell references described here:
# https://bioconductor.org/packages/devel/bioc/vignettes/SingleR/inst/doc/SingleR.html#1_Introduction

## subset RENAL object to inculde only TEC, not podocytes. ----
RENAL <- SetIdent(RENAL, value = "renalcelltype")
RENAL_TEC <- subset(RENAL, idents = c("PDC", "PDC_inj"), invert = T)

## load reference, subset to include only cell types of interest ----
ref <- HumanPrimaryCellAtlasData()
table(ref$label.main)
ctoi <- c("Embryonic_stem_cells", "Endothelial_cells", "Epithelial_cells", "Fibroblasts", "HSC_-G-CSF", "HSC_CD34+", "iPS_cells", "Keratinocytes", "MSC", "Smooth_muscle_cells", "Tissue_stem_cells")
sub_ref <- ref[,which(ref$label.main %in% ctoi)]

## retrieve data, lognormalize expression matrix for singleR ----
test <- GetAssayData(RENAL_TEC, assay = 'RNA', slot = 'data')
test <- LogNormalize(test)
 
## use singleR to compare HPCA data with our RENAL SO ----
singler_results_kidney_ref <- SingleR::SingleR(
  test = test,
  ref = sub_ref,
  labels = sub_ref$label.main
)

## annotate urine TEC in RENAL SO via new metadata $ref_labels ----
RENAL@meta.data$ref_labels <- singler_results_kidney_ref@listData$labels

## visualize in heatmap ----
RENAL.tbl <- RENAL_TEC@meta.data
RENAL.tbl$cell <- rownames(RENAL.tbl)
RENAL.tbl <- RENAL.tbl %>% arrange(renalcelltype)
RENAL.tbl$renalcelltype

names(renalcelltype_colors) <- levels(RENAL.tbl$renalcelltype)
annot.colors <- list(urine =  renalcelltype_colors)
annot.col <- data.frame(
  urine = RENAL.tbl$renalcelltype)
rownames(annot.col) <- rownames(RENAL.tbl)

p <- SingleR::plotScoreHeatmap(
  singler_results_kidney_ref,
  show.labels = F,
  order.by = "labels",
  annotation_col = annot.col,
  row.names = levels(singler_results_kidney_ref), 
  annotation_colors = annot.colors,
  fontsize = 18,
  annotation_legend = T,
  main = "Referencing of urine TEC with 'HumanPrimaryCellAtlasData'")

png("SupplFig8a.png", width=2000, height=1000, type="cairo-png")
p
dev.off()

# Suppl. Fig. 8B, reference annotation in UMAP ----

## identify annotated cells to highlight ----
RENAL_TEC <- SetIdent(RENAL_TEC, value="ref_labels")
esc <- WhichCells(RENAL_TEC, idents = c("Embryonic_stem_cells"))
tsc <- WhichCells(RENAL_TEC, idents = c("Tissue_stem_cells"))
msc <- WhichCells(RENAL_TEC, idents = c("MSC"))
ips <- WhichCells(RENAL_TEC, idents = c("iPS_cells"))
fib <- WhichCells(RENAL_TEC, idents = c("Fibroblasts"))
end <- WhichCells(RENAL_TEC, idents = c("Endothelial_cells"))
smc <- WhichCells(RENAL_TEC, idents = c("Smooth_muscle_cells"))
epc <- WhichCells(RENAL_TEC, idents = c("Epithelial_cells"))

## plot these cells in UMAP ----
a <- DimPlot(RENAL_TEC, cells.highlight=esc, label=F, pt.size = 1, sizes.highlight = 2)+theme_void()+ggtitle("Embryonic_stem_cells")
b <- DimPlot(RENAL_TEC, cells.highlight=tsc, label=F, pt.size = 1, sizes.highlight = 2)+theme_void()+ggtitle("Tissue_stem_cells")
c <- DimPlot(RENAL_TEC, cells.highlight=msc, label=F, pt.size = 1, sizes.highlight = 2)+theme_void()+ggtitle("Mesenchymal_stem_cells")
d <- DimPlot(RENAL_TEC, cells.highlight=ips, label=F, pt.size = 1, sizes.highlight = 2)+theme_void()+ggtitle("iPS_cells")
e <- DimPlot(RENAL_TEC, cells.highlight=fib, label=F, pt.size = 1, sizes.highlight = 2)+theme_void()+ggtitle("Fibroblasts")
f <- DimPlot(RENAL_TEC, cells.highlight=end, label=F, pt.size = 1, sizes.highlight = 2)+theme_void()+ggtitle("Endothelial_cells")
g <- DimPlot(RENAL_TEC, cells.highlight=smc, label=F, pt.size = 1, sizes.highlight = 2)+theme_void()+ggtitle("Smooth_muscle_cells")
h <- DimPlot(RENAL_TEC, cells.highlight=epc, label=F, pt.size = 1, sizes.highlight = 2)+theme_void()+ggtitle("Epithelial_cells")

## compile and save ----
png("SupplFig8b.png", width=2000, height=1000, type="cairo-png")
a+b+c+d+e+f+g+h+plot_layout(widths=c(1,1,1,1))&
  NoLegend()&
  scale_color_manual(values= c("lightgrey", "#558B2FFF"))&
  theme(plot.title = element_text(size=36, face="bold"))
dev.off()

# Suppl. Fig. 8C, violin plots ----

## compile and save ----
pw <- VlnPlot(RENAL_TEC, group.by="renalcellgroup", c("EPCAM", "CRYAB", "PLCG2", "SOX4", "HES1", "CITED2", "DDX17", "ELF3", "IER2"), ncol =1) & 
  scale_fill_manual(values=renalcellgroup_colors) 

for (i in c(1:4,6:9)) { pw[[i]] <- pw[[i]] + theme(axis.title.y = element_blank()
)}
for (i in c(1:8)) { pw[[i]] <- pw[[i]] + theme(axis.title.x = element_blank(),
                                               axis.ticks.x = element_blank(),
                                               axis.text.x = element_blank(),
                                               axis.text.y = element_text(size=24, face="bold"),
                                               plot.title = element_text(size=24, face="bold"))}
pw[[9]] <- pw[[9]] + theme(axis.title.x = element_blank(),
                           axis.text.x = element_text(size=24, face="bold"),
                           axis.text.y = element_text(size=24, face="bold"),
                           plot.title = element_text(size=24, face="bold"))
pw[[5]] <- pw[[5]] + ylab("Expression Level\n")+theme(axis.title.y = element_text(size=24, face="bold", vjust = 0))

png("SupplFig8c.png", width=800, height=2000, type="cairo-png")
pw
dev.off()


# Suppl. Fig. 17 ----

## subset URINE object to include only one cell subset ----
URINE@active.assay <- "RNA"
URINE <- NormalizeData(URINE)
URINE <- SetIdent(URINE, value="cellgroup")
URINE_PDC <- subset(URINE, idents = "PDC", invert=F)
URINE_TEC <- subset(URINE, idents = "TEC", invert=F)
URINE_LEUK <- subset(URINE, idents = "LEUK", invert=F)
URINE_UGEC <- subset(URINE, idents = "UGEC", invert=F)



## find DEG COVID vs. non-COVID in podocytes and TEC subsets ----
COVID.markers.PDC <- FindMarkers(URINE_PDC, ident.1= "yes", ident.2="no", group.by="COVID19_infection", logfc.threshold = 0.01)
COVID.markers.TEC <- FindMarkers(URINE_TEC, ident.1= "yes", ident.2="no", group.by="COVID19_infection", logfc.threshold = 0.01)
COVID.markers.LEUK <- FindMarkers(URINE_LEUK, ident.1= "yes", ident.2="no", group.by="COVID19_infection", logfc.threshold = 0.01)
COVID.markers.UGEC <- FindMarkers(URINE_UGEC, ident.1= "yes", ident.2="no", group.by="COVID19_infection", logfc.threshold = 0.01)
COVID.markers <- FindMarkers(URINE, ident.1= "yes", ident.2="no", group.by="COVID19_infection", logfc.threshold = 0.01)

## Suppl Fig 17A ----

### filter DEG for top DEG ----
subset.DE.top <- COVID.markers %>% filter(p_val < 1e-50) %>% top_n(n = 20, wt = avg_log2FC)
subset.DE.low <- COVID.markers %>% filter(p_val < 1e-50) %>% top_n(n = -20, wt = avg_log2FC)

### volcano plot ----
v <- EnhancedVolcano(COVID.markers,
                     lab = rownames(COVID.markers),
                     x = 'avg_log2FC',
                     y = 'p_val_adj',
                     title = "DEG non-COVID AKI vs. COVID AKI",
                     selectLab = c(rownames(subset.DE.top), rownames(subset.DE.low)),
                     xlab = bquote(~Log[2]~ 'fold change'),
                     ylab = bquote(~-Log[10]~ '(adj. p)'),
                     pCutoff = 10e-50,
                     FCcutoff = 0.5,
                     pointSize = 1.0,
                     labSize = 3.0,
                     labCol = 'black',
                     boxedLabels = F,
                     colAlpha = 3/5,
                     legendPosition = 'right',
                     legendLabSize = 10,
                     legendIconSize = 2.0,
                     drawConnectors = TRUE,
                     widthConnectors = 1.0,
                     colConnectors = 'black')

### save ----

ggsave(plot = v, "supplfig17a.png", w = 3000, h = 2000, units = "px", scale = 1, type = "cairo-png")

## Suppl Fig 17B ----

### filter for top COVID DEG ----
PDC <- COVID.markers.PDC %>% 
  mutate(gene = rownames(.)) %>%
  filter(p_val_adj < 10e-10) %>%
  slice_max(avg_log2FC, n=10) %>%
  arrange(avg_log2FC)
TEC <- COVID.markers.TEC %>% 
  mutate(gene = rownames(.)) %>%
  filter(p_val_adj < 10e-10) %>%
  slice_max(avg_log2FC, n=10) %>%
  arrange(avg_log2FC)
LEUK <- COVID.markers.LEUK %>% 
  mutate(gene = rownames(.)) %>%
  filter(p_val_adj < 10e-10) %>%
  slice_max(avg_log2FC, n=10) %>%
  arrange(avg_log2FC)
UGEC <- COVID.markers.UGEC %>% 
  mutate(gene = rownames(.)) %>%
  filter(p_val_adj < 10e-10) %>%
  slice_max(avg_log2FC, n=10) %>%
  arrange(avg_log2FC)

### Barplots of top COVID DEG ----
a <- ggplot(TEC, aes(x=gene, y=avg_log2FC))+
  geom_bar(stat = "identity", fill = "blue")+
  scale_x_discrete(limits = TEC$gene)+
  ggtitle("TEC")+
  theme_classic()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(face= "bold"))+
  coord_flip()
b <- ggplot(PDC, aes(x=gene, y=avg_log2FC))+
  geom_bar(stat = "identity", fill = "blue")+
  scale_x_discrete(limits = PDC$gene)+
  ggtitle("Podocytes")+
  theme_classic()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(face= "bold"))+
  coord_flip()
c <-ggplot(LEUK, aes(x=gene, y=avg_log2FC))+
  geom_bar(stat = "identity", fill = "blue")+
  scale_x_discrete(limits = LEUK$gene)+
  ggtitle("Immune cells")+
  theme_classic()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(face= "bold"))+
  coord_flip()
d <- ggplot(UGEC, aes(x=gene, y=avg_log2FC))+
  geom_bar(stat = "identity", fill = "blue")+
  scale_x_discrete(limits = UGEC$gene)+
  ggtitle("UGEC")+ 
  theme_classic()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(face= "bold"))+
  coord_flip()

### compile and save ----

ggsave(plot = a+b+c+d+plot_layout(widths = c(1,1,1,1)), "supplfig17b.png", w = 3000, h = 700, units = "px", scale = 1, type = "cairo-png")
