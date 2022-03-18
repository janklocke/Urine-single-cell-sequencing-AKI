# description ----

# This is the R script for generating Figure 1 and Supplemental Figure 4 in 
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

## loading data ----
URINE <- readRDS("~/SO_all_urine_cells.rds")
URINE0 <- readRDS("~/SO_all_urine_cells_noQC.rds")

## loading color schemes ----
celltype_colors <- colorRampPalette(paletteer_d("palettesForR::Pastels"))(19)
cellgroup_colors <- brewer.pal(5, "Paired")
AKI_type_colors <- paletteer_d("unikn::pal_signal")
names(AKI_type_colors) <- levels(factor(URINE$AKI_type))[c(1,3,2)]
sampling_colors <- paletteer_d("ggsci::alternating_igv")
names(sampling_colors) <- levels(factor(URINE$sampling))
gating_colors <- paletteer_d("nbapalettes::pacers_classic")
names(gating_colors) <- levels(factor(URINE$gating))

custom_colors <- list()
colors_dutch <- c(
  '#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67',
  '#F79F1F','#A3CB38','#1289A7','#D980FA','#B53471',
  '#EE5A24','#009432','#0652DD','#9980FA','#833471',
  '#EA2027','#006266','#1B1464','#5758BB','#6F1E51')
colors_spanish <- c(
  '#40407a','#706fd3','#f7f1e3','#34ace0','#33d9b2',
  '#2c2c54','#474787','#aaa69d','#227093','#218c74',
  '#ff5252','#ff793f','#d1ccc0','#ffb142','#ffda79',
  '#b33939','#cd6133','#84817a','#cc8e35','#ccae62')
custom_colors$discrete <- c(colors_dutch, colors_spanish)
custom_colors$cell_cycle <- setNames(
  c('#45aaf2', '#f1c40f', '#e74c3c', '#7f8c8d'),
  c('G1',      'S',       'G2M',     '-'))

## functions ----
number_ticks <- function(n) {function(limits) pretty(limits, n)}

# Dimensional plots (UMAP) (Figure 1A-C, Supplemental Fig. 4 A-C)----

## Urine cell clusters (Figure 1A) ----
p1 <- DimPlot(URINE, reduction = "umap", label = T, repel = F, label.size = 10, group.by = "celltype_nr", cols = celltype_colors) +
  labs(title = "Urine cell clusters") +
  theme_void() +
  theme(plot.title = element_text(face= "bold", size = 20, hjust=0.5, vjust = 1)) + 
  NoLegend()

## Gender distribution (Figure 1B) ----
p2 <- DimPlot(URINE, reduction = "umap", shuffle = T, label.size = 4, group.by = "gender") +
  theme_void() + 
  scale_color_manual(values = c("red", "lightblue", "lightgrey")) + 
  ggtitle("Sex") + 
  theme(plot.title = element_text(face= "bold", size = 18, hjust=0.5)) + 
  NoLegend()

## Distribution by disease type (Figure 1C) ----
URINE <- SetIdent(URINE, value = "AKI_type")
p5 <- DimPlot(URINE, reduction = "umap", shuffle = T, group.by = "AKI_type", cols = AKI_type_colors) +
  theme_void() + 
  scale_color_manual(values=AKI_type_colors)+
  ggtitle("AKI type") + 
  theme(plot.title = element_text(face= "bold", size = 18, hjust=0.5),
        legend.text = element_text(face= "bold", size = 10), 
        legend.title = element_text(face= "bold", size = 10), 
        legend.position='topleft') + 
  NoLegend()

## Distribution of broad cell groups (Suppl. Fig. 4 B) ----
p3 <- DimPlot(URINE, reduction = "umap", shuffle = T, label.size = 5, label = F, repel = T, group.by = "cellgroup") +
  theme_void() + 
  scale_color_brewer(palette = "Paired") +
  ggtitle("Celltype") + 
  theme(plot.title = element_text(face= "bold", size = 30, hjust=0.5),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12), 
        legend.position='right') + 
  NoLegend()

## Distribution by single samples (Suppl. Fig. 4C)
p4 <- DimPlot(URINE, reduction = "umap", shuffle = T, label.size = 5, label = F, repel = T, group.by = "orig.ident") +
  theme_void() + 
  ggtitle("Sample") + 
  theme(plot.title = element_text(face= "bold", size = 30, hjust=0.5)) + 
  NoLegend()

# Dotplots of marker gene expression (Figure 1D) ----

## assigning marker genes ----
URINEMarkerGenes <- c("NPHS2", "PODXL", "GPX3", "ALDOB", "CRYAB", "EPCAM", "AQP2", "FXYD4", "MMP7", "IGFBP7", 
                      "LAMC2", "GCLM", "TXNRD1", "MKI67", "UBE2S", "SOX4", "PLCG2",
                      "TIMP1", "IL1B", "C1QB", "CD74", "MT1G", "MT1H", "HSPH1", "SLC11A1",
                      "SPP1", "ALOX5AP", "MNDA", "CSF3R", "CD3E", "NKG7", "CD79A", "IGKC", 
                      "HBB", "HBA1", "KRT13", "PSCA", "UPK2")

## creating dotplot ----
p <- DotPlot(URINE, features = unique(URINEMarkerGenes), group.by = "celltype", col.min = 0, dot.scale = 5) + 
  scale_y_discrete(limits = rev(levels(URINE$celltype)), labels = paste(rev(levels(URINE$celltype_nr)), rev(levels(URINE$celltype)))) +
  scale_x_discrete(position = "top")  +
  theme_ipsum() +
  theme(axis.text.y = element_text(size = 16, color = rev(celltype_colors), face = "bold", hjust = 0), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 16, face = "bold", angle = 90, hjust= 0.05, vjust=0.2),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position='right') +
  scale_colour_viridis(name = 'avg. expr.', option = "viridis", direction = -1) +
  annotate("rect", xmin = 0.5, xmax = 2.5, ymin = 18.5, ymax = 19.5,
           alpha = .0, color = cellgroup_colors[3], linetype = "dashed", size = 1.2) +
  annotate("rect", xmin = 2.5, xmax = 17.5, ymin = 11.5, ymax = 18.5,
           alpha = .0, color = cellgroup_colors[4], linetype = "dashed", size = 1.2) +
  annotate("rect", xmin = 17.5, xmax = 33.5, ymin = 3.5, ymax = 11.5,
           alpha = .0, color = cellgroup_colors[2], linetype = "dashed", size = 1.2) +
  annotate("rect", xmin = 33.5, xmax = 35.5, ymin = 2.5, ymax = 3.5,
           alpha = .0, color = cellgroup_colors[1], linetype = "dashed", size = 1.2) +
  annotate("rect", xmin = 35.5, xmax = 38.5, ymin = 0.5, ymax = 2.5,
           alpha = .0, color = cellgroup_colors[5], linetype = "dashed", size = 1.2) 

# Barplots of AKI etiology and gender by cluster (Figure 1D) ----

## Barplot of AKI etiology by cluster ----
p6 <- URINE@meta.data%>% group_by(celltype, AKI_type) %>%
  filter(AKI_type != "na") %>%
  summarise(count = n()) %>%
  spread(AKI_type, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('celltype', 'total_cell_count', everything())) %>%
  arrange(factor(celltype, levels = levels(URINE$celltype))) %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'celltype') %>%
  ggplot(aes(celltype, value)) +
  geom_bar(aes(fill = variable), position = 'fill', stat = 'identity') +
  scale_x_discrete(limits = rev(levels(URINE$celltype))) + 
  coord_cartesian(clip = 'off') +
  scale_fill_manual(name = "AKI type", values = c(AKI_type_colors[c(1,3,2)], "lightgrey"), limits = c("CS", "pneumonia", "prerenal", "na"))+
  coord_flip() + 
  theme_void() 

## Barplot of gender by cluster ----
p7 <- URINE@meta.data%>% group_by(celltype, gender) %>%
  filter(gender != "na") %>%
  summarise(count = n()) %>%
  spread(gender, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('celltype', 'total_cell_count', everything())) %>%
  arrange(factor(celltype, levels = levels(URINE$celltype))) %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'celltype') %>%
  ggplot(aes(celltype, value)) +
  geom_bar(aes(fill = variable), position = 'fill', stat = 'identity') +
  scale_x_discrete(limits = rev(levels(URINE$celltype))) + 
  coord_cartesian(clip = 'off') +
  scale_fill_manual(name = "sex", values = c("red", "lightblue", "lightgrey")) + 
  coord_flip() + 
  theme_void()

p6 <- p6+ggtitle("AKI type") + theme(plot.title=element_text(size= 18, angle=90))
p7 <- p7+ggtitle("sex") + theme(plot.title=element_text(size= 18, angle=90))

# Assembly of Figure 1 plots ----

## Figure 1 A-C ----
layout <- "
AAB
AAC
"
pw <- p1+p2+p5+plot_layout(design = layout, guides = "collect")
pw <- pw & theme(legend.text = element_text(size=12),
                 legend.title = element_text(size=24))
ggsave(plot = pw, "fig1A.png", w = 4000, h = 2000, units = "px", scale = 1, type = "cairo-png")

## Figure 1 D ----
layout2 <- 
  "BBBBBBBBBBBBBBCD
BBBBBBBBBBBBBBCD
BBBBBBBBBBBBBBCD
BBBBBBBBBBBBBBCD
BBBBBBBBBBBBBBCD
BBBBBBBBBBBBBBCD"
qw <- p+p6+p7+plot_layout(design = layout2, guides = "collect")
qw <- qw & theme(legend.text = element_text(size=12),
                 legend.title = element_text(size=12))
ggsave(plot = qw, "fig1D.png", w = 4200, h = 2400, units = "px", scale = 1, type = "cairo-png")

# Bar- and violinplots for Suppl. Fig. 4D ----

## Create a meta.data table of URINE and URINE0 ----

URINE0_tbl <- URINE0@meta.data
URINE_tbl <- URINE@meta.data

## bar- and violinplots of Suppl. Fig. 4 plots ----

### violin plot nFeatures ----

q1 <- ggplot(URINE0_tbl, aes(x=patient, y=nFeature_RNA, fill = patient)) +
  geom_violin(scale = 'width', trim = T) +
  theme_classic() + 
  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, face= "bold", size=10),
        axis.title.x = element_text(face= "bold", size=18),
        axis.text.y = element_text(vjust = 0.5, hjust=0, size=18, face= "bold"), 
        axis.title.y = element_blank()) + 
  scale_y_continuous(name = "Mean\n Genes", breaks=number_ticks(2), labels = scales::label_number_si()) +
  scale_color_manual(values = c(custom_colors[[1]], custom_colors[[2]])) +
  NoLegend() + 
  geom_hline(yintercept = 200, linetype = "solid", color = "red", size = 1, alpha = 0.2) + 
  geom_hline(yintercept = 4000, linetype = "solid", color = "red", size = 1, alpha = 0.2)+ 
  coord_flip()

### violin plot nCount ----

q2 <- ggplot(URINE0_tbl, aes(x=patient, y=nCount_RNA, fill = patient)) +
  geom_violin(scale = 'width', trim = T) +
  theme_classic() + 
  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, face= "bold", size=10),
        axis.title.x = element_text(face= "bold", size=18),
        axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5), 
        axis.title.y = element_text(face= "bold", size=12)) + 
  scale_y_continuous(name = "Mean \n Transcripts", breaks=number_ticks(2), limits=c(0,150000), labels = scales::label_number_si()) +
  scale_color_manual(values = c(custom_colors[[1]], custom_colors[[2]])) +
  NoLegend()+ 
  geom_hline(yintercept = 500, linetype = "solid", color = "red", size = 1, alpha = 0.2) + 
  geom_hline(yintercept = 50000, linetype = "solid", color = "red", size = 1, alpha = 0.2)+ 
  coord_flip()

### violin plot percent mitochondrial RNA ----

q3 <- ggplot(URINE0_tbl, aes(x=patient, y=percent.mt, fill = patient)) +
  geom_violin(scale = 'width', trim = T) +
  theme_classic() + 
  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, face= "bold", size=10),
        axis.title.x = element_text(face= "bold", size=18),
        axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5), 
        axis.title.y = element_text(face= "bold", size=12)) + 
  scale_y_continuous(name = "% Mt \n Transcripts", breaks=number_ticks(3), limits=c(0,100)) +
  scale_color_manual(values = c(custom_colors[[1]], custom_colors[[2]])) +
  NoLegend()+ 
  geom_hline(yintercept = 20, linetype = "solid", color = "red", size = 1, alpha = 0.2) + 
  geom_hline(yintercept = 10, linetype = "longdash", color = "red", size = 1, alpha = 0.2)+ 
  coord_flip()

### bar plot cell type ----

q4 <- URINE_tbl%>% group_by(patient, celltype) %>%
  summarise(count = n()) %>%
  spread(celltype, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('patient', 'total_cell_count', everything())) %>%
  arrange(factor(patient, levels = levels(URINE$patient))) %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'patient') %>%
  ggplot(aes(patient, value)) +
  geom_bar(aes(fill = variable), position = 'fill', stat = 'identity') +
  scale_fill_manual(name = 'Cluster', values = celltype_colors) +
  scale_y_continuous(name = 'Cluster \n proportion', breaks=number_ticks(1)) +
  coord_cartesian(clip = 'off') +
  theme_classic() +
  theme(
    legend.position = 'left',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(vjust = 0.5, hjust=1, face= "bold", size=10),
    axis.title.x = element_text(face= "bold", size=18),
    axis.text.y = element_text(vjust = 0.5, hjust=0, face= "bold", size=12), 
    axis.title.y = element_blank(),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
  ) +
  NoLegend()+ 
  coord_flip()

### bar plot cell group ----

q5 <- URINE_tbl%>% group_by(patient, cellgroup) %>%
  summarise(count = n()) %>%
  spread(cellgroup, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('patient', 'total_cell_count', everything())) %>%
  arrange(factor(patient, levels = levels(URINE$patient))) %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'patient') %>%
  ggplot(aes(patient, value)) +
  geom_bar(aes(fill = variable), position = 'fill', stat = 'identity') +
  scale_fill_brewer(palette = "Paired", labels = c("ERY", "LEUK", "PDC", "TEC", "UGEC")) +
  scale_y_continuous(name = 'Celltype \n', breaks=number_ticks(1)) +
  coord_cartesian(clip = 'off') +
  theme_classic() +
  theme(
    legend.position = 'left',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(vjust = 0.5, hjust=1, face= "bold", size=16),
    axis.title.x = element_text(face= "bold", size=18),
    axis.text.y = element_text(vjust = 0.5, hjust=0, face= "bold", size=16), 
    axis.title.y = element_blank(), 
    plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
  ) +
  coord_flip()

### bar plot gender ----

q6 <- URINE_tbl%>% group_by(patient, gender) %>%
  summarise(count = n()) %>%
  spread(gender, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('patient', 'total_cell_count', everything())) %>%
  arrange(factor(patient, levels = levels(URINE$patient))) %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'patient') %>%
  ggplot(aes(patient, value)) +
  geom_bar(aes(fill = variable), position = 'fill', stat = 'identity') +
  scale_y_continuous(name = 'sex \n') +
  scale_fill_manual(name = 'sex', values = c("red", "lightblue", "lightgrey")) + 
  coord_cartesian(clip = 'off') +
  theme_classic() +
  theme(
    legend.position = 'right',
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face= "bold", size=14, angle = 90, vjust = 0.5),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
  ) +
  coord_flip()

### bar plot AKI etiology ----

q7 <- URINE_tbl %>% group_by(patient, AKI_type) %>%
  summarise(count = n()) %>%
  spread(AKI_type, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('patient', 'total_cell_count', everything())) %>%
  arrange(factor(patient, levels = levels(URINE$patient))) %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'patient') %>%
  ggplot(aes(patient, value)) +
  geom_bar(aes(fill = variable), position = 'fill', stat = 'identity') +
  scale_fill_manual(name = 'AKI type', values = AKI_type_colors, labels=c("CS", "prerenal", "pneumonia")) +
  scale_y_continuous(name = 'AKI type') +
  coord_cartesian(clip = 'off') +
  theme_classic() +
  theme(
    legend.position = 'right',
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face= "bold", size=14, angle = 90, vjust = 0.5),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
  ) +
  coord_flip()

### bar plot sampling ----

q8 <- URINE_tbl%>% group_by(patient, sampling) %>%
  summarise(count = n()) %>%
  spread(sampling, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('patient', 'total_cell_count', everything())) %>%
  arrange(factor(patient, levels = levels(URINE$patient))) %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'patient') %>%
  ggplot(aes(patient, value)) +
  geom_bar(aes(fill = variable), position = 'fill', stat = 'identity') +
  scale_fill_manual(name = 'Pooling', values =  sampling_colors) +
  scale_y_continuous(name = 'pooled') +
  coord_cartesian(clip = 'off') +
  theme_classic() +
  theme(
    legend.position = 'right',
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face= "bold", size=14, angle = 90, vjust = 0.5),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
  ) +
  coord_flip()

### bar plot gating ----

q9 <- URINE_tbl%>% group_by(patient, gating) %>%
  summarise(count = n()) %>%
  spread(gating, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('patient', 'total_cell_count', everything())) %>%
  arrange(factor(patient, levels = levels(URINE$patient))) %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'patient') %>%
  ggplot(aes(patient, value)) +
  geom_bar(aes(fill = variable), position = 'fill', stat = 'identity') +
  scale_fill_manual(name = 'Gating', values = gating_colors) +
  scale_y_continuous(name = 'gating') +
  coord_cartesian(clip = 'off') +
  theme_classic() +
  theme(
    legend.position = 'right',
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face= "bold", size=14, angle = 90, vjust = 0.5),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
  ) + 
  coord_flip()

### bar plot doublets ----

q10 <- URINE0_tbl%>% group_by(patient, multiplet_class) %>%
  summarise(count = n()) %>%
  spread(key = multiplet_class, value = count) %>%
  mutate(perc_doublet = (doublet/(doublet + singlet)*100)) %>%
  replace_na(list(perc_doublet = 0)) %>%
  ggplot(aes(x=patient, y=perc_doublet, fill = patient)) +
  geom_bar(aes(fill = patient), stat = 'identity') +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 10, face = "bold", vjust = 0.5, hjust=1),
        axis.text.y = element_text(vjust = 0.5, hjust=0.5), 
        axis.title.y = element_text(face= "bold", size=12), 
        axis.title.x = element_text(face= "bold", size=18)) + 
  scale_y_continuous(name = "% doublets", breaks=number_ticks(3)) +
  scale_x_discrete(name = "Sample") +
  coord_cartesian(clip = 'off') +
  NoLegend()+ 
  coord_flip(clip = 'off')

### bar plot single-cell count ----

q11 <- URINE_tbl%>% group_by(patient) %>%
  summarise(count = n()) %>%
  ggplot(aes(x=patient, y=count, fill = patient)) +
  geom_bar(aes(fill = patient), stat = 'identity') +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 10, face = "bold", vjust = 0.5, hjust=1),
        axis.text.y = element_text(vjust = 0.5, hjust=0.5), 
        axis.title.y = element_text(face= "bold", size=12), 
        axis.title.x = element_text(face= "bold", size=18)) + 
  scale_y_continuous(name = "single cells \n post-QC", labels = scales::label_number_si()) +
  scale_x_discrete(name = "Sample") +
  coord_cartesian(clip = 'off') +
  NoLegend()+ 
  coord_flip()


## Assembly of Suppl. Fig. 4 plots ----

layout3 <- "
LLLLLLLLLLMMMMMMMMMMNNNNNNNNNN##
LLLLLLLLLLMMMMMMMMMMNNNNNNNNNN##
LLLLLLLLLLMMMMMMMMMMNNNNNNNNNN##
LLLLLLLLLLMMMMMMMMMMNNNNNNNNNN##
LLLLLLLLLLMMMMMMMMMMNNNNNNNNNN##
LLLLLLLLLLMMMMMMMMMMNNNNNNNNNN##
AAAABBBBCCCCDDDDEEEEFGHIJJJJKKKK
AAAABBBBCCCCDDDDEEEEFGHIJJJJKKKK
AAAABBBBCCCCDDDDEEEEFGHIJJJJKKKK
AAAABBBBCCCCDDDDEEEEFGHIJJJJKKKK
AAAABBBBCCCCDDDDEEEEFGHIJJJJKKKK
AAAABBBBCCCCDDDDEEEEFGHIJJJJKKKK
AAAABBBBCCCCDDDDEEEEFGHIJJJJKKKK
AAAABBBBCCCCDDDDEEEEFGHIJJJJKKKK
AAAABBBBCCCCDDDDEEEEFGHIJJJJKKKK
AAAABBBBCCCCDDDDEEEEFGHIJJJJKKKK
AAAABBBBCCCCDDDDEEEEFGHIJJJJKKKK
AAAABBBBCCCCDDDDEEEEFGHIJJJJKKKK
"

rw <-q1 + q2 + q3 + q4 + q5 + q6 + q7 + q8 + q9 + q10 + q11 + p2 + p3 + p4 + plot_layout(design = layout3, guides = "collect")
for (i in 2:11) { rw[[i]] <- rw[[i]] + theme(axis.text.y = element_blank(),
                                             axis.title.y = element_blank() ) }
for (i in c(1:5,10:11)) { rw[[i]] <- rw[[i]] + theme(axis.text.x = element_text(size = 20, face = "bold", vjust = 0.5, hjust=0.5),
                                                     axis.title.x = element_text(size = 24)) }
for (i in c(6:9)) { rw[[i]] <- rw[[i]] + theme(axis.title.x = element_text(size = 20)) }

png("supplfig1.png", width = 2000, height = 1500, type="cairo-png")
rw & theme(legend.text = element_text(size=24),
           legend.title = element_text(size=24))
dev.off()
