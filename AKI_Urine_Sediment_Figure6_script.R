# description ----

# This is the R script for generating Figure 6 in 
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
library(ggpubr)

## loading data ----
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

renalcellgroup_colors <- c(renalcelltype_colors[c(1,8,14,17,2,21)], "grey")
names(renalcellgroup_colors) <- c(levels(factor(RENAL$renalcellgroup)), "total_cell_count")

## functions ----
number_ticks <- function(n) {function(limits) pretty(limits, n)}


# Fig 6AB ----

## functions for Fig6AB ----
## comparison_plots_time() creates a boxplot of absolute counts of a chosen kidney cell type per patient and
## compares them over three timeframes (d0-5, d6-10, d11-21):
comparison_plots_time <- function(SO, xa, xn, ygroup, ya) {
  groups <- levels(factor(SO@meta.data[[xn]]))
  my_comparisons <- list(groups[1:2] ,groups[2:3] , groups[c(1,3)])
  a <-SO@meta.data %>% group_by({{ygroup}}, patient, {{xa}}) %>%
    dplyr::summarize(count = n()) %>%
    spread({{ygroup}}, count, fill = 0) %>%
    ungroup() %>%
    mutate(total_cell_count = rowSums(.[c(3:ncol(.))])) %>%
    filter({{xa}} != "pool") %>%
    ggplot(aes(x={{xa}}, y={{ya}}, color={{xa}}, fill={{xa}})) +
    geom_boxplot(alpha=0.7, color="black", size=0.25) +
    geom_jitter(size=1, pch=21, stroke=0.25, color= "black", alpha = 0.8, width =0.25, height=0)+
    scale_y_log10(name= "abs. cell count", limits=c(1,50000))+
    scale_x_discrete(name= "time sections after AKI onset", limits=c("d0-5", "d6-10", "d11-21"))+
    scale_color_manual(values=rep(renalcellgroup_colors[[deparse(substitute(ya))]],3))+
    scale_fill_manual(values=rep(renalcellgroup_colors[[deparse(substitute(ya))]],3))+
    theme_minimal()+
    theme(legend.position = "bottom", 
          plot.title = element_text(size=16, vjust=0.5, hjust=0.5, face="bold"),
          axis.title.y = element_text(size=16,hjust=0.5), 
          axis.title.x =  element_text(size=16,hjust=0.5, vjust=-1), 
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10), 
          axis.text = element_text(vjust = 0.5, hjust=0.5, face= "bold", size=12))+
    stat_compare_means(comparisons = my_comparisons,  label = "p.signif", size= 4, stroke=2, hide.ns=F)
  return(a)
}

## group_barB() generates a geom_area plot for displaying time shifts in kidney cell type proportions. 
## patients who contribute less than 50 cells to the analysis are excluded. 
group_barB <- function(SO) { 
  p <- SO@meta.data %>% group_by(timepoint2, renalcellgroup) %>%
    summarise(count = n()) %>%
    filter(renalcellgroup != "PDC") %>%
    spread(renalcellgroup, count, fill = 0) %>%
    ungroup() %>%
    mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
    mutate_at(colnames(.[c(2:(ncol(.)-1))]), function(x) x/.$total_cell_count) %>%
    filter(total_cell_count > 49, timepoint2 < 13) %>%
    dplyr::select(c('timepoint2', 'total_cell_count', everything())) %>%
    arrange(factor(timepoint2, levels = levels(SO$timepoint2))) %>%
    dplyr::select(-c('total_cell_count')) %>%
    reshape2::melt(id.vars = 'timepoint2') %>%
    ungroup()
  p$timepoint <- as.numeric(p$timepoint)
  
  p <-ggplot(p, aes(x=timepoint2, y=value, fill=variable)) +
    geom_area(colour= "black", alpha = 0.6, size = 0.25) +
    ylab("cell type\nproportion [%]") +
    xlab("time since AKI onset [d]") +
    scale_y_continuous(breaks=number_ticks(1)) +
    scale_x_continuous(breaks = seq(from = 0, to = 12, by = 4))+
    scale_fill_manual(values = renalcellgroup_colors, name = "cell type") +
    theme_minimal() +
    ggtitle("relative celltype proportions") +
    theme(legend.position = "right", 
          plot.title = element_text(size=16, vjust=0.5, hjust=0.5, face="bold"),
          axis.title.y = element_text(size=16,hjust=0.5), 
          axis.title.x = element_text(size=16,hjust=0.5, vjust=-1), 
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10), 
          axis.text = element_text(vjust = 0.5, hjust=0.5, face= "bold", size=12))+
    geom_segment(x=0, y=0, yend=1, xend=0, size=0.25) +
    geom_segment(x=12, y=0, yend=1, xend=12, size=0.25) +
    geom_segment(x=0, y=0, yend=0, xend=12, size=0.25)
  
  return(p)
}

## preparation of RENAL seurat object ----

## To visualize an overall trend of shifting TEC cell proportions over time, patients are binned in 3d intervals for Fig 6B. 
RENAL$timepoint2 <- NA
RENAL@meta.data[RENAL$AKI_timepoint.fine %in% c("0", "-2"),]$timepoint2 <- 0
RENAL@meta.data[RENAL$AKI_timepoint.fine %in% c("1", "2", "3"),]$timepoint2 <- 3
RENAL@meta.data[RENAL$AKI_timepoint.fine %in% c("4", "5", "6"),]$timepoint2 <- 6
RENAL@meta.data[RENAL$AKI_timepoint.fine %in% c("7", "8", "9"),]$timepoint2 <- 9
RENAL@meta.data[RENAL$AKI_timepoint.fine %in% c("10", "11", "12"),]$timepoint2 <- 12
RENAL@meta.data[RENAL$AKI_timepoint.fine %in% c("13", "14", "16"),]$timepoint2 <- 15

## Patients sample sizes are downsampled to max. 100 cells/sample on order to account for varying sample sizes 
RENAL <- SetIdent(RENAL, value= "patient")
RENAL_sub <- subset(RENAL, downsample = 100)

group_barB(RENAL_sub)

## create plots, compile and save ----
pw <- comparison_plots_time(RENAL, AKI_timepoint.main, "AKI_timepoint.main", renalcellgroup, total_cell_count)+ggtitle("total")+
  comparison_plots_time(RENAL, AKI_timepoint.main, "AKI_timepoint.main", renalcellgroup, TEC_emt)+ggtitle("TEC_emt")+
  comparison_plots_time(RENAL, AKI_timepoint.main, "AKI_timepoint.main", renalcellgroup, TEC_str)+ggtitle("TEC_str")+
  comparison_plots_time(RENAL, AKI_timepoint.main, "AKI_timepoint.main", renalcellgroup, TEC_prlf)+ggtitle("TEC_prlf")+
  comparison_plots_time(RENAL, AKI_timepoint.main, "AKI_timepoint.main", renalcellgroup, TEC_prg)+ggtitle("TEC_prg")+
  group_barB(RENAL_sub)+
  plot_layout(widths=c(1,1,1,1,1,1.5))

for (i in c(2:5)) {pw[[i]] <- pw[[i]] + theme(axis.title.y = element_blank(),
                                              axis.text.y = element_blank()) }
for (i in c(1:2,4:5)) {pw[[i]] <- pw[[i]] + theme(axis.title.x = element_blank()) }

ggsave(plot = pw &NoLegend(), "fig6ab.png", w = 4200, h = 870, units = "px", scale = 1, type = "cairo-png")

# Fig. 6C, individual patients ----

## subset RENAL object to individual patients with multiple measurements ----
RENAL <- SetIdent(RENAL, value = "patient")
RENAL_P17 <- subset(RENAL, ident = c("P017.1", "P017.2"))
RENAL_P18 <- subset(RENAL, ident = c("P018.1", "P018.2"))
RENAL_P19 <- subset(RENAL, ident = c("P019.1", "P019.2"))
RENAL_P23 <- subset(RENAL, ident = c("P023.1", "P023.2", "P023.3", "P023.4"))
RENAL_P24 <- subset(RENAL, ident = c("P024.1", "P024.2"))
RENAL_P02 <- subset(RENAL, ident = c("P002.1", "P002.2"))

## functions for Fig 6C ----

## UC_plot_uo creates a plot for urinary output after AKI onset for a single patient based on the data provided in UC_data.df 
UC_plot_uo <- function(patientname){
  
  p1 <- UC_data.df %>% filter(patient == patientname) %>% 
    ggplot(aes(x=day, y=urinary_output)) +
    geom_col(width = 1, size=5, fill="#69b3a2", alpha=0.6) +
    theme_minimal() +
    scale_y_continuous(name = "urinary\noutput\n[ml]\n", limits = c(0, 6500), labels = scales::label_number_si()) +
    geom_hline(yintercept = 2000, linetype = "solid", color = "black", size = 1, alpha = 0.1) + 
    xlim(min(UC_data.df$day), 13) +
    theme(axis.title.y = element_text(color = "#69b3a2", size=16), 
          axis.title.x = element_text(size=16), 
          axis.text = element_text(vjust = 0.5, hjust=0, face= "bold", size=12))
  
  return(p1)}

## UC_plot_cr creates a plot for serum creatinine levels after AKI onset for a single patient based on the data provided in UC_data.df 
UC_plot_cr <- function(patientname){
  
  p2<- UC_data.df %>% filter(patient == patientname) %>% 
    ggplot(aes(x=day, y=serum_creatinine)) +
    geom_area(fill="lightgoldenrod3", alpha=0.5) +
    geom_line(color="lightgoldenrod3", size=1) +
    geom_point(color="lightgoldenrod3", size=1) +
    theme_minimal()+
    scale_y_continuous(name = "serum\ncreatinine\n[mg/dl]\n", limits = c(0, 5)) +
    xlim(min(UC_data.df$day), 13) +
    geom_hline(yintercept = 1, linetype = "solid", color = "black", size = 1, alpha = 0.1) +
    labs(title = paste0("P0", str_remove_all(patientname, "[:alpha:]"))) +
    theme(plot.title = element_text(size=16, hjust=0.5, face="bold"),
          axis.title.y = element_text(color = "lightgoldenrod3", size=16), 
          axis.title.x = element_text(size=16),
          axis.text = element_text(vjust = 0.5, hjust=0, face= "bold", size=12))
  
  
  if (patientname == "P017") {  
    p2<- UC_data.df %>% filter(patient == patientname) %>% 
      ggplot(aes(x=day, y=serum_creatinine)) +
      geom_area(fill="lightgoldenrod3", alpha=0.5) +
      geom_line(color="lightgoldenrod3", size=1) +
      geom_point(color="lightgoldenrod3", size=1) +
      theme_minimal()+
      scale_y_continuous(name = "serum\ncreatinine\n[mg/dl]\n", limits = c(0, 5)) +
      xlim(min(UC_data.df$day), 13) +
      geom_hline(yintercept = 1, linetype = "solid", color = "black", size = 1, alpha = 0.1) +
      geom_segment(aes(x = 3.5, y = 0, xend = 6.5, yend = 0), size=0.5, color = "red")+ 
      geom_text(label= "CRRT", x=5, y=0.4, size=5, color ="red") +
      labs(title = paste0("P0", str_remove_all(patientname, "[:alpha:]"))) +
      theme(plot.title = element_text(size=16, hjust=0.5, face="bold"),
            axis.title.y = element_text(color = "lightgoldenrod3", size=16), 
            axis.title.x = element_text(size=16),
            axis.text = element_text(vjust = 0.5, hjust=0, face= "bold", size=12))
    
  }
  if (patientname == "P018") {  
    p2<- UC_data.df %>% filter(patient == patientname) %>% 
      ggplot(aes(x=day, y=serum_creatinine)) +
      geom_area(fill="lightgoldenrod3", alpha=0.5) +
      geom_line(color="lightgoldenrod3", size=1) +
      geom_point(color="lightgoldenrod3", size=2) +
      theme_minimal()+
      scale_y_continuous(name = "serum\ncreatinine\n[mg/dl]\n", limits = c(0, 5)) +
      xlim(min(UC_data.df$day), 13) +
      geom_hline(yintercept = 1, linetype = "solid", color = "black", size = 1, alpha = 0.1) +
      geom_segment(aes(x = 1.5, y = 0, xend = 3.5, yend = 0), size=0.5, color = "red")+ 
      geom_segment(aes(x = 6.5, y = 0, xend = 12.5, yend = 0), size=0.5, color = "red")+ 
      geom_text(label= "CRRT", x=7.5, y=0.4, size=5, color ="red") +
      labs(title = paste0("P0", str_remove_all(patientname, "[:alpha:]"))) +
      theme(plot.title = element_text(size=16, hjust=0.5, face="bold"),
            axis.title.y = element_text(color = "lightgoldenrod3", size=16), 
            axis.title.x = element_text(size=16),
            axis.text = element_text(vjust = 0.5, hjust=0, face= "bold", size=12))
    
  } 
  
  return(p2)}

## UC_plot_sc creates a plot for single-cells captured each sampling after AKI onset for a single patient based on the data provided in UC_data.df 
UC_plot_sc <- function(patientname){
  
  p3<- UC_data.df %>% filter(patient == patientname) %>% 
    ggplot(aes(x=day, y=single_cell_count)) +
    geom_col(fill ="red", size=3, alpha =0.4) +
    theme_minimal()+
    scale_y_continuous(name = "hq single\ncell trans-\ncriptomes\n", limits = c(0, 6000), labels = scales::label_number_si()) +
    xlim(min(UC_data.df$day), 13) +
    theme(axis.title.y = element_text(color = "red", size=16), 
          axis.title.x = element_text(size=16),
          axis.text = element_text(vjust = 0.5, hjust=0, face= "bold", size=12))
  
  
  return(p3)}

# group_bar() creates a plot similar to Fig 6B, but for a single patient. 
group_bar <- function(SO) { 
  p <- SO@meta.data %>% group_by(AKI_timepoint.fine, renalcellgroup) %>%
    dplyr::summarize(count = n()) %>%
    filter(renalcellgroup != "PDC") %>%
    spread(renalcellgroup, count, fill = 0) %>%
    ungroup() %>%
    mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
    mutate_at(colnames(.[c(2:(ncol(.)-1))]), function(x) x/.$total_cell_count) %>%
    dplyr::select(c('AKI_timepoint.fine', 'total_cell_count', everything())) %>%
    arrange(factor(AKI_timepoint.fine, levels = levels(SO$AKI_timepoint.fine))) %>%
    dplyr::select(-c('total_cell_count')) %>%
    reshape2::melt(id.vars = 'AKI_timepoint.fine') %>%
    ungroup()
  p$AKI_timepoint.fine <- as.numeric(p$AKI_timepoint.fine)
  
  p <-ggplot(p, aes(x=AKI_timepoint.fine, y=value, fill=variable)) +
    geom_area(colour= "black", alpha = 0.6, size=0.25) +
    scale_x_continuous(limits=c(min(UC_data.df$day), 13)) +
    ylab("Renal\ncell type\nproportion [%]\n") +
    xlab("\nTime since AKI onset [d]\n") +
    scale_y_continuous(breaks=number_ticks(1)) +
    scale_fill_manual(values = renalcellgroup_colors, name = "renal cell type") +
    theme_minimal() +
    theme(legend.position = "bottom",   
          axis.title.y = element_text(size=16), 
          axis.title.x = element_text(size=16, hjust = 2),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16), 
          axis.text = element_text(vjust = 0.5, hjust=0, face= "bold", size=12)) +
    geom_segment(x=min(p$AKI_timepoint.fine), y=0, yend=1, xend=min(p$AKI_timepoint.fine), size=0.25) +
    geom_segment(x=max(p$AKI_timepoint.fine), y=0, yend=1, xend=max(p$AKI_timepoint.fine), size=0.25) +
    geom_segment(x=min(p$AKI_timepoint.fine), y=0, yend=0, xend=max(p$AKI_timepoint.fine), size=0.25)
  
  return(p)
}

## load clinical data ----
UC_data.df <- read.delim("~/Urine_AKI_urine_output_crea_data.csv", sep=";")
colnames(UC_data.df)[1] <- "patient"

## create plots, compile and save ----
layout <- "ABCD
EFGH
IJKL
MNOP"

pw <- UC_plot_cr("P023") +  UC_plot_cr("P017") + UC_plot_cr("P024") + UC_plot_cr("P002") +
  UC_plot_uo("P023") +  UC_plot_uo("P017") + UC_plot_uo("P024") + UC_plot_uo("P002") +
  UC_plot_sc("P023") +  UC_plot_sc("P017") + UC_plot_sc("P024") + UC_plot_sc("P002") +
  group_bar(RENAL_P23) + group_bar(RENAL_P17) + group_bar(RENAL_P24) + group_bar(RENAL_P02) +
  plot_layout(design = layout, guides = "collect") & theme(legend.position = 'bottom')

for (i in 1:12) { pw[[i]] <- pw[[i]] + theme(axis.text.x = element_blank(),
                                             axis.ticks.x = element_blank(),
                                             axis.title.x = element_blank() ) }
for (i in c(2:4, 6:8, 10:12, 14:16)) { pw[[i]] <- pw[[i]] + theme(axis.text.y = element_blank(),
                                                                  axis.ticks.y = element_blank(),
                                                                  axis.title.y = element_blank())+ 
  NoLegend() }
for (i in c(1:13,15:16)) { pw[[i]] <- pw[[i]] + theme(axis.title.x = element_blank())}

ggsave(plot = pw , "fig6c.png", w = 4200, h = 2200, units = "px", scale = 1, type = "cairo-png")
