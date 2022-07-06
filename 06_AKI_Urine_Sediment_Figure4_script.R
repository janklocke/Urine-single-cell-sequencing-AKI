# description ----

# This is the R script for generating Figure 4 in 
# "Urinary single-cell sequencing captures intrarenal injury and repair processes 
# in human acute kidney injury" by Klocke et al. 

# in here, we use the proximal tubular cell ischemia reperfusion dataset by Kirita et al. to train a multinomial model 
# using marker genes of the pt cell states and aply the model on the urinary kidney cells to predict their cell state. 
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
library(caret)
library(glmnet)
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

## loading functions ----
# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(x){
  human = useMart(host="https://dec2021.archive.ensembl.org", "ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart(host="https://dec2021.archive.ensembl.org", "ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , 
                   mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  return(genesV2)
}

# functions to generate prediction matrices 
make.conf.matrix <- function(orig,pred){
  ret = matrix(0, ncol = length(unique(pred))+1, nrow = length(unique(orig)))
  colnames(ret) = c(sort(unique(pred)),"error")
  rownames(ret) = sort(unique(orig))
  for(clust in rownames(ret)){
    for(clust2 in colnames(ret)){
      print(c(clust,clust2))
      ret[clust,clust2] = length(orig[orig==clust&pred==clust2])
    }
    ret[clust,] = round(ret[clust,]/sum(ret[clust,]),2)
    ret[clust,"error"] = 1-ret[clust,clust]
  }
  return(ret)
}

make.pred.matrix <- function(orig,pred){
  ret = matrix(0, ncol = length(unique(pred)), nrow = length(unique(orig)))
  colnames(ret) = c(sort(unique(pred)))
  rownames(ret) = sort(unique(orig))
  for(clust in rownames(ret)){
    for(clust2 in colnames(ret)){
      print(c(clust,clust2))
      ret[clust,clust2] = length(orig[orig==clust&pred==clust2])
    }
    ret[clust,] = round(ret[clust,]/sum(ret[clust,]),2)
  }
  return(ret)
}


## loading reference data ----
# the following SO "tmp" and metadata contain all proximal tubular cells from 
# the ischemia reperfusion dataset by Kirita et al. and was kindly provided by the authors.

#SO
tmp= readRDS("/hum.PT.complete.orig.data.rds")
#meta data
meta.hum.pt = read.table("/Yuhei_pt_subclust.csv",
                         sep=",", header = T, row.names = 1)


## rename clusters 
ClusterNames <- c("Inj", "Inj", "healthy", "healthy", "failed-repair", "healthy", "repair", "Sev.Inj")
tmp <- SetIdent(tmp, value ="seurat_clusters")
names(ClusterNames) <- levels(tmp)
tmp <- RenameIdents(tmp, ClusterNames)
tmp$ct <- tmp@active.ident

## Find marker genes in tmp ----
markers.pt.hum = FindAllMarkers(tmp, min.pct = .25, only.pos = T, verbose = T, group.by = "ct")


# Fig 4AB ----

cols <- c(
  paletteer_d("ggsci::light_green_material")[7], 
  paletteer_d("ggsci::red_material")[c(2,5)], 
  paletteer_d("ggsci::light_blue_material")[5], 
  paletteer_d("ggsci::amber_material")[8])  


## Visualization of the mouse PT IRI data as UMAP dimplot ----
b <-  DimPlot(tmp) +   scale_color_manual(values = cols, name = "PT cell states", limits= c("healthy", "Inj", "Sev.Inj", "repair", "failed-repair"), 
                                          labels = c("healthy PT S1-3", "injured PT", "severe injury", "repairing PT", "failed-repair PT")) + ggtitle("Mouse PT ischemia\nreperfusion\ninjury (Kirita et al.)")

## and of key marker genes for healthy, injured and reactive cell states ----
a <-  DotPlot(tmp, features = c("Slc5a12","Adra1a","Bcat1", "Slc22a30","Slc7a13", "Myc","Havcr1","Krt20","Hspa1a", "Top2a",
                                "Vcam1", "Dcdc2a","Sema5a","Ccl2"),
              scale.max=30) + 
  scale_y_discrete(name = "PT cell states", limits= c("healthy", "Inj", "Sev.Inj", "repair", "failed-repair"), 
                   labels = c("healthy PT S1-3", "injured PT", "severe injury", "repairing PT", "failed-repair PT"))+
  coord_flip()+
  scale_color_viridis(direction=-1)+
  theme(axis.text.x = element_text(angle = 30, hjust=1))

## compile and save ----
ggsave(plot = b/a+plot_layout(heights = c(1,3)), "fig4ab.png", scale = 1, w = 1500, h = 3000, units = "px", type = "cairo-png")


# create a list of shared marker genes between tmp (mouse pt data) and RENAL (urine kidney cells) ----
tmp.class = sort(setdiff(unique(Idents(tmp)),"missing"))
genes.tmp = list()
genes.tmp[["hum"]] = c()
unique(markers.pt.hum$cluster)
#!
#a. marker genes or
for(clust in tmp.class){
  print(clust)
  genes.tmp2 = markers.pt.hum[markers.pt.hum$cluster==clust & markers.pt.hum$p_val<.05,]$gene
  #genes.tmp2 = markers.pt.hum[markers.pt.hum$cluster==clust & markers.pt.hum$p_val_adj<.05,]$gene
  print(length(genes.tmp2))
  genes.tmp[["hum"]] = unique(c(genes.tmp[["hum"]], genes.tmp2))
}

#b. variable features 
tmp = FindVariableFeatures(tmp, nfeatures = 2000)
genes.tmp[["hum"]] = tmp@assays$RNA@var.features

#c. variable features expressed in at least 10% of some subcluster
length(genes.tmp[["hum"]])

genes.tmp[["hum.conv"]] = convertMouseGeneList(genes.tmp[["hum"]])
genes.tmp[["hum.conv"]] = genes.tmp[["hum.conv"]][!duplicated(genes.tmp[["hum.conv"]]$MGI.symbol),]
genes.tmp[["hum.conv"]] = genes.tmp[["hum.conv"]][!duplicated(genes.tmp[["hum.conv"]]$HGNC.symbol),]
genes.tmp[["hum.human.conv"]] = genes.tmp[["hum.conv"]]$HGNC.symbol
length(genes.tmp[["hum.human.conv"]])

RENAL <- SetIdent(RENAL, value="renalcelltype")
RENALsub <- subset(RENAL, idents = c("PDC", "PDC_inj"), invert=T)
DefaultAssay(RENALsub) = "RNA"
RENALsub <- NormalizeData(RENALsub)
Idents(RENALsub) = "renalcelltype"

markers.RENAL = FindAllMarkers(RENALsub, min.pct = .25, only.pos = T, verbose = T)

genes.tmp[["urine"]] = c()
urine.clusts = unique(Idents(RENALsub))
for(clust in urine.clusts){
  print(clust)
  genes.tmp2 = markers.RENAL[markers.RENAL$cluster==clust & markers.RENAL$p_val<.05,]$gene
  print(length(genes.tmp2))
  genes.tmp[["urine"]] = unique(c(genes.tmp[["urine"]], genes.tmp2))
}
length(genes.tmp[["urine"]])

#!
genes.tmp[["hum.human.conv"]] = genes.tmp[["hum.conv"]]$HGNC.symbol
length(genes.tmp[["hum.human.conv"]])
genes.tmp[["inters"]] = intersect(genes.tmp[["urine"]],genes.tmp[["hum.human.conv"]])

all(genes.tmp$hum.conv$MGI.symbol %in% rownames(tmp))
all(genes.tmp$hum.conv$HGNC.symbol %in% rownames(RENALsub))


length(genes.tmp[["inters"]])

names(genes.tmp)
head(genes.tmp[["hum.conv"]])

genes.tmp2 = genes.tmp[["hum.conv"]][genes.tmp[["hum.conv"]]$HGNC.symbol %in% genes.tmp[["inters"]],]
length(genes.tmp2$MGI.symbol)
rownames(genes.tmp2) = genes.tmp2$MGI.symbol
head(genes.tmp2)

for(clust in tmp.class){
  print(clust)
  print(length(intersect(genes.tmp2$MGI.symbol,
                         markers.pt.hum[markers.pt.hum$cluster==clust & markers.pt.hum$p_val_adj<.05,]$gene)))
}

# randomly select cells of tmp to be assigned as training (2/3) and test (1/3) data ----
set.seed(0)
cells.tmp = list("train"=c(),"test"=c())
tmp.class
tmp.max=100000000
tmp.min=tmp.max
test.len = 1000000
#cat = "seurat_clusters"
cat = "ct"
for(clust in tmp.class){
  print(clust)
  cells.tmp2 = rownames(tmp@meta.data[tmp@meta.data[[cat]]==clust,])
  print(length(cells.tmp2))
  num = min(.66*length(cells.tmp2),tmp.max)
  cells.tmp3 = sample(cells.tmp2, size=num,replace=F)
  print(length(cells.tmp3))
  if(length(cells.tmp3)<tmp.min){
    tmp.min = length(cells.tmp3)
  }
  cells.tmp[["train"]] = c(cells.tmp[["train"]], cells.tmp3)
  cells.tmp[["test"]] = c(cells.tmp[["test"]], 
                          setdiff(cells.tmp2, 
                                  cells.tmp3)[1:min(test.len,length(setdiff(cells.tmp2, 
                                                                            cells.tmp3)))])
}
tmp.min
length(cells.tmp[["train"]])
length(cells.tmp[["test"]])

## prepare mouse matrix train----
all(tmp@meta.data[cells.tmp[["train"]],]$orig.cellname %in% rownames(meta.hum.pt))
length(rownames(genes.tmp2))
all(genes.tmp2$MGI.symbol %in% rownames(tmp@assays$RNA))
all(genes.tmp2$HGNC.symbol %in% rownames(RENAL@assays$RNA))
all(cells.tmp[["train"]] %in% colnames(tmp@assays$RNA))

x = as.matrix(tmp@assays$RNA[genes.tmp2$MGI.symbol,cells.tmp[["train"]]])
genes.include = apply(x,1,max)>0
genes.include = names(genes.include[genes.include])
genes.tmp2 = genes.tmp2[genes.include,]
x = as.matrix(tmp@assays$RNA[genes.tmp2$MGI.symbol,cells.tmp[["train"]]])
rownames(x) = genes.tmp2[rownames(x),]$HGNC.symbol
length(rownames(x))==length(unique(rownames(x)))
x = scale(t(x))

#y clusters
y = tmp@meta.data[rownames(x),][[cat]]
unique(y)
table(y)
tmp.class
y=factor(y, levels=tmp.class)

## prepare mouse matrix test----
x.test = as.matrix(tmp@assays$RNA[genes.tmp2$MGI.symbol,cells.tmp[["test"]]])
rownames(x.test) = genes.tmp2[rownames(x.test),]$HGNC.symbol
# x.test = t(x.test)
goi = apply(x.test,1,max)==0
goi=goi[goi]
goi
x.test = scale(t(x.test))
x.test[,names(goi)] = 0
y.test = tmp@meta.data[rownames(x.test),][[cat]]
table(y.test)


# fit GLM ----


library(doParallel)
registerDoParallel(2)

fit = glmnet(x, y, family = "multinomial", type.multinomial = "ungrouped", standardize = T,
             parallel=T)#, weights = we.ind)


glm.lambdasearch=data.frame(num.comp=numeric(),acc.train=numeric(),
                            acc.test=numeric(),lam=numeric())

i=1
for(la in fit$lambda){
  print(i)
  coef=coef(fit, s = la, exact = TRUE)
  num.comp=length(as.vector(coef[[1]])[as.vector(coef[[1]])!=0])
  pred.train=predict(fit, newx = x, s = la, type = "class")
  pred.train=factor(x=pred.train,levels = tmp.class)
  pred.train=data.frame(pred=pred.train, cell=rownames(x),orig=y)
  rownames(pred.train) = pred.train$cell
  
  all.orig=table(pred.train$orig)
  abs.length=length(pred.train$orig)
  
  fails=unique(pred.train[pred.train$pred!=pred.train$orig,]$cell)
  fails.orig = table(pred.train[fails,]$orig)
  acc=1-length(fails)/abs.length
  
  #test data
  pred.test=predict(fit, newx = x.test, s = la, type = "class")
  acc.test=1-length(y.test[y.test!=pred.test])/length(y.test)
  fails.test=y.test[y.test!=pred.test]
  fails.test.orig=table(fails.test)
  all.test.orig = table(y.test)
  tmp.frame=data.frame(num.comp=num.comp,acc.train=round(acc,2),
                       acc.test=round(acc.test,2),lam=la)
  glm.lambdasearch=rbind(glm.lambdasearch,tmp.frame)
  i=i+1
}

glm.lambdasearch

## determine lambda with highest accuracy ----
la=glm.lambdasearch[100,]$lam

pred.train=predict(fit, newx = x, s = la, type = "class")
make.conf.matrix(y,pred.train)
pred.test=predict(fit, newx = x.test, s = la, type = "class")
make.conf.matrix(y.test,pred.test)

1-length(y.test[y.test!=pred.test])/length(y.test)
fails.test=y.test[y.test!=pred.test]
fails.test.orig=table(fails.test)
all.test.orig = table(y.test)
fails.test.orig/all.test.orig

# predict human urine kidney cells ----
all(genes.tmp2$HGNC.symbol %in% rownames(RENALsub))
goi=setdiff(genes.tmp2$HGNC.symbol,rownames(RENALsub))
x.test = as.matrix(RENALsub@assays$RNA[genes.tmp2$HGNC.symbol,])
#x.test = scale(t(x.test))
x.test = scale(t(x.test))
y.test = RENALsub@meta.data[rownames(x.test),]$renalcelltype
table(y.test)
length(y.test)
pred.test=predict(fit, newx = x.test, s = la, type = "class")
table(pred.test)
length(pred.test)
tmp.conf.mat = make.pred.matrix(y.test,pred.test[,1])
tmp.conf.mat


## rename clusters in prediction matrix ----
x <- levels(factor(RENAL$renalcelltype))
x <- str_replace(x, "PCT", "PT")
x <- str_replace(x, "DTL", "TEC_VCAM1")

tmp.conf.mat <- tmp.conf.mat[,c(5,3,2,4,1)]
colnames(tmp.conf.mat) <- c("severe injury", "injury", "healthy", "repair", "failed-repair")
rownames(tmp.conf.mat) <- x[x != "PDC" & x != "PDC_inj"]


# Fig 4C, heatmap of predicted cell states per cluster ----
hm <- Heatmap(t(tmp.conf.mat), column_km =5, col = viridis(10, option = "D"), row_order = colnames(tmp.conf.mat), 
              column_title = "urine AKI TEC clusters", row_title = "mouse IRI\ncell states*",     
              row_title_gp = gpar(fontsize = 60, fontface = "bold"),
              column_title_gp = gpar(fontsize = 60, fontface = "bold"),
              row_names_gp = gpar(fontsize = 40, fontface = "bold"),
              column_names_gp = gpar(fontsize = 40, fontface = "bold"), 
              column_dend_height = unit(100, "mm"),
              column_dend_gp = gpar(lwd = 2),
              column_names_rot =45,
              width = unit(63, "cm"), height = unit(15, "cm"),
              heatmap_legend_param = list(labels_gp = gpar(fontsize  = 40), 
                                          title_gp = gpar(fontsize  = 40), 
                                          title = "fraction\nof cluster", 
                                          title_position = "leftcenter-rot", 
                                          legend_height = unit(8, "cm"))
)

## save ----
png("heatmap.png", width = 2500, height = 1100, type="cairo-png")
hm
dev.off()

# Fig 4D, visualize cell state prediction in UMAP dim red. from Fig 2. ----

## extract cell names for each cell state ----
RENALsub$origcellname <- rownames(RENALsub@meta.data)
y.test = RENALsub@meta.data[rownames(x.test),]$origcellname
table(y.test)
length(y.test)
pred.test=predict(fit, newx = x.test, s = la, type = "class")
table(pred.test)
length(pred.test)
tmp.conf.mat = make.pred.matrix(y.test,pred.test[,1])
tmp.conf.mat
tmp.conf.mat <- tmp.conf.mat[rownames(RENALsub@meta.data),]


RENALsub$mousepred_FR <- as.character(tmp.conf.mat[,"failed-repair"])
RENALsub$mousepred_Rep <- as.character(tmp.conf.mat[,"repair"])
RENALsub$mousepred_Hea <- as.character(tmp.conf.mat[,"healthy"])
RENALsub$mousepred_Inj <- as.character(tmp.conf.mat[,"Inj"])
RENALsub$mousepred_Sev <- as.character(tmp.conf.mat[,"Sev.Inj"])

## highlight each cell state in dim plot ----
cols <- c(viridis(10)[10], makeTransparent(viridis(10)[1], 20))
a <- DimPlot(RENALsub, group.by = "mousepred_FR", pt.size=0.5, order =T)  + ggtitle("failed repair")& 
  scale_color_manual(name = "cell state\nprediction", values = cols, labels = c("negative", "positive")) & 
  theme_void() &
  theme(title = element_text(hjust = 0.5))
b <- DimPlot(RENALsub, group.by = "mousepred_Rep", pt.size=0.5, order =T) + ggtitle("repair")& 
  scale_color_manual(name = "cell state\nprediction", values = cols, labels = c("negative", "positive")) & 
  theme_void() &
  theme(title = element_text(hjust = 0.5))
c <- DimPlot(RENALsub, group.by = "mousepred_Hea", pt.size=0.5, order =T) + ggtitle("healthy")& 
  scale_color_manual(name = "cell state\nprediction", values = cols, labels = c("negative", "positive")) & 
  theme_void() &
  theme(title = element_text(hjust = 0.5))
d <- DimPlot(RENALsub, group.by = "mousepred_Inj", pt.size=0.5, order =T) + ggtitle("injury")& 
  scale_color_manual(name = "cell state\nprediction", values = cols, labels = c("negative", "positive")) & 
  theme_void() &
  theme(title = element_text(hjust = 0.5))
e <- DimPlot(RENALsub, group.by = "mousepred_Sev", pt.size=0.5, order =T) + ggtitle("severe injury")& 
  scale_color_manual(name = "cell state\nprediction", values = cols, labels = c("negative", "positive")) & 
  theme_void() &
  theme(title = element_text(hjust = 0.5))

f <-DimPlot(RENALsub, label=T, repel=T, label.size=4.5, pt.size=0.5)+ scale_color_manual(values=makeTransparent(renalcelltype_colors,70))+theme_void()+NoLegend()

## compile and save ----
layout <- "AAABC
AAADE
AAAFG"
pw <- f+guide_area()+c+d+e+a+b+plot_layout(design=layout, guides="collect") &theme(text=element_text(size=16, vjust=0.5))
ggsave(plot = pw, "fig4d.png", w = 2500, h = 2000, units = "px", type = "cairo-png")
