# Cell Annotation

rm(list = ls())
gc()

library(Seurat)
library(dplyr)
library(Seurat)
library(cowplot)
library(ggplot2)
library(reshape2)
# library(DoubletFinder)
library(dittoSeq)
library(scales)
library(MySeuratWrappers)
library(tidydr)
library(stringr)
library(viridis)
library(scCustomize)


sce = qs::qread("./02_Cluster/sce_res0.15.qs")

project_dir = getwd()
anno_dir=paste(project_dir,"03_annotation",sep="/")

if(!file.exists(anno_dir)){
  dir.create(anno_dir)
}


#------------------------ singleR

library(celldex)
# hpca.se = HumanPrimaryCellAtlasData()
# # saveRDS(hpca.se, file = "d:/ProjectR/celldex/HumanPrimaryCellAtlasData.rds")
# # table(hpca.se@colData@listData$label.main)
# # table(hpca.se@colData@listData$label.fine)
# 
# # rm(list = ls())
hpca.se = readRDS("d:/ProjectR/celldex/HumanPrimaryCellAtlasData.rds")



library(SingleR)
seu_obj4SingleR <- GetAssayData(sce, slot = "data")  
cluster_seu_obj = sce@meta.data$seurat_clusters
table(cluster_seu_obj)

pred.hpca <- SingleR(test = seu_obj4SingleR, ref = hpca.se,
                     labels = hpca.se$label.main,clusters = cluster_seu_obj)  
table(pred.hpca$labels)

pdf(paste(anno_dir,"singleR.scoreHeatmap.pdf",sep="/"), width = 12, height = 8)
plotScoreHeatmap(pred.hpca)
dev.off()

celltype.singler = data.frame(clusters = rownames(pred.hpca), celltype = pred.hpca$labels)
colnames(celltype.singler)= c("seurat_clusters","CellType.SingleR")

write.table(celltype.singler,paste(anno_dir,"singleR.cell.anno.xls",sep="/"),sep="\t",quote=F,row.names=F)



#------------------- Manual
DefaultAssay(sce) = "RNA"


markers = c("CD74", "CD14", "FCGR3A",
            "COL1A1", "RUNX2", "IBSP",
            "LUM", "DCN", "FGF7", "COL3A1",
            "IL7R", "CD3D", "NKG7",
            "TOP2A", "PCNA", "MKI67", 
            "CTSK", "MMP9", "ACP5",
            "ACAN", "COL2A1", "SOX9", 
            "PECAM1", "VWF",
            "ACTA2", "RGS5", 
            "MS4A1", "CD79A",
            "MYLPF", "MYL1"
            # "CXCL12", "SFRP2", "MME",
            )

# source("./scripts/Dotplot_anno.R")
# celltypes = c("B cells", "T cells", "Epithelial cells", "Macrophages", "Fibroblasts", "Mast cells", "Endothelial cells")
# cellCols = CellChat::scPalette(7)
# names(cellCols) = sort(celltypes)
# Dotplot_anno(sce, features = markers, celltype_color = cellCols,level = c(0:25),
#              group = c(rep('B cells',3), rep('T cells',2),rep('Epithelial cells',3),
#                        rep('Macrophages',3), rep('Fibroblasts',3),rep("Mast cells",1), rep('Endothelial cells',2)),
#              color = colorRampPalette(c("gray", "red"))(50),
#              order = T)
DimPlot(sce, reduction = "tsne", group.by = "seurat_clusters",label = TRUE, raster = FALSE)

DotPlot(sce, features = markers, group.by = "CellType",cols = c("gray", "red"))+ RotatedAxis()
ggsave(paste(anno_dir,"Dotpplot_markers.pdf",sep="/"),width = 10,height = 7)
ggsave(paste(anno_dir,"Dotpplot_markers.png",sep="/"),width = 10,height = 7)

markers = c("LYZ", "CD68", 
            "CD2", "CD3D", "CD3E", "GNLY", "NKG7", "KLRB1", 
            "COL1A1", "ACTA2", "VIM", 
            "IGHG1", "MZB1", 
            "ALPL", "RUNX2", "IBSP",
            "ACP5", "CTSK", 
            "EGFL7", "PLVAP", 
            "MS4A1", "CD79A"
)




genes = c("COL1A1",
          "LUM",
          "RUNX2",
          "SOX9",
          "CTSK",
          "CD3D",
          "CD74",
          "CD79A",
          "PECAM1",
          "MYLPF",
          "RGS5"
)

plots=list()
for (i in 1:length(genes)){
  plots[[i]]=FeaturePlot_scCustom(seurat_object = sce, reduction = "tsne",
                                  colors_use = colorRampPalette(c("gray", "red" ))(50), 
                                  features = genes[i],raster = F)+NoAxes() + NoLegend()#+
  #theme(panel.border = element_rect(fill = NA,color = "black",
  #size=1.5,linetype = "solid"))
}

library(patchwork)
p<-wrap_plots(plots, ncol = 4);p


library(ComplexHeatmap)
library(circlize)

col_fun = colorRamp2(c(0,2,4), c("gray", "#DF5D5D", "#FF0000"))
lgd3 = Legend(at = c(0,2,4), col_fun = col_fun, 
              nrow = 1,#border = "black",direction = "horizontal",
              legend_height = unit(4, "cm")
)
draw(lgd3,x = unit(0.80, "npc"), y = unit(0.2, "npc"),
     just = c("right", "bottom"))

export::graph2pdf(file = paste(anno_dir,"FeaturePlot_markers.pdf",sep="/"),width = 18,height = 13)
export::graph2png(file = paste(anno_dir,"FeaturePlot_markers.png",sep="/"),width = 18,height = 13,dpi = 300)

marker_dir = paste0(anno_dir, "/featureplot_markers")
if(!file.exists(marker_dir)){
  dir.create(marker_dir)
}

for(gene in markers){
  FeaturePlot(sce, features = gene,reduction = "tsne", cols = c("gray", "red" ), raster = FALSE)
  ggsave(paste0(marker_dir, "/",gene, "_Featureplot.pdf"), width = 6,height = 5)
  ggsave(paste0(marker_dir, "/",gene, "_Featureplot.png"), width = 6,height = 5)
}


#---------- 

anno = readxl::read_xlsx("./03_annotation/Cell_Annotation_Manual.xlsx")

new.cluster.ids <- anno$CellType
Idents(sce) <- sce@meta.data$seurat_clusters
names(new.cluster.ids) <- levels(sce)
sce <- RenameIdents(sce,new.cluster.ids)
DimPlot(sce, reduction = "tsne", label = TRUE, pt.size = 1)
sce$CellType = Idents(sce)

# saveRDS(sce, file = paste0(anno_dir, "/sce_anno.rds"))
qs::qsave(sce, file = paste0(anno_dir, "/sce_anno.qs"))


sce_sub = subset(sce, downsample = 3000)

qs::qsave(sce_sub, file = paste0(anno_dir, "/sce_anno_sample3000.qs"))

#----------------violin markers
cells_numbers <- length(unique(FetchData(sce,"CellType")[,1]))
celltypes = as.character(unique(FetchData(sce,"CellType")[,1]))
cols = CellChat::scPalette(cells_numbers)
names(cols) = sort(celltypes)

VlnPlot(sce, features = markers,stacked = T,pt.size = 0,cols = cols,direction = "horizontal", x.lab = "", y.lab = '') + 
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
ggsave(paste(anno_dir,"VnPlot_markers.pdf",sep="/"),width = 10,height = 7)
ggsave(paste(anno_dir,"VnPlot_markers.png",sep="/"),width = 10,height = 7)

#---------- umap
library(patchwork)
library(Seurat)
library(ggplot2)
library(data.table)
library(tidyverse)
library(ggsci)
library(patchwork)
library(reshape2)



# all
pc12 <- Embeddings(object = sce,reduction = 'umap') %>%
  data.frame()

# check
head(pc12,3)
lower <- floor(min(min(pc12$UMAP_1),min(pc12$UMAP_2)))-1
# get relative line length
linelen <- abs(0.3*lower) + lower
# mid point
mid <- abs(0.3*lower)/2 + lower
# axies data
axes <- data.frame(x = c(lower,lower,lower,linelen),y = c(lower,linelen,lower,lower),
                   group = c(1,1,2,2),
                   label = rep(c('UMAP_2','UMAP_1'),each = 2))
# axies label
label <- data.frame(lab = c('UMAP_2','UMAP_1'),angle = c(90,0),
                    x = c(lower - 0.5,mid),y = c(mid,lower -0.5))

p <- DimPlot(sce,reduction = 'umap',label = T) +
  NoAxes() + NoLegend() +
  theme(aspect.ratio = 1) +scale_color_manual(values = cols)+
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab),fontface  = 'italic') + 
  ggtitle("All")+ theme(plot.title = element_text(hjust = 0.5))
p


sce_normal = subset(sce, Tissue == "Normal")

pc12 <- Embeddings(object = sce_normal,reduction = 'umap') %>%
  data.frame()

# check
head(pc12,3)
lower <- floor(min(min(pc12$UMAP_1),min(pc12$UMAP_2)))-1
# get relative line length
linelen <- abs(0.3*lower) + lower
# mid point
mid <- abs(0.3*lower)/2 + lower
# axies data
axes <- data.frame(x = c(lower,lower,lower,linelen),y = c(lower,linelen,lower,lower),
                   group = c(1,1,2,2),
                   label = rep(c('UMAP_2','UMAP_1'),each = 2))
# axies label
label <- data.frame(lab = c('UMAP_2','UMAP_1'),angle = c(90,0),
                    x = c(lower - 0.5,mid),y = c(mid,lower -0.5))

p1 <- DimPlot(sce_normal,reduction = 'umap',label = T) +
  NoAxes() + NoLegend() +
  theme(aspect.ratio = 1) +scale_color_manual(values = cols)+
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab),fontface  = 'italic') + 
  ggtitle("Normal")+ theme(plot.title = element_text(hjust = 0.5))
p1


sce_tumor = subset(sce, Tissue == "Tumor")

pc12 <- Embeddings(object = sce_tumor,reduction = 'umap') %>%
  data.frame()

# check
head(pc12,3)
lower <- floor(min(min(pc12$UMAP_1),min(pc12$UMAP_2)))-1
# get relative line length
linelen <- abs(0.3*lower) + lower
# mid point
mid <- abs(0.3*lower)/2 + lower
# axies data
axes <- data.frame(x = c(lower,lower,lower,linelen),y = c(lower,linelen,lower,lower),
                   group = c(1,1,2,2),
                   label = rep(c('UMAP_2','UMAP_1'),each = 2))
# axies label
label <- data.frame(lab = c('UMAP_2','UMAP_1'),angle = c(90,0),
                    x = c(lower - 0.5,mid),y = c(mid,lower -0.5))

p2 <- DimPlot(sce_tumor,reduction = 'umap',label = T) +
  NoAxes() + NoLegend() +
  theme(aspect.ratio = 1) +scale_color_manual(values = cols)+
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab),fontface  = 'italic') + 
  ggtitle("Tumor")+ theme(plot.title = element_text(hjust = 0.5))
p2

p + p1 + p2
ggsave(paste0(anno_dir, "/umap_anno_group.pdf"),width = 18, height = 6)
ggsave(paste0(anno_dir, "/umap_anno_group.png"),width = 18, height = 6)

p
ggsave(paste0(anno_dir, "/umap_anno_All.pdf"),width = 7, height = 7)
ggsave(paste0(anno_dir, "/umap_anno_All.png"),width = 7, height = 7)

p1
ggsave(paste0(anno_dir, "/umap_anno_Normal.pdf"),width = 7, height = 7)
ggsave(paste0(anno_dir, "/umap_anno_Normal.png"),width = 7, height = 7)

p2
ggsave(paste0(anno_dir, "/umap_anno_Tumor.pdf"),width = 7, height = 7)
ggsave(paste0(anno_dir, "/umap_anno_Tumor.png"),width = 7, height = 7)

#---------- tsne
if(F){
  
  # all
  pc12 <- Embeddings(object = sce,reduction = 'tsne') %>%
    data.frame()
  
  # check
  head(pc12,3)
  lower <- floor(min(min(pc12$tSNE_1),min(pc12$tSNE_2)))-1
  # get relative line length
  linelen <- abs(0.3*lower) + lower
  # mid point
  mid <- abs(0.3*lower)/2 + lower
  # axies data
  axes <- data.frame(x = c(lower,lower,lower,linelen),y = c(lower,linelen,lower,lower),
                     group = c(1,1,2,2),
                     label = rep(c('tSNE_2','tSNE_1'),each = 2))
  # axies label
  label <- data.frame(lab = c('tSNE_2','tSNE_1'),angle = c(90,0),
                      x = c(lower - 0.5,mid),y = c(mid,lower -0.5))
  
  p <- DimPlot(sce,reduction = 'tsne',label = T) +
    NoAxes() + NoLegend() +
    theme(aspect.ratio = 1) +scale_color_manual(values = cols)+
    geom_line(data = axes,
              aes(x = x,y = y,group = group),
              arrow = arrow(length = unit(0.1, "inches"),
                            ends="last", type="closed")) +
    geom_text(data = label,
              aes(x = x,y = y,angle = angle,label = lab),fontface  = 'italic') + 
    ggtitle("All")+ theme(plot.title = element_text(hjust = 0.5))
  p
  
  
  sce_normal = subset(sce, Tissue == "Normal")
  
  pc12 <- Embeddings(object = sce_normal,reduction = 'tsne') %>%
    data.frame()
  
  # check
  head(pc12,3)
  lower <- floor(min(min(pc12$tSNE_1),min(pc12$tSNE_2)))-1
  # get relative line length
  linelen <- abs(0.3*lower) + lower
  # mid point
  mid <- abs(0.3*lower)/2 + lower
  # axies data
  axes <- data.frame(x = c(lower,lower,lower,linelen),y = c(lower,linelen,lower,lower),
                     group = c(1,1,2,2),
                     label = rep(c('tSNE_2','tSNE_1'),each = 2))
  # axies label
  label <- data.frame(lab = c('tSNE_2','tSNE_1'),angle = c(90,0),
                      x = c(lower - 0.5,mid),y = c(mid,lower -0.5))
  
  p1 <- DimPlot(sce_normal,reduction = 'tsne',label = T) +
    NoAxes() + NoLegend() +
    theme(aspect.ratio = 1) +scale_color_manual(values = cols)+
    geom_line(data = axes,
              aes(x = x,y = y,group = group),
              arrow = arrow(length = unit(0.1, "inches"),
                            ends="last", type="closed")) +
    geom_text(data = label,
              aes(x = x,y = y,angle = angle,label = lab),fontface  = 'italic') + 
    ggtitle("Normal")+ theme(plot.title = element_text(hjust = 0.5))
  p1
  
  
  sce_tumor = subset(sce, Tissue == "Tumor")
  
  pc12 <- Embeddings(object = sce_tumor,reduction = 'tsne') %>%
    data.frame()
  
  # check
  head(pc12,3)
  lower <- floor(min(min(pc12$tSNE_1),min(pc12$tSNE_2)))-1
  # get relative line length
  linelen <- abs(0.3*lower) + lower
  # mid point
  mid <- abs(0.3*lower)/2 + lower
  # axies data
  axes <- data.frame(x = c(lower,lower,lower,linelen),y = c(lower,linelen,lower,lower),
                     group = c(1,1,2,2),
                     label = rep(c('tSNE_2','tSNE_1'),each = 2))
  # axies label
  label <- data.frame(lab = c('tSNE_2','tSNE_1'),angle = c(90,0),
                      x = c(lower - 0.5,mid),y = c(mid,lower -0.5))
  
  p2 <- DimPlot(sce_tumor,reduction = 'tsne',label = T) +
    NoAxes() + NoLegend() +
    theme(aspect.ratio = 1) +scale_color_manual(values = cols)+
    geom_line(data = axes,
              aes(x = x,y = y,group = group),
              arrow = arrow(length = unit(0.1, "inches"),
                            ends="last", type="closed")) +
    geom_text(data = label,
              aes(x = x,y = y,angle = angle,label = lab),fontface  = 'italic') + 
    ggtitle("Tumor")+ theme(plot.title = element_text(hjust = 0.5))
  p2
  
  p + p1 + p2
  ggsave(paste0(anno_dir, "/tsne_anno_group.pdf"),width = 18, height = 6)
  ggsave(paste0(anno_dir, "/tsne_anno_group.png"),width = 18, height = 6)
  
  p
  ggsave(paste0(anno_dir, "/tsne_anno_All.pdf"),width = 7, height = 7)
  ggsave(paste0(anno_dir, "/tsne_anno_All.png"),width = 7, height = 7)
  
  p1
  ggsave(paste0(anno_dir, "/tsne_anno_Normal.pdf"),width = 7, height = 7)
  ggsave(paste0(anno_dir, "/tsne_anno_Normal.png"),width = 7, height = 7)
  
  p2
  ggsave(paste0(anno_dir, "/tsne_anno_Tumor.pdf"),width = 7, height = 7)
  ggsave(paste0(anno_dir, "/tsne_anno_Tumor.png"),width = 7, height = 7)}


#------------ Cell Percentage

cell_cols = cols


cell_counts <- FetchData(sce, vars = c ("CellType", "sample")) 
theme_set(theme_cowplot())

ggplot(data = cell_counts, aes(x = Tissue, fill = CellType)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = cell_cols) +
  coord_flip() + ylab("Percentage")+
  scale_y_reverse()
ggsave(paste0(anno_dir, "/CellPercentage_by_group.pdf"),width = 9, height = 4)
ggsave(paste0(anno_dir, "/CellPercentage_by_group.png"),width = 9, height = 4)

ggplot(data = cell_counts, aes(x = sample, fill = CellType)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = cell_cols) +
  coord_flip() + ylab("Percentage")+
  scale_y_reverse()
ggsave(paste0(anno_dir, "/CellPercentage_by_PatientID.pdf"),width = 9, height = 6)
ggsave(paste0(anno_dir, "/CellPercentage_by_PatientID.png"),width = 9, height = 6)

ggplot(data = cell_counts, aes(x = CellType, fill = Tissue)) +
  geom_bar(position = "fill") +
  ggsci::scale_fill_nejm() +  
  # scale_fill_manual(values = cell_cols) +
  coord_flip() + ylab("Percentage")+
  scale_y_reverse()
ggsave(paste0(anno_dir, "/CellPercentage_by_CellType_group.pdf"),width = 9, height = 7)
ggsave(paste0(anno_dir, "/CellPercentage_by_CellType_group.png"),width = 9, height = 7)

# pcols = CellChat::scPalette(19)[8:19]
# 
# ggplot(data = cell_counts, aes(x = CellType, fill = PatientID)) +
#   geom_bar(position = "fill") +
#   scale_fill_manual(values = pcols) +
#   coord_flip() + ylab("Percentage")+
#   scale_y_reverse()
# ggsave(paste0(anno_dir, "/CellPercentage_by_CellType_PatientID.pdf"),width = 9, height = 7)
# ggsave(paste0(anno_dir, "/CellPercentage_by_CellType_PatientID.png"),width = 9, height = 7)


