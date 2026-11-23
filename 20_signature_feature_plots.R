rm(list = ls())
library(Seurat)
library(ggplot2)
library(scCustomize)
library(MySeuratWrappers)

project_dir = getwd()
feature_dir=paste(project_dir,"10_signature",sep="/")
if(!file.exists(feature_dir)){
  dir.create(feature_dir)
}


sce = qs::qread("./03_annotation/sce_anno.qs")

markers = read.csv("./04_Cell_Markers/all_markers.xls", sep = "\t")
os.genes = readRDS("../bulkRNAseq/Rdata/os_realated_genes.rds")
markers = markers[markers$gene %in% os.genes,]
markers$difference = markers$pct.1 - markers$pct.2
markers = markers[markers$cluster %in% unique(sce$CellType)[c(1,3,4,5,6,7,9)],]

top5 = markers %>% group_by(cluster)%>% top_n(n = 10, wt = avg_log2FC) %>% top_n(n = 5, wt = difference)
top5 = top5[top5$difference > 0.49,]
top5 = top5$gene
top5[top5 == "HLA-DQA1"] = "LY96"


source("./feature_plot_density.R")
plot_density(sce,marker = top5, dim = "TSNE", rdc = 2, size = .5, ncol = 4)
export::graph2pdf(file = paste(feature_dir, "OS.PD.related.genes.Featureplot.pdf",sep = "/"), width = 12, height = 9)
export::graph2png(file = paste(feature_dir, "OS.PD.related.genes.Featureplot.png",sep = "/"), width = 12, height = 9)

plot_density(sce,marker = "HTRA1", dim = "TSNE", rdc = 2, size = .5)
export::graph2pdf(file = paste(feature_dir, "HTRA1.Featureplot.pdf",sep = "/"), width = 5.5, height = 5)
export::graph2png(file = paste(feature_dir, "HTRA1.Featureplot.png",sep = "/"), width = 5.5, height = 5)



features = readRDS("../bulkRNAseq/Rdata/Features_selected_LASSO.rds")

plot_density(sce,marker = features, dim = "TSNE", rdc = 2, size = .5, ncol = 5)
export::graph2pdf(file = paste(feature_dir, "Model.genes.Featureplot.pdf",sep = "/"), width = 14, height = 9)
export::graph2png(file = paste(feature_dir, "Model.genes.Featureplot.png",sep = "/"), width = 14, height = 9)



plots=list()
for (i in 1:length(features)){
  plots[[i]]=FeaturePlot_scCustom(seurat_object = sce, 
                                  colors_use = colorRampPalette(c("gray", "blue" ))(50), 
                                  features = features[i])+NoAxes() + NoLegend()#+
  #theme(panel.border = element_rect(fill = NA,color = "black",
  #size=1.5,linetype = "solid"))
}

library(patchwork)
p<-wrap_plots(plots, ncol = 5);p


library(ComplexHeatmap)
library(circlize)

col_fun = colorRamp2(c(0,2,4), colorRampPalette(c("gray", "blue" ))(50)[c(1,25,50)])
lgd3 = Legend(at = c(0,2,4), col_fun = col_fun, 
              nrow = 1,#border = "black",direction = "horizontal",
              legend_height = unit(4, "cm")
)
draw(lgd3,x = unit(1, "npc"), y = unit(0.5, "npc"),
     just = c("right", "bottom"))
export::graph2pdf(file = paste(feature_dir,"FeaturePlot_markers.pdf",sep="/"),width = 20,height = 10)
export::graph2png(file = paste(feature_dir,"FeaturePlot_markers.png",sep="/"),width = 20,height = 10)



FeaturePlot(sce, features = features,ncol = 5, raster = F, cols = c("grey", "blue"), reduction = "tsne")
ggsave(paste(feature_dir,"Features_umap.pdf",sep="/"),width = 20,height = 10)
ggsave(paste(feature_dir,"Features_umap.png",sep="/"),width = 20,height = 10)



#--------------------
cells_numbers <- length(unique(FetchData(sce,"CellType")[,1]))
celltypes = as.character(unique(FetchData(sce,"CellType")[,1]))
cols = CellChat::scPalette(cells_numbers)
names(cols) = sort(celltypes)

VlnPlot(sce, features = features,stacked = T,pt.size = 0,cols = cols,direction = "horizontal", x.lab = "", y.lab = '') + 
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
ggsave(paste(feature_dir,"VnPlot_markers.pdf",sep="/"),width = 10,height = 7)
ggsave(paste(feature_dir,"VnPlot_markers.png",sep="/"),width = 10,height = 7)


feature_dir=paste(project_dir,"10_signature/features",sep="/")
if(!file.exists(feature_dir)){
  dir.create(feature_dir)
}

for (i in 1:length(features)) {
  VlnPlot(sce, features = features[i], col = cols)+
    theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1)) + 
    NoLegend() + labs(x = '')
  ggsave(paste0(feature_dir, "/", features[i],"_vln.pdf"),width = 8,height = 4.5)
  ggsave(paste0(feature_dir, "/", features[i],"_vln.png"),width = 8,height = 4.5)
  FeaturePlot(sce,features[i])
  ggsave(paste0(feature_dir, "/", features[i],"_umap.pdf"),width = 5.5,height = 5)
  ggsave(paste0(feature_dir, "/", features[i],"_umap.png"),width = 5.5,height = 5)
  
}




#----------------
geneset = list(features)

sce <- Seurat::AddModuleScore(object = sce,features = geneset,name = "SignatureScore")
colnames(sce@meta.data)[10] = "SignatureScore"
colnames(sce@meta.data)

score_dat = sce@meta.data
write_excel_csv(score_dat, file = paste0(feature_dir, "/SignatureScore_AddModuleScore.csv"))


library(ggpubr)

cell_cols = CellChat::scPalette(14)
names(cell_cols) = sort(unique(sce@meta.data$CellType))

ggviolin(sce@meta.data, x = "CellType", y = "SignatureScore",
         color = "CellType",add = 'mean_sd',fill = 'CellType',
         add.params = list(color = "black")) + 
  # stat_compare_means(comparisons = my_comparisons,label = "p.signif") +
  scale_color_manual(values = cell_cols) + 
  scale_fill_manual(values = cell_cols) +
  theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1)) + 
  NoLegend() + labs(x = '')
ggsave(paste0(feature_dir, "/SignatureScore_AddModuleScore_vln.pdf"),width = 8,height = 4.5)
ggsave(paste0(feature_dir, "/SignatureScore_AddModuleScore_vln.png"),width = 8,height = 4.5)



FeaturePlot(sce, "SignatureScore",reduction = "tsne",raster = F, cols = c("skyblue", "red"))
ggsave(paste0(feature_dir, "/SignatureScore__AddModuleScore_tsne.pdf"),width = 5.5,height = 5)
ggsave(paste0(feature_dir, "/SignatureScore__AddModuleScore_tsne.png"),width = 5.5,height = 5)
