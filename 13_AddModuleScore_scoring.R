rm(list = ls())
project_dir = getwd()
genescore_dir=paste(project_dir,"05_Genescore",sep="/")
if(!file.exists(genescore_dir)){
  dir.create(genescore_dir)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggradar)
  library(dplyr)
  library(scales)
  library(tibble)
  #library(CellChat)
  library(ggsci)
})
sce = qs::qread("./03_annotation/sce_anno_sample3000.qs")
# Idents(sce) <- "Celltype"

DefaultAssay(sce) <- "RNA"


os.genes = readRDS("../bulkRNAseq/Rdata/os_realated_genes.rds")

geneset = list(os.genes)
sce <- Seurat::AddModuleScore(object = sce,features = geneset,name = "os.PD")
# colnames(sce@meta.data)[12:ncol(sce@meta.data)] = names(genelist)
colnames(sce@meta.data)

score_dat = sce@meta.data
readr::write_excel_csv(score_dat, file = paste0(genescore_dir, "/PCD_GeneScore_AddModuleScore.csv"))




#------------- Vln + Featureplot
cells_numbers <- length(unique(FetchData(sce,"CellType")[,1]))
celltypes = as.character(unique(FetchData(sce,"CellType")[,1]))
cell_cols = CellChat::scPalette(cells_numbers)
names(cell_cols) = sort(celltypes)

my_comparisons = list()
for (cell in levels(sce@meta.data$CellType)[-3]) {
  my_comparisons = c(list(c(celltypes[1], cell)), my_comparisons)
}

library(ggpubr)
ggviolin(sce@meta.data, x = "CellType", y = "os.PD1",
         color = "CellType",add = 'mean_sd',fill = 'CellType',ylab = "Score",
         add.params = list(color = "black")) + 
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") +
  scale_color_manual(values = cell_cols) + 
  scale_fill_manual(values = cell_cols) +
  theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1)) + 
  NoLegend() + labs(x = '')
ggsave(paste0(genescore_dir, "/","os.PD_Score_AddModuleScore_vln.pdf"),width = 8,height = 6)
ggsave(paste0(genescore_dir, "/","os.PD_Score_AddModuleScore_vln.png"),width = 8,height = 6)


df<- data.frame(sce@meta.data, sce@reductions$tsne@cell.embeddings)
class_avg <- df %>%
  group_by(CellType) %>%
  summarise(
    tSNE_1 = median(tSNE_1),
    tSNE_2 = median(tSNE_2)
  )

ggplot(df, aes(tSNE_1, tSNE_2))  +
  geom_point(aes(colour  = os.PD1)) + scale_colour_gradientn(colors = c("skyblue","white", "#e4002b")) +
  ggrepel::geom_label_repel(aes(label = CellType),
                            data = class_avg,
                            size = 2.5,
                            label.size = 1,
                            segment.color = NA
  )+   theme_bw()+
  labs(title = "OS PD-related GeneSet", color = "Score")+
  theme(plot.title = element_text(hjust = 0.5),panel.grid = element_blank(),legend.position = c(0.9,0.85))
ggsave(paste0(genescore_dir, "/","os.PD_Score_AddModuleScore_tsne.pdf"),width = 6,height = 5.5)
ggsave(paste0(genescore_dir, "/","os.PD_Score_AddModuleScore_tsne.png"),width = 6,height = 5.5)



library(ggpubr)

umap_dir = paste0(genescore_dir,"/GeneScore_umap")
if(!dir.exists(umap_dir)){
  dir.create(umap_dir)
}

vln_dir = paste0(genescore_dir,"/GeneScore_vln")
if(!dir.exists(vln_dir)){
  dir.create(vln_dir)
}


# com_dat = combn(as.character(unique(sort(sce@meta.data$CellType))), 2) %>% as.data.frame()
# my_comparisons = as.list(com_dat)

for (i in 1:length(genelist)) {
  
  ggviolin(sce@meta.data, x = "CellType", y = names(genelist[i]),
           color = "CellType",add = 'mean_sd',fill = 'CellType',
           add.params = list(color = "black")) + 
    # stat_compare_means(comparisons = my_comparisons,label = "p.signif") +
    scale_color_manual(values = cell_cols) + 
    scale_fill_manual(values = cell_cols) +
    theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1)) + 
    NoLegend() + labs(x = '')
  ggsave(paste0(vln_dir, "/", names(genelist[i]),"_Score_AddModuleScore_vln.pdf"),width = 8,height = 4.5)
  ggsave(paste0(vln_dir, "/", names(genelist[i]),"_Score_AddModuleScore_vln.png"),width = 8,height = 4.5)
  FeaturePlot(sce,names(genelist[i]), raster = FALSE)
  ggsave(paste0(umap_dir, "/", names(genelist[i]),"_Score_AddModuleScore_umap.pdf"),width = 5.5,height = 5)
  ggsave(paste0(umap_dir, "/", names(genelist[i]),"_Score_AddModuleScore_umap.png"),width = 5.5,height = 5)
  
}


#heatmap

cells_numbers <- length(unique(FetchData(sce,"CellType")[,1]))
celltypes = as.character(unique(FetchData(sce,"CellType")[,1]))
cell_cols = CellChat::scPalette(cells_numbers)
names(cell_cols) = sort(celltypes)

rt = score_dat[,11:ncol(score_dat)]
rt = aggregate(rt[,2:ncol(rt)], list(CellType = rt$CellType), mean) 
rt = tibble::column_to_rownames(rt, "CellType")
rt = rt %>% scale(center = T) %>% t()

ann_cols = data.frame(CellType = colnames(rt), row.names = colnames(rt))
ann_colors = list(
  CellType = cell_cols
)

library(RColorBrewer)
colors = colorRampPalette(brewer.pal(9, "RdPu"))(25)
# colors = rev(colorRampPalette(brewer.pal(9, "RdBu"))(25))

library(pheatmap)

pheatmap(rt, cluster_cols = F, cluster_rows = T,show_colnames = F,color = colors,
         annotation_col = ann_cols,
         annotation_colors = ann_colors,
         cellwidth = 25,angle_col = 0)
dev.copy2pdf(file = paste0(genescore_dir, "/Heatmap.PCD_GeneScore_AddModuleScore.pdf"), width = 8, height = 6)
