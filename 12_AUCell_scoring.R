library(AUCell)
library(Seurat)
library(SeuratData)
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggplot2)
rm(list = ls())

project_dir = getwd()
aucell_dir=paste(project_dir,"05_AUCell",sep="/")
if(!file.exists(aucell_dir)){
  dir.create(aucell_dir)
}


sce = qs::qread("./03_annotation/sce_anno_sample3000.qs")

# DimPlot(sce, reduction = "umap", label = TRUE, pt.size = 0.5,group.by = "CellType") + NoLegend()

# Idents(seu_imm) <- "cell_type"

exprMatrix <- sce@assays$RNA@data
cells_rankings <- AUCell_buildRankings(exprMatrix, splitByBlocks = T)

export::graph2pdf(file = paste0(aucell_dir, "/Gene_detected.pdf"), width = 7,height = 5)
export::graph2png(file = paste0(aucell_dir, "/Gene_detected.png"), width = 7,height = 5)


os.genes = readRDS("../bulkRNAseq/Rdata/os_realated_genes.rds")


cells_AUC_PD <- AUCell_calcAUC(os.genes, cells_rankings, nCores = 10,aucMaxRank=nrow(cells_rankings)*0.05)


set.seed(12345)
cells_assignment_PCD <- AUCell_exploreThresholds(cells_AUC_PD, plotHist=TRUE, assign=TRUE)

export::graph2pdf(file = paste0(aucell_dir, "/PD_AUC_histogram.pdf"), width = 7,height = 6)
export::graph2png(file = paste0(aucell_dir, "/PD_AUC_histogram.png"), width = 7,height = 6)


aucMat_PCD <- getAUC(cells_AUC_PD)

sce@meta.data$AUC_PCD = as.numeric(aucMat_PCD)
sce@meta.data$AUC_group = ifelse(sce@meta.data$AUC_PCD > 0.094, "High", "Low")

# saveRDS(sce,file = paste(aucell_dir, "sce.AUCell.rds", sep = "/") )
qs::qsave(sce, file = paste0(aucell_dir, "/sce.AUCell.qs"))

meta.data = sce@meta.data
write_excel_csv(meta.data,file = paste(aucell_dir, "PD_AUCscore.csv", sep = "/"))


df<- data.frame(sce@meta.data, sce@reductions$tsne@cell.embeddings)
class_avg <- df %>%
  group_by(CellType) %>%
  summarise(
    tSNE_1 = median(tSNE_1),
    tSNE_2 = median(tSNE_2)
    )
colorPal <- grDevices::colorRampPalette(c("gray", "#87CEFA"))(5)

ggplot(df, aes(tSNE_1, tSNE_2))  +
  geom_point(aes(colour  = AUC_PCD)) + scale_colour_gradientn(colors = c("skyblue","white", "red")) +
  ggrepel::geom_label_repel(aes(label = CellType),
                            data = class_avg,
                            size = 2.5,
                            label.size = 1,
                            segment.color = NA
                            )+   theme_bw()+
  labs(title = "PD Geneset AUC Score", color = "AUC")+
  theme(plot.title = element_text(hjust = 0.5),panel.grid = element_blank(),legend.position = c(0.9,0.85))
ggsave(paste(aucell_dir, "tSNE.PD.AUC.Score.pdf", sep = "/"),width = 7,height = 6.5)
ggsave(paste(aucell_dir, "tSNE.PD.AUC.Score.png", sep = "/"),width = 7,height = 6.5)


ggplot(data.frame(sce@meta.data, sce@reductions$tsne@cell.embeddings), aes(tSNE_1, tSNE_2, color=AUC_group)
) + geom_jitter( size=0.5
) +  scale_color_manual(values = c("#E35050","lightgrey"))  + theme_light(base_size = 15)+labs(title = "PD Geneset AUC Group", color = "")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),panel.grid = element_blank(),legend.position = c(0.9,0.85))
ggsave(paste(aucell_dir, "tSNE.PD.AUC.Group.pdf", sep = "/"),width = 7,height = 6.5)
ggsave(paste(aucell_dir, "tSNE.PD.AUC.Group.png", sep = "/"),width = 7,height = 6.5)



cells_numbers <- length(unique(FetchData(sce,"CellType")[,1]))
celltypes = as.character(unique(FetchData(sce,"CellType")[,1]))
cell_cols = CellChat::scPalette(cells_numbers)
names(cell_cols) = sort(celltypes)

my_comparisons = list()
for (cell in levels(sce@meta.data$CellType)[-3]) {
  my_comparisons = c(list(c(celltypes[1], cell)), my_comparisons)
}

library(ggpubr)
ggviolin(sce@meta.data, x = "CellType", y = "AUC_PCD",
         color = "CellType",add = 'mean_sd',fill = 'CellType',ylab = "Score",
         add.params = list(color = "black")) + 
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") +
  scale_color_manual(values = cell_cols) + 
  scale_fill_manual(values = cell_cols) +
  theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1)) + 
  NoLegend() + labs(x = '')
ggsave(paste0(aucell_dir, "/","os.PD_AUC.Score_vln.pdf"),width = 8,height = 6)
ggsave(paste0(aucell_dir, "/","os.PD_AUC.Score_vln.png"),width = 8,height = 6)

