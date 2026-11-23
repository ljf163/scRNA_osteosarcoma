rm(list = ls())

library(Seurat)

project_dir = getwd()
scenic_dir=paste(project_dir,"09_scenic",sep="/")
if(!file.exists(scenic_dir)){
  dir.create(scenic_dir)
}


sce = qs::qread("./03_annotation/sce_anno_sample3000.qs")
sce = subset(sce, downsample = 1000)

qs::qsave(sce, file = "./09_scenic/sce_sample1000.qs")

write.csv(t(as.matrix(sce@assays$RNA@counts)),file = paste(scenic_dir, "for.scenic.data.csv", sep = "/"))





#-----------------------
rm(list=ls())
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(scRNAseq)
library(patchwork)
library(ggplot2) 
library(stringr)
library(circlize)


project_dir = getwd()
scenic_dir=paste(project_dir,"09_scenic",sep="/")
if(!file.exists(scenic_dir)){
  dir.create(scenic_dir)
}



loom <- open_loom('./09_scenic/os.scenic/out_SCENIC.loom') 

regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons_incidMat[1:4,1:4] 
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings <- get_embeddings(loom)  
close_loom(loom)

rownames(regulonAUC)
names(regulons)

#
seurat.data = qs::qread("./09_scenic/sce_sample1000.qs")



sub_regulonAUC <- regulonAUC[,match(colnames(seurat.data),colnames(regulonAUC))]
dim(sub_regulonAUC)
seurat.data
identical(colnames(sub_regulonAUC), colnames(seurat.data))

cellClusters <- data.frame(row.names = colnames(seurat.data), 
                           seurat_clusters = as.character(seurat.data$seurat_clusters))
cellTypes <- data.frame(row.names = colnames(seurat.data), 
                        celltype = seurat.data$CellType)
head(cellTypes)
head(cellClusters)
sub_regulonAUC[1:4,1:4] 

save(sub_regulonAUC,cellTypes,cellClusters,seurat.data,
     file = paste(scenic_dir, 'for_rss_and_visual.Rdata', sep = "/"))


regulonsToPlot = c('VDR(+)','HOXC6(+)')
regulonsToPlot %in% row.names(sub_regulonAUC)
seurat.data@meta.data = cbind(seurat.data@meta.data ,t(assay(sub_regulonAUC[regulonsToPlot,])))


# Vis
p1 = DotPlot(seurat.data, features = unique(regulonsToPlot)) + RotatedAxis()
p2 = RidgePlot(seurat.data, features = regulonsToPlot , ncol = 2) 
p3 = VlnPlot(seurat.data, features = regulonsToPlot,pt.size = 0)
p4 = FeaturePlot(seurat.data,features = regulonsToPlot,reduction = "tsne")

wrap_plots(p1,p2,p3,p4)



# Split the cells by cluster:
selectedResolution <- "celltype" # select resolution
cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes[,selectedResolution])

sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 
dim(sub_regulonAUC)

# Calculate average expression:
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))

# Scale expression. 
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 
dim(regulonActivity_byGroup_Scaled)
#[1] 209   9
regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)


Heatmap(
  regulonActivity_byGroup_Scaled,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = TRUE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)


rss <- calcRSS(AUC=getAUC(sub_regulonAUC), 
               cellAnnotation=cellTypes[colnames(sub_regulonAUC), selectedResolution]) 
rss=na.omit(rss) 
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)



rss=regulonActivity_byGroup_Scaled
head(rss)
df = do.call(rbind,
             lapply(1:ncol(rss), function(i){
               dat= data.frame(
                 path  = rownames(rss),
                 cluster =   colnames(rss)[i],
                 sd.1 = rss[,i],
                 sd.2 = apply(rss[,-i], 1, median)  
               )
             }))
df$fc = df$sd.1 - df$sd.2
top5 <- df %>% group_by(cluster) %>% top_n(5, fc)
rowcn = data.frame(path = top5$cluster) 
n = rss[top5$path,] 
#rownames(rowcn) = rownames(n)
pheatmap(n,
         annotation_row = rowcn,
         show_rownames = T)


sce_regulons <- read.csv("./09_scenic/os.scenic/reg.csv")
sce_regulons <- sce_regulons[-2, ]
colnames(sce_regulons) <- sce_regulons[1,]
sce_regulons <- sce_regulons[-1, ]
colnames(sce_regulons) <- c("TF","ID","AUC","NES","MotifSimilarityQvalue","OrthologousIdentity",
                            "Annotation","Context","TargetGenes","RankAtMax")

sig_regulons <- sce_regulons[which(sce_regulons$AUC>0.1),]
sig_regulons <- sig_regulons[, c("TF","TargetGenes")]
sig_regulons$TargetGenes <-gsub("\\[","",sig_regulons$TargetGenes)
sig_regulons$TargetGenes <-gsub("\\]","",sig_regulons$TargetGenes)
sig_regulons$TargetGenes <-gsub("\\(","",sig_regulons$TargetGenes)
sig_regulons$TargetGenes <-gsub("\\)","",sig_regulons$TargetGenes)
sig_regulons$TargetGenes <-gsub("\\'","",sig_regulons$TargetGenes)

split_sig_regulons<-str_split(sig_regulons$TargetGenes,",")


TF_target = do.call(rbind,
                    lapply(1:length(split_sig_regulons), function(i){
                      tmp = as.data.frame(split_sig_regulons[[i]])
                      targets <- tmp[seq(1,nrow(tmp),2), ]
                      scores <- tmp[seq(0,nrow(tmp),2), ]
                      dat <- data.frame(TF = sig_regulons$TF[i],Target = targets, Score = scores)
                    }))

TF_target <- distinct(TF_target)


pdrgs = readRDS("../bulkRNAseq/Rdata/os_realated_genes.rds")


hub_tfs = TF_target[TF_target$Target %in% pdrgs,]

node_info = data.frame(Type = c(rep("TF", nrow(hub_tfs)), rep("Target", nrow(hub_tfs))),
                       Node = c(hub_tfs$TF, hub_tfs$Target))

node_info = distinct(node_info)

write.csv(hub_tfs, file = paste(scenic_dir, "hub_TFs.csv", sep = "/"), quote = F, row.names = F)
write.csv(node_info, file = paste(scenic_dir, "TF_Target_node_info.csv", sep = "/"), quote = F, row.names = F)


rss=regulonActivity_byGroup_Scaled
head(rss)
df = do.call(rbind,
             lapply(1:ncol(rss), function(i){
               dat= data.frame(
                 path  = rownames(rss),
                 cluster =   colnames(rss)[i],
                 sd.1 = rss[,i],
                 sd.2 = apply(rss[,-i], 1, median)  
               )
             }))
df$fc = df$sd.1 - df$sd.2
df$TF = gsub("\\(.*", "", df$path)

hubs = df[df$TF %in% hub_tfs$TF,]
hubs = subset(hubs, cluster %in% c("Fibroblasts", "Chondroblastic cells", "Osteoblastic cells"))

# top10 <- hubs %>% group_by(cluster) %>% top_n(10, fc)
# rowcn = data.frame(path = top10$cluster) 
n = rss[unique(hubs$path),] 
pheatmap(n,
         # annotation_row = rowcn,
         show_rownames = T)
export::graph2pdf(file = paste(scenic_dir, "TF_heatmap.pdf", sep = "/"), width = 9, height = 9)
export::graph2png(file = paste(scenic_dir, "TF_heatmap.png", sep = "/"), width = 9, height = 9)
