library(Seurat)
library(dplyr)
rm(list = ls())
project_dir = getwd()

cell_marker_dir=paste(project_dir,"04_Cell_Markers",sep="/")
if(!file.exists(cell_marker_dir)){
  dir.create(cell_marker_dir)
}

sce = qs::qread("./03_annotation/sce_anno.qs")


#find marker
DefaultAssay(sce) <- "RNA"
# Idents(sce) <- "CellType"
sce <- ScaleData(sce, feature=rownames(sce),verbose = FALSE)

markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers,paste(cell_marker_dir,"all_markers.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)

# markers_0.5 <- markers[markers$avg_log2FC > 0.5,]
# write.table(markers_0.5,paste(cell_marker_dir,"all_markers_cutoff_0.5.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)


#cluster summary
marker_number<-table(markers$cluster)%>% reshape2::melt()
colors=hue_pal()(nrow(marker_number))
marker_number$Cluster <- as.factor(marker_number$Var1)
ggplot(data=marker_number,mapping=aes(x=Cluster,y=value,fill=Cluster))+
  geom_bar(stat="identity",width=0.5)+theme_classic()+ylab("Number")+
  geom_text(aes(label = value),size = 3,hjust = 0.5,vjust = -0.5, position = "stack")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave(paste(cell_marker_dir,"marker_number.pdf",sep="/"),width = 9,height = 6)
ggsave(paste(cell_marker_dir,"marker_number.png",sep="/"),width = 9,height = 6)


cells_numbers <- length(unique(FetchData(sce,"CellType")[,1]))
celltypes = as.character(unique(FetchData(sce,"CellType")[,1]))
cols = CellChat::scPalette(cells_numbers)
names(cols) = celltypes
cols = c("#E41A1C", "#BC9DCC", "#4DAF4A", "#F29403", "#984EA3", "#F781BF",  "#377EB8" )
cluster_top10_markers=markers %>% group_by(cluster)%>% top_n(n = 5, wt = avg_log2FC)
markers = c("CD79A","MZB1", # B cell
            "CD3D","CD3E",         # T cell
            "EPCAM", "KRT18","KRT8", # Epi
            "CD68","CD163","CD14", # Macro
            "COL1A1","COL1A2","DCN", # fibroblasts 
            "TPSAB1",                # Mast Cell
            "PECAM1","CLDN5"       # Endo
)

library(scRNAtoolVis)
pdf(paste(cell_marker_dir,"cell_top5_markers_heatmap.pdf",sep="/"),width = 10,height = 8)
AverageHeatmap(object = sce,annoCol = T, myanCol = cols,
               markerGene = cluster_top10_markers$gene,
               row_title = "Cluster top 10 Marker genes",
               # markGenes = markers,
               # showRowNames = F,
               column_names_rot = 30)
dev.off()

png(paste(cell_marker_dir,"cell_top5_markers_heatmap.png",sep="/"),width = 900,height = 800)
AverageHeatmap(object = sce,annoCol = T, myanCol = cols,
               markerGene = cluster_top10_markers$gene,
               row_title = "Cluster top 10 Marker genes",
               # markGenes = markers,
               # showRowNames = F,
               column_names_rot = 30)
dev.off()

cluster_top5_markers=markers %>% group_by(cluster)%>% top_n(n = 5, wt = avg_log2FC)
library(scRNAtoolVis)
pdf(paste(cell_marker_dir,"cell_top5_markers_heatmap.pdf",sep="/"),width = 10,height = 8)
AverageHeatmap(object = sce,
               markerGene = cluster_top5_markers$gene,
               row_title = "Cluster top 5 Marker genes",
               column_names_rot = 30)
dev.off()

vol.markers <- FindAllMarkers(sce, only.pos = FALSE,
                              min.pct = 0.25,
                              logfc.threshold = 0)

cells_numbers <- length(unique(FetchData(sce,"CellType")[,1]))
celltypes = as.character(unique(FetchData(sce,"CellType")[,1]))
cols = CellChat::scPalette(cells_numbers)
names(cols) = sort(celltypes)

jjVolcano(diffData = vol.markers,
          tile.col = cols,
          size  = 3.5,
          fontface = 'italic',
          polar = T) +
  ylim(-10,12)
ggsave(paste(cell_marker_dir,"cell_markers_manhattan.pdf",sep="/"), width = 11, height = 10)
ggsave(paste(cell_marker_dir,"cell_markers_manhattan.png",sep="/"), width = 11, height = 10)


library(scRNAtoolVis)
jjVolcano(diffData = vol.markers,
          log2FC.cutoff = 0.7,
          col.type = "adjustP",
          topGeneN = 5,tile.col =cols,legend.position = c(0.92,0.95))+
  theme(axis.text.x = element_text(size = 8,face = "bold")) 
ggsave(paste(cell_marker_dir,"cell_markers_manhattan.pdf",sep="/"), width = 12, height = 7)
ggsave(paste(cell_marker_dir,"cell_markers_manhattan.png",sep="/"), width = 12, height = 7)




for( clust_num in  unique(Idents(sce))){
  cluster_dir=paste(cell_marker_dir,paste("cluster",clust_num,sep="_"),sep="/")
  if(!file.exists(cluster_dir)){
    dir.create(cluster_dir)
  }
  # cluster_dir_enrich=paste(cluster_dir,"enrichment",sep="/")
  # if(!file.exists(cluster_dir_enrich)){dir.create(cluster_dir_enrich)}
  cluster_markers=subset(markers,cluster==clust_num)
  rownames(cluster_markers)<-cluster_markers$gene
  if(nrow(cluster_markers)>1){
    genelist=cluster_markers$gene
    # enrichment(species=opt$type,outDir=cluster_dir_enrich,geneList=genelist)
    write.table(cluster_markers,paste(cluster_dir,paste("cluster",clust_num,"markers.xls",sep="_"),sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
    top10_markers=cluster_markers %>%  top_n(n = 10, wt = avg_log2FC)
    VlnPlot(sce, features = top10_markers$gene,pt.size = 0.1 ,ncol=5)
    
    ggsave(paste(cluster_dir,"top10_vilion.pdf",sep="/"),width =20,height = 7)
    ggsave(paste(cluster_dir,"top10_vilion.png",sep="/"),width =20,height = 7)
    
    FeaturePlot(sce, features = top10_markers$gene, min.cutoff = "q9",ncol=5)
    ggsave(paste(cluster_dir,"top10_cell_exp_distribution.pdf",sep="/"),width = 20,height = 7)
    ggsave(paste(cluster_dir,"top10_cell_exp_distribution.png",sep="/"),width = 20,height = 7)
    if(length(sample_list)>1){
      DotPlot(sce, features = rev(unique(top10_markers$gene)), cols = rainbow(length(sample_list))) + RotatedAxis()}else{
        DotPlot(sce, features = rev(unique(top10_markers$gene)))}
    
    ggsave(paste(cluster_dir,"top10_exp_pct.pdf",sep="/"),width = 15,height =15)
    ggsave(paste(cluster_dir,"top10_exp_pct.png",sep="/"),width = 15,height =15)
    
  }}
