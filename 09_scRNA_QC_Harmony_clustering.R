
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
options(stringsAsFactors = FALSE)


#------------------ QC
project_dir = getwd()
### GSE152048
sce = qs::qread("./00_Rawdata/sce_raw_GSE152048.qs")

cell_num_raw=table(sce$sample)
##QC ,filter low Q cell 
seurat_qc_dir=paste(project_dir,"01_QC_Filter",sep="/")

if(!file.exists(seurat_qc_dir)){
  dir.create(seurat_qc_dir)
}
#save basic merge rds 
# saveRDS(sce,paste0(seurat_qc_dir,"/Basicinfo_inter.rds"))

#get mt percentage 
sce[["percent.MT"]] <- PercentageFeatureSet(sce, pattern = "^MT-")


##output experiment raw information
out_gene_umi=as.matrix(sce@meta.data)
row_name=c("Cell",rownames(out_gene_umi))
out_gene_umi=rbind(colnames(out_gene_umi),out_gene_umi)
out_gene_umi=as.data.frame(cbind(as.matrix(row_name),out_gene_umi))
write.table(out_gene_umi,file=paste(seurat_qc_dir,"GSE152048_Basicinfo_nGene_nUMI_mito.txt",sep="/"),sep="\t",quote=F,row.names=F,col.names=F)

#filter 
cellqc.raw=VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA","percent.MT"), ncol = 3,pt.size = 0,group.by="sample")
sce <- subset(sce, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.MT <10)
cellqc.filter=VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA","percent.MT"), ncol = 3,pt.size = 0,group.by="sample")


#plot_grid(cellqc.raw,cellqc.filter,nrow=2)
cellqc.raw/cellqc.filter
ggsave(paste(seurat_qc_dir,"GSE152048_cells_qc_filter.pdf",sep="/"),width = 17,height = 8)
ggsave(paste(seurat_qc_dir,"GSE152048_cells_qc_filter.png",sep="/"),width = 17,height = 8)

cell_num_filter=table(sce$sample)

filter_summary=cbind(cell_num_raw,cell_num_filter)
filter_summary=data.frame("cell_num_raw"=filter_summary[,1],"cell_num_filter"=filter_summary[,2])
filter_summary[nrow(filter_summary)+1,]=apply(filter_summary,2,sum)
rownames(filter_summary)[nrow(filter_summary)]="total"
filter_summary$percent=apply(filter_summary,1,function(x){round(x[2]*100/x[1],2)})

write.table(filter_summary,paste(seurat_qc_dir,"GSE152048_cells_filter_stat.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
qs::qsave(sce, file = paste(seurat_qc_dir, "sce_qc_GSE152048.qs",sep="/"))

####GSE162454
sce = qs::qread("./00_Rawdata/sce_raw_GSE162454.qs")
cell_num_raw=table(sce$sample)
##QC ,filter low Q cell 
seurat_qc_dir=paste(project_dir,"01_QC_Filter",sep="/")

if(!file.exists(seurat_qc_dir)){
  dir.create(seurat_qc_dir)
}
#save basic merge rds 
# saveRDS(sce,paste0(seurat_qc_dir,"/Basicinfo_inter.rds"))

#get mt percentage 
sce[["percent.MT"]] <- PercentageFeatureSet(sce, pattern = "^MT-")


##output experiment raw information
out_gene_umi=as.matrix(sce@meta.data)
row_name=c("Cell",rownames(out_gene_umi))
out_gene_umi=rbind(colnames(out_gene_umi),out_gene_umi)
out_gene_umi=as.data.frame(cbind(as.matrix(row_name),out_gene_umi))
write.table(out_gene_umi,file=paste(seurat_qc_dir,"GSE162454_Basicinfo_nGene_nUMI_mito.txt",sep="/"),sep="\t",quote=F,row.names=F,col.names=F)

#filter 
cellqc.raw=VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA","percent.MT"), ncol = 3,pt.size = 0,group.by="sample")
sce <- subset(sce, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.MT <10)
cellqc.filter=VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA","percent.MT"), ncol = 3,pt.size = 0,group.by="sample")


#plot_grid(cellqc.raw,cellqc.filter,nrow=2)
cellqc.raw/cellqc.filter
ggsave(paste(seurat_qc_dir,"GSE162454_cells_qc_filter.pdf",sep="/"),width = 17,height = 8)
ggsave(paste(seurat_qc_dir,"GSE162454_cells_qc_filter.png",sep="/"),width = 17,height = 8)

cell_num_filter=table(sce$sample)

filter_summary=cbind(cell_num_raw,cell_num_filter)
filter_summary=data.frame("cell_num_raw"=filter_summary[,1],"cell_num_filter"=filter_summary[,2])
filter_summary[nrow(filter_summary)+1,]=apply(filter_summary,2,sum)
rownames(filter_summary)[nrow(filter_summary)]="total"
filter_summary$percent=apply(filter_summary,1,function(x){round(x[2]*100/x[1],2)})

write.table(filter_summary,paste(seurat_qc_dir,"GSE162454_cells_filter_stat.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
qs::qsave(sce, file = paste(seurat_qc_dir, "sce_qc_GSE162454.qs",sep="/"))








#----------------- remove batch ---harmoney

rm(list = ls())
gc()


project_dir = getwd()
seurat_exp_cluster_dir=paste(project_dir,"02_Cluster",sep="/")
if(!file.exists(seurat_exp_cluster_dir)){dir.create(seurat_exp_cluster_dir)}

GSE152048 = qs::qread("./01_QC_Filter/sce_qc_GSE152048.qs")
GSE162454 = qs::qread("./01_QC_Filter/sce_qc_GSE162454.qs")
sce = merge(GSE152048, GSE162454)

qc_summary <- sce@meta.data %>%
  group_by(sample) %>%
  summarise(
    total_cells = n(),
    mean_genes = mean(nFeature_RNA) %>% round(),
    mean_UMI = mean(nCount_RNA)%>% round(),
    mean_percent.MT = mean(percent.MT)%>% round(digits = 2),
    median_genes = median(nFeature_RNA)%>% round(),
    median_UMI = median(nCount_RNA)%>% round(),
    median_percent.MT = median(percent.MT)%>% round(digits = 2)
  )

openxlsx::write.xlsx(qc_summary, file = paste(seurat_exp_cluster_dir, "qc_summary.xlsx", sep = "/"))

sce <- sce %>%
  Seurat::NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData()
sce <- RunPCA(sce, npcs = 50, verbose = FALSE)

sce.pca <- sce

library(harmony)
sce=sce %>% RunHarmony(c("sample", "Batch"), plot_convergence = TRUE)


sce.int <-sce
#ElbowPlot
ElbowPlot(sce.int,ndims = 30)
# ggsave(paste(seurat_exp_cluster_dir,"pca.ElbowPlot.pdf",sep="/"),width = 7,height = 6)
# ggsave(paste(seurat_exp_cluster_dir,"pca.ElbowPlot.png",sep="/"),width = 7,height = 6)

sce.int <- RunUMAP(sce.int, reduction = "harmony", dims = 1:30)
sce.int <- RunTSNE(sce.int, reduction = "harmony", dims = 1:30)
sce.int <- FindNeighbors(sce.int, reduction = "harmony", dims = 1:30)
sce.int = FindClusters(sce.int, resolution = 0.1)

library(ggplot2)
sce.int <- RunUMAP(sce.int, reduction = "pca", dims = 1:30,reduction.name = "umap.unintegrated")
sce.int <- RunTSNE(sce.int, reduction = "pca", dims = 1:30,reduction.name = "tsne.unintegrated")
library(patchwork)
DimPlot(sce.int, reduction = "tsne.unintegrated", group.by = "Batch", raster = FALSE,pt.size = .1) + NoAxes() +
  labs(title = "Batch effect correction before") +
(DimPlot(sce.int, reduction = "tsne", group.by = "Batch", raster = FALSE, pt.size = .1) + NoAxes()+
  labs(title = "Batch effect correction after")) +plot_layout(guides = "collect")

ggsave(paste(seurat_exp_cluster_dir,"Batch effect correction.pdf",sep="/"),width = 13,height = 6)
ggsave(paste(seurat_exp_cluster_dir,"Batch effect correction.png",sep="/"),width = 13,height = 6)



VlnPlot(sce.int, features = c("nFeature_RNA", "nCount_RNA","percent.MT"), ncol = 3,pt.size = 0,group.by="seurat_clusters")

# res_list = list()
# for(i in seq(0.2,1, 0.1)){
#   sce = FindClusters(sce.int, resolution = i)
#   res_list = c(res_list,sce)
# }
# saveRDS(res_list, file = paste(seurat_exp_cluster_dir,"sce_res_0.2to1.rds",sep="/"))


# resolution=0.8

DimPlot(sce.int, reduction = "umap", group.by = "Batch", raster = FALSE) + NoAxes()
ggsave(paste(seurat_exp_cluster_dir,"UMAP_Batch_harmony.pdf",sep="/"),width = 6.5,height = 6)
ggsave(paste(seurat_exp_cluster_dir,"UMAP_Batch_harmony.png",sep="/"),width = 6.5,height = 6)

p1 <- DimPlot(sce.int, reduction = "umap", group.by = "sample", raster = FALSE)
p2 <- DimPlot(sce.int, reduction = "umap", group.by = "seurat_clusters",label = TRUE, raster = FALSE)
plot_grid(p1, p2)
ggsave(paste(seurat_exp_cluster_dir,"UMAP_sample_cluster.pdf",sep="/"),width = 14,height = 7)
ggsave(paste(seurat_exp_cluster_dir,"UMAP_sample_cluster.png",sep="/"),width = 14,height = 7)

p2
ggsave(paste(seurat_exp_cluster_dir,"UMAP_cluster.pdf",sep="/"),width = 6.5,height = 6)
ggsave(paste(seurat_exp_cluster_dir,"UMAP_cluster.png",sep="/"),width = 6.5,height = 6)
p1
ggsave(paste(seurat_exp_cluster_dir,"UMAP_sample.pdf",sep="/"),width = 6.5,height = 6)
ggsave(paste(seurat_exp_cluster_dir,"UMAP_sample.png",sep="/"),width = 6.5,height = 6)


#find marker
# DefaultAssay(sce.int) <- "RNA"
# sce.int <- ScaleData(sce.int, feature=rownames(sce.int),verbose = FALSE)

markers <- FindAllMarkers(sce.int, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers,paste(seurat_exp_cluster_dir,"all_markers.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)


qs::qsave(sce.int, file = paste(seurat_exp_cluster_dir,"sce_res0.1.qs",sep="/"))


sce <- qs::qread("./02_Cluster/sce_res0.1.qs")
DimPlot(sce, label = T)
sce2 <- qs::qread("./02_Cluster/sce_res0.15.qs")
DimPlot(sce2, label = T)
library(cowplot)
library(Seurat)
library(ggplot2)
library(clustree)
resolution_values <- seq(0.1, 0.9, by = 0.1)

sce <- FindClusters(sce, resolution = 0.15)

plots <- list()
for (resolution in resolution_values) {
  sce <- FindClusters(sce, resolution = resolution)
  plot <- DimPlot(sce, reduction = "umap")
  plot <- plot + ggtitle(paste("Resolution =", resolution))
  plots[[length(plots) + 1]] <- plot
}
# final_plot <- plot_grid(plotlist = plots, ncol = 3)
# final_plot
# ggsave(plot = final_plot,filename = paste(dir,"umap_resolution.pdf",sep = "/"),width = 14,height = 14)
# 

# pdf(file = paste(dir,"clustree_plot.pdf",sep = "/"),height = 18,width = 10,onefile = F)
colnames(sce@meta.data)[17] <- "RNA_snn_res.0.1"
colnames(sce@meta.data)[7] <- "res.0.1"
p =clustree(sce@meta.data,prefix = "RNA_snn_res.",layout = c("sugiyama"),node_label_size = 2,node_size_range = c(8, 15))
print(p)

# dev.off()
export::graph2pdf(file = paste(seurat_exp_cluster_dir,"clustree_plot.pdf",sep = "/"),height = 10,width = 15)
export::graph2png(file = paste(seurat_exp_cluster_dir,"clustree_plot.png",sep = "/"),height = 10,width = 15)

