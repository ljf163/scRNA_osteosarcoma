# R version 4.3.2
suppressPackageStartupMessages({
  library(infercnv)
  library(rjags)
  library(Seurat) 
  library(optparse)
  library(tidyverse)
})

option_list <- list(
  make_option(c("-r", "--rds"), help = " ", default = "./sce_anno_sample3000.qs"),
  make_option(c("-o","--out"),help = "out dir",default = "05_infercnv"),
  make_option(c("-c","--core"),help = "Core",default = 10),
  make_option(c("-v","--reverse"),help = "Reverse the biological process ",default = "NULL")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

################################################################################
sce <- qs::qread(opt$rds)
table(sce$CellType)
sce <- subset(sce,CellType %in% c("Chondroblastic cells", "Osteoblastic cells", "Myeloid cells"))
count <- as.matrix(sce@assays$RNA@counts)

metadata <- sce@meta.data
head(metadata)
metadata$barcode <- rownames(metadata)


meta <- data.frame(barcode = rownames(metadata),celltype = metadata$CellType)
meta <- na.omit(meta)
write.table(meta,"meta.txt",quote = T,row.names = F,col.names = F)

meta_new <- read.table("./meta.txt",
                       header = FALSE, row.names = 1,stringsAsFactors = FALSE,
                       colClasses = c("character", "character"))

ref <- read.table("gencode_v21_gen_pos.complete.txt")
head(ref)
ref <- dplyr::distinct(ref,V1,.keep_all = T)
rownames(ref) <- ref$V1
ref <- ref[,-1]
head(ref)



infercnv_obj = CreateInfercnvObject(raw_counts_matrix = count,
                                    annotations_file = meta_new,
                                    gene_order_file = ref,
                                    ref_group_names = c("Myeloid cells"),
                                    chr_exclude=c("chrY", "chrM"))


infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff= 0.1,
                             out_dir="05_infercnv",
                             cluster_by_groups=TRUE,
                             denoise = TRUE,
                             HMM = FALSE,
                             num_threads = 20,
                             write_expr_matrix = T,
                             min_cells_per_gene = 10,
                             no_prelim_plot = T
)
################################################################################

library(RColorBrewer)

infercnv_obj <- readRDS("./05_infercnv/run.final.infercnv_obj")
infercnv::plot_cnv(infercnv_obj,
                   output_filename = "inferCNV_heatmap",
                   output_format = "pdf",
                   custom_color_pal =  color.palette(c("#2067AE","white","#B11F2B")))


library(ComplexHeatmap)
library(colorRamps)
expr <- infercnv_obj@expr.data

gene_pos <- read.delim("./gencode_v21_gen_pos.complete.txt", header = F)
gene_pos <- gene_pos[gene_pos$V1 %in% rownames(expr),]

new_cluster <- unique(gene_pos$V2)

top_color <- HeatmapAnnotation(cluster = anno_block(labels = gsub("chr", "", new_cluster),
                                                    gp = gpar(col = "white"),
                                                    labels_gp = gpar(cex = 1, col = "black"),
                                                    height = unit(5,"mm")))


ref_cell <- c(infercnv_obj@reference_grouped_cell_indices$`Myeloid cells`)

obs_cell <- c(infercnv_obj@observation_grouped_cell_indices$`Chondroblastic cells`, 
              infercnv_obj@observation_grouped_cell_indices$`Osteoblastic cells`)

cell_anno =data.frame(cell_id =c(colnames(expr)[ref_cell],colnames(expr)[obs_cell]),
                      group =c(rep("Myeloid",length(ref_cell)),
                               rep("Chondroblastic",length(infercnv_obj@observation_grouped_cell_indices$`Chondroblastic cells`)),
                               rep("Osteoblastic",length(infercnv_obj@observation_grouped_cell_indices$`Osteoblastic cells`))
                               ))



set.seed(123)
kmeans.result <- kmeans(t(expr), 5)
kmeans_df <- data.frame(kmeans.result$cluster)
colnames(kmeans_df) <- "k_cluster"
kmeans_df <- as_tibble(cbind(cell_id = rownames(kmeans_df), kmeans_df))
kmeans_df=kmeans_df%>%inner_join(cell_anno,by="cell_id")%>%arrange(k_cluster)
kmeans_df$k_cluster=as.factor(kmeans_df$k_cluster)




#method3------------------------------------------------------------------------

obs <- read.table("./05_infercnv/infercnv.observations.txt", header=T, check.names=F)
ref <- read.table("./05_infercnv/infercnv.references.txt", header=T, check.names=F)

estimateCNV <- function(obs,
                        ref, 
                        score_threshold, 
                        cor_threshold
){
  cnv_obs <- colMeans((obs-1)^2)
  cnv_ref <- colMeans((ref-1)^2)
  cell_top <- names(sort(cnv_obs, decreasing=T))[1:round(length(cnv_obs)*0.5)]
  cnv_top <- rowMeans(obs[, cell_top])
  cor_obs <- apply(obs, 2, function(x)cor(x, cnv_top))
  cor_ref <- apply(ref, 2, function(x)cor(x, cnv_top))

  cnv <- data.frame(score=c(cnv_obs, cnv_ref),
                    cor=c(cor_obs, cor_ref),
                    barcode=c(colnames(obs), colnames(ref)))


  cnv$type <- 'Other'
  cnv$type[cnv$score>score_threshold & cnv$cor>cor_threshold] <- 'Malignant'
  cnv$type[cnv$score<score_threshold & cnv$cor<cor_threshold] <- 'Not Malignant'
  return(cnv)
}

CNV_score_cor <- estimateCNV(obs, ref,
                             score_threshold= 0.005,
                             cor_threshold= 0.4)
colnames(CNV_score_cor)[3] <- "cell_id"
CNV_score_cor=CNV_score_cor%>%inner_join(kmeans_df,by="cell_id")

CNV_score_cor_Epi <- subset(CNV_score_cor, group %in% c("Chondroblastic", "Osteoblastic"))

library(ggplot2)
ggplot(CNV_score_cor_Epi, aes(x=score,y=cor,color=type))+
  geom_point()+
  theme_bw()+
  theme(axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title = element_text(color = 'black', size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(x="CNV Score",
       y="CNV Correlation")+
  scale_color_manual(values =c("#009E73","grey","black") )
table(CNV_score_cor_Epi$k_cluster,CNV_score_cor_Epi$type)
ggsave("CNV_Score.pdf", width = 7, height = 5)
ggsave("CNV_Score.png", width = 7, height = 5)

dir.create("./Malignant_identified")
saveRDS(CNV_score_cor_Epi,"./CNV_score_cor.rds")
write.csv(CNV_score_cor_Epi, file = "./CNV_score_cor.csv", row.names = F, quote = F)
################################################################################
rm(list = ls())
