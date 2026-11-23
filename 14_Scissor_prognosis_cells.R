# Authored by LiuHS
# The current species can only be set to Homo sapiens
# install and library pkgs
for (i in c(1:10)) {
	suppressPackageStartupMessages({
		if (!require(devtools)){
			install.packages("devtools")
		}
		if (!require(BiocManager)){
			install.packages("BiocManager")
		}
		if (!require(Scissor)) {
			devtools::install_github('sunduanchen/Scissor')
		}
		if (!require(Seurat)) {
			install.packages("Seurat")
		}
		if (!require(tidyverse)) {
			install.packages("tidyverse")
		}
		if (!require(optparse)) {
			install.packages("optparse")
		}
		if (!require(data.table)) {
			install.packages("data.table")
		}
		if (!require(clusterProfiler)) {
			BiocManager::install("clusterProfiler")
		}
		if (!require(org.Hs.eg.db)) {
			BiocManager::install("org.Hs.eg.db")
		}
		
		library(Scissor)
		library(Seurat)
		library(tidyverse)
		library(clusterProfiler)
		library(optparse)
		library(data.table)
		library(org.Hs.eg.db)
	})
}
# make Parameters Selection
option_list <- list(
	make_option(c("-w", "--workdir"), help = "script work directory ,defualt is run directory",default = "./07_scissor"),
	make_option(c("-s","--rds"),help = "seurat rds object",default = "./03_annotation/sce_anno_sample3000.qs"),
	make_option(c("-r","--res"),help = "resolution for Findcluster",default = 0.15),
	make_option(c("-e","--exp"),help = "load TCGA exp FPKM or counts",default = "./exp.csv"),
	make_option(c("-m","--meta"),help = "metadata for TCGA samples",default = "./metadata.csv"),
	make_option(c("-c","--cutoff",help = "cutoff cells used for scissor,0.03,0.05"),default = 0.05),
	make_option(c("-a","--alpha",help = "Alpha cells used for scissor,0.01,0.05,0.1,0.2"),default = 0.1),
	make_option(c("-n","--name",help = "you condition name"),default = "Condition")
)
opt_parser=OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (!dir.exists(opt$workdir)) {
	dir.create(opt$workdir)
}
result_dir <- opt$workdir

sce <- qs::qread(opt$rds)
# sce <- subset(sce, downsample = 1000)
# sc_dataset <- sce@assays$RNA@counts # need counts matrix of cancer cells
print("Successfully imported single cell matrix")
# run Seurat for readuction
# sc_dataset <- Seurat_preprocessing(sc_dataset, resolution = opt$res,min.cells = 400,min.features = 3)
sc_dataset <- sce
# load TCGA data
print("load TCGA data")
bulk_dataset <- readRDS("../bulkRNAseq/Rdata/exprSet_TARGET_OS_count.rds")
bulk_dataset <- as.matrix(bulk_dataset)
# rownames <- bulk_dataset$V1
# bulk_dataset <- as.matrix(bulk_dataset[,-1])
# rownames(bulk_dataset) <- rownames
# load meta data
metadata <- readRDS("../bulkRNAseq/Rdata/survData_TARGET_OS.rds")
colnames(metadata) = c("samples", "status", "time")
metadata = metadata[,c(1,3,2)]

identical(colnames(bulk_dataset), metadata$samples)

# bulk_dataset <- bulk_dataset[,metadata$samples]

if (identical(metadata$samples,colnames(bulk_dataset))) {
	if (colnames(metadata)[2] == "time") {
		print("Now run survival model")
		phenotype <- metadata[,2:3]
		colnames(phenotype) <- c("time", "status")
		infos1 <- Scissor(bulk_dataset, sc_dataset, phenotype, alpha = opt$alpha, 
											family = "cox", 
											Save_file = paste(result_dir,'Scissor_survival.RData',sep = "/"))
		# infos2 <- Scissor(bulk_dataset, sc_dataset, phenotype, alpha = NULL, cutoff = opt$cutoff, 
		# 									family = "cox", Load_file = paste(result_dir,'Scissor_survival.RData',sep = "/"))
		# infos1 <- Scissor(bulk_dataset, sc_dataset, phenotype, alpha = infos2$para$alpha, 
		# 									family = "cox", cutoff = opt$cutoff,
		# 									Save_file = paste(result_dir,'Scissor_survival.RData',sep = "/"))
		# Reliability Significance Test
		# load(paste(result_dir,'Scissor_survival.RData',sep = "/"))
		# result1 <- reliability.test(X, Y, network, alpha = 0.05, family = "cox", cell_num = numbers, n = 100, nfold = 10)
		# evaluate_summary <- evaluate.cell(paste(result_dir,'Scissor_survival.RData',sep = "/"), infos1, FDR_cutoff = 0.05, bootstrap_n = 100)
		# 
		# # write out results
		# write.csv(evaluate_summary, paste(result_dir,'evaluate_summary.csv',sep = "/"))
		# write(paste("Reliability_test p value =",result1$p),paste(result_dir,'Reliability_test.txt',sep = "/"))
		
		# umap figures
		Scissor_select <- rep("none", ncol(sc_dataset))
		names(Scissor_select) <- colnames(sc_dataset)
		Scissor_select[infos1$Scissor_pos] <- "Scissor+"
		Scissor_select[infos1$Scissor_neg] <- "Scissor-"
		sc_dataset <- AddMetaData(sc_dataset, metadata = Scissor_select, col.name = "scissor")
		sc_dataset$scissor <- factor(sc_dataset$scissor,levels = c("Scissor+","Scissor-","none"))
		DimPlot(sc_dataset, reduction = 'tsne', group.by = 'scissor', cols = c("indianred1","royalblue","lightgray"), pt.size = .1)
		ggsave(paste(result_dir,"tsne_scissor.png",sep = "/"), width = 6, height = 5)
		ggsave(paste(result_dir,"tsne_scissor.pdf",sep = "/"), width = 6, height = 5)
		
		
		Scissor_select <- rep("none", ncol(sc_dataset))
		names(Scissor_select) <- colnames(sc_dataset)
		Scissor_select[infos1$Scissor_pos] <- "Favorable"
		Scissor_select[infos1$Scissor_neg] <- "Unfavorable"
		sc_dataset <- AddMetaData(sc_dataset, metadata = Scissor_select, col.name = "prognosis")
		sc_dataset$prognosis <- factor(sc_dataset$prognosis,levels = c("Favorable","Unfavorable","none"))
		DimPlot(sc_dataset, reduction = 'tsne', group.by = 'prognosis', cols = c("indianred1","royalblue","lightgray"), pt.size = .1)
		ggsave(paste(result_dir,"tsne_scissor_prognosis.png",sep = "/"), width = 6, height = 5)
		ggsave(paste(result_dir,"tsne_scissor_prognosis.pdf",sep = "/"), width = 6, height = 5)
		## scissor cells
		sce_result <- subset(sc_dataset,scissor %in% c("Scissor+","Scissor-"))
		saveRDS(sce_result,paste(result_dir,"sce_result.rds",sep = "/"))
		saveRDS(sc_dataset,paste(result_dir,"sce_result_all.rds",sep = "/"))
		### T test DEGs
		DefaultAssay(sce_result) <- "RNA"
		table(sce_result$scissor)
		Idents(sce_result) <- "scissor"
		DEG <- Seurat::FindMarkers(sce_result,ident.1 = "Scissor+",ident.2 = "Scissor-")
		DEG <- DEG[which(DEG$p_val_adj <= 0.05),]
		write.csv(DEG,paste(result_dir,"DEG_result.csv",sep = "/"))
	} else if (colnames(metadata)[2] == "status") {
		print("Now run  Classification Model")
		phenotype <- metadata$status
		names(phenotype) <- metadata$samples
		#identical(names(phenotype),colnames(bulk_dataset))
		tag <- c(0:(length(unique(metadata$status))-1))
		infos1 <- Scissor(bulk_dataset, sc_dataset, phenotype, tag = tag, alpha = opt$alpha, 
											family = "binomial", 
											Save_file = paste(result_dir,"Scissor_Classify_Model.RData",sep = "/"))
		
		infos2 <- Scissor(bulk_dataset, sc_dataset, phenotype, tag = tag, cutoff = opt$cutoff, 
											family = "binomial", Load_file = paste(result_dir,"Scissor_Classify_Model.RData",sep = "/"))
		infos1 <- Scissor(bulk_dataset, sc_dataset, phenotype, tag = tag, alpha = infos2$para$alpha,cutoff = opt$cutoff, 
											family = "binomial", 
											Save_file = paste(result_dir,"Scissor_Classify_Model.RData",sep = "/"))
		
		# Reliability Significance Test
		load(paste(result_dir,'Scissor_Classify_Model.RData',sep = "/"))
		numbers <- length(infos1$Scissor_pos) + length(infos1$Scissor_neg)
		result1 <- reliability.test(X, Y, network, alpha = 0.05, family = "binomial", cell_num = numbers, n = 100, nfold = 10)
		evaluate_summary <- evaluate.cell(paste(result_dir,'Scissor_Classify_Model.RData',sep = "/"), infos1, FDR_cutoff = 0.05, bootstrap_n = 100)
		
		# write out results
		write.csv(evaluate_summary, paste(result_dir,'evaluate_summary.csv',sep = "/"))
		write(paste("Reliability_test p value =",result1$p),paste(result_dir,'Reliability_test.txt',sep = "/"))
		
		# umap figures
		Scissor_select <- rep("none", ncol(sc_dataset))
		names(Scissor_select) <- colnames(sc_dataset)
		Scissor_select[infos1$Scissor_pos] <- "Scissor+"
		Scissor_select[infos1$Scissor_neg] <- "Scissor-"
		sc_dataset <- AddMetaData(sc_dataset, metadata = Scissor_select, col.name = "scissor")
		sc_dataset$scissor <- factor(sc_dataset$scissor,levels = c("Scissor+","Scissor-","none"))
		DimPlot(sc_dataset, reduction = 'umap', group.by = 'scissor', cols = c("indianred1","royalblue","grey"), pt.size = 1)
		ggsave(paste(result_dir,"umap_scissor.png",sep = "/"))
		ggsave(paste(result_dir,"umap_scissor.pdf",sep = "/"))
		
		
		Scissor_select <- rep("none", ncol(sc_dataset))
		names(Scissor_select) <- colnames(sc_dataset)
		Scissor_select[infos1$Scissor_pos] <- paste(opt$name,"+",sep = "")
		Scissor_select[infos1$Scissor_neg] <- paste(opt$name,"-",sep = "")
		sc_dataset <- AddMetaData(sc_dataset, metadata = Scissor_select, col.name = "Condition")
		
		sc_dataset$Condition <- factor(sc_dataset$Condition,levels = c(paste(opt$name,"+",sep = ""),paste(opt$name,"-",sep = ""),"none"))
		DimPlot(sc_dataset, reduction = 'umap', group.by = 'Condition', cols = c("indianred1","royalblue","grey"), pt.size = 1)
		ggsave(paste(result_dir,"umap_scissor_Condition.png",sep = "/"))
		ggsave(paste(result_dir,"umap_scissor_Condition.pdf",sep = "/"))
		## scissor cells
		sce_result <- subset(sc_dataset,scissor %in% c("Scissor+","Scissor-"))
		saveRDS(sce_result,paste(result_dir,"sce_result.rds",sep = "/"))
		saveRDS(sc_dataset,paste(result_dir,"sce_result_all.rds",sep = "/"))
		### T test DEGs
		DefaultAssay(sce_result) <- "RNA"
		table(sce_result$scissor)
		# Idents(sce_result) <- "scissor"
		Idents(sc_dataset) <- "scissor"
	
		DEG <- Seurat::FindMarkers(sc_dataset,ident.1 = "Scissor+")
		DEG <- DEG[which(DEG$p_val_adj <= 0.05),]
		write.csv(DEG,paste(result_dir,"DEG_result.csv",sep = "/"))
	} else {
		stop("Please check you metadata colnames must be samples time status")
	}
} else {
	stop("Please check you metadata and bulk data,Ensure consistent sample order!")
}


library(gplots)
balloonplot(table(sc_dataset$scissor,sc_dataset$CellType))


export::graph2pdf(file = paste(result_dir,"Scissor_ballon_plot.pdf",sep = "/"), width = 8, height = 6)
export::graph2png(file = paste(result_dir,"Scissor_ballon_plot.png",sep = "/"), width = 8, height = 6)

########################## KEGG and Go function analysis
annotationDB = "org.Hs.eg.db"
orgnm = "hsa"
Genelist <- read.csv(paste(result_dir,"DEG_result.csv",sep = "/"),row.names = 1)
Genelist <- rownames(Genelist)
gene = bitr(Genelist,
						fromType = "SYMBOL",
						toType = "ENTREZID",
						OrgDb = annotationDB)
go <-
	enrichGO(
		gene          = gene$ENTREZID,
		#universe = all_entrezid$ENTREZID,
		OrgDb         = annotationDB,
		ont           = "ALL",
		pAdjustMethod = "BH",
		pvalueCutoff  = 0.05,
		qvalueCutoff  = 0.05,
		readable      = TRUE
	)
go_res = go@result
barplot(go, showCategory = 10,color = "qvalue")
ggsave(paste(result_dir,"GO_pathways_barplot.pdf",sep = "/"),
			 width = 8,
			 height = 6)
ggsave(paste(result_dir,"GO_pathways_barplot.png",sep = "/"),
			 width = 8,
			 height = 6)


dotplot(go, showCategory = 10)


ggsave(paste(result_dir,"GO_pathways_dotplot.pdf",sep = "/"),
			 width = 8,
			 height = 6)
ggsave(paste(result_dir,"GO_pathways_dotplot.png",sep = "/"),
			 width = 8,
			 height = 6)

write.csv(go_res,paste(result_dir,"go_res.csv",sep = "/"))


kk = enrichKEGG(
	gene = gene$ENTREZID,
	organism = orgnm,
	#universe = all_entrezid$ENTREZID,
	pAdjustMethod = "BH",
	pvalueCutoff = 0.05,
	qvalueCutoff = 1
)

reskk = kk@result
reskk = reskk[reskk$p.adjust < 0.05, ]

barplot(kk, showCategory = 10)

ggsave(paste(result_dir,"KEGG_Pathways_barplot.pdf",sep = "/"),
			 width = 8,
			 height = 6)

ggsave(paste(result_dir,"KEGG_Pathways_barplot.png",sep = "/"),
			 width = 8,
			 height = 6)
dotplot(kk, showCategory = 10)
ggsave(paste(result_dir,"KEGG_Pathways_dotplot.pdf",sep = "/"),
			 width = 8,
			 height = 6)

ggsave(paste(result_dir,"KEGG_Pathways_dotplot.png",sep = "/"),
			 width = 8,
			 height = 6)

write.csv(reskk,paste(result_dir,"KEGG_result.csv",sep = "/"))
