
rm(list = ls())

library(limma)

pd = readRDS("./Rdata/phenoData_GSE87624.rds")

exprSet = readRDS("./Rdata/exprSet_GSE87624_tpm.rds")
# exprSet = log2(exprSet + 1)
group_sample = pd[pd$state %in% c("metastasis", "primary"),]
table(group_sample$state)


exprSet = exprSet[, match(group_sample$sampleID, colnames(exprSet))]

identical(colnames(exprSet), group_sample$sampleID)


project_dir = getwd()
result_dir = paste0(project_dir, "/results/DEG/")
if(!dir.exists(result_dir)){
  dir.create(result_dir)
}


logFC_cutoff = 0.58
pvalue_cutoff = 0.05



exprSet=normalizeBetweenArrays(exprSet)
# boxplot(exprSet,outline=FALSE, notch=T,col=group, las=2)
exprSet <- log2(exprSet+1)



design <- model.matrix(~0+factor(group_sample$state))  
colnames(design) <- c("metastasis", "primary")
row.names(design) <- colnames(exprSet)

cont.matrix <- makeContrasts(contrasts="metastasis-primary", 
                             levels = design)


df.fit <- lmFit(exprSet, design)
fit <- contrasts.fit(df.fit, cont.matrix)
fit <- eBayes(fit)
DEG <- topTable(fit,n = Inf, adjust = "fdr", coef="metastasis-primary")

k1 = (DEG$P.Value < pvalue_cutoff) & (DEG$logFC < -logFC_cutoff)
k2 = (DEG$P.Value < pvalue_cutoff) & (DEG$logFC > logFC_cutoff)
DEG$change = ifelse(k1, "DOWN", ifelse(k2, "UP", "NOT"))
table(DEG$change)


saveRDS(DEG, file =  "./Rdata/DEG_metastasis-primary.rds" )
write.csv(DEG, file = paste0(result_dir, "DEG_metastasis-primary.csv"), quote = F)


library(ggVolcano)

DEG$gene <- row.names(DEG)
data <- add_regulate(DEG, 
                     log2FC_name = "logFC",
                     fdr_name = "P.Value",
                     log2FC = logFC_cutoff, 
                     fdr = pvalue_cutoff)

ggvolcano(data, x = "log2FoldChange", y = "padj",log2FC_cut = logFC_cutoff,
          label = "gene", label_number = 20,add_label = T,output = FALSE,legend_position="DR")+
  theme(panel.grid = element_blank())

ggsave(paste0(result_dir, "vocalo_metastasis-primary.pdf"), width = 8, height = 6)
ggsave(paste0(result_dir, "vocalo_metastasis-primary.png"), width = 8, height = 6)
          
