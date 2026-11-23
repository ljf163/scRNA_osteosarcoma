rm(list = ls())
library(sva)
library(tidyverse)

exp_tpm = readRDS("./Rdata/exprSet_TARGET_OS_tpm.rds")
tcga=as.data.frame(log2(exp_tpm+1))
tcga = limma::normalizeBetweenArrays(tcga)
boxplot(tcga)


GSE21257 = readRDS("./Rdata/exprSet_GSE21257.rds")
GSE16091 = readRDS("./Rdata/exprSet_GSE16091.rds")


cohort =data.frame(sample = c(colnames(tcga),colnames(GSE16091),colnames(GSE21257)),
                   projectid = c(rep("tcga", ncol(tcga)),
                                 rep("GSE16091", ncol(GSE16091)),
                                 rep("GSE21257", ncol(GSE21257))
                                 ))


tcga = data.frame(genes = rownames(tcga),  tcga)
GSE16091 = data.frame(genes = rownames(GSE16091),  GSE16091)
GSE21257 = data.frame(genes = rownames(GSE21257),  GSE21257)


exp_list = list(tcga, GSE16091, GSE21257)


my_merge <- function(df1, df2){                                # Create own merging function
  merge(df1, df2, by = "genes")
}

exp = Reduce(my_merge, exp_list) 

rownames(exp) = exp[,1]
exp = exp[,-1]

expr_combat <- ComBat(dat = exp, batch = cohort$projectid,ref.batch = "tcga") %>% as.data.frame()

cols <- ggsci::pal_nejm()(3)
names(cols) <- unique(cohort$projectid)

sample_cols <- cols[cohort$projectid]

boxplot(exp,outline = FALSE, notch = FALSE, las = 2, xaxt = "n", col = sample_cols, main = "Batch effect correction before") 
export::graph2pdf(file = "Batch effect correction before.pdf", width = 10 ,height = 6)
export::graph2png(file = "Batch effect correction before.png", width = 10 ,height = 6)

boxplot(expr_combat,outline = FALSE, notch = FALSE, las = 2, xaxt = "n", col = sample_cols, main = "Batch effect correction after")
export::graph2pdf(file = "Batch effect correction after.pdf", width = 10 ,height = 6)
export::graph2png(file = "Batch effect correction after.png", width = 10 ,height = 6)


colnames(expr_combat) = gsub("\\.", "-", colnames(expr_combat))
tcga = expr_combat[,cohort$sample[cohort$projectid == "tcga"]]
GSE16091 = expr_combat[,cohort$sample[cohort$projectid == "GSE16091"]]
GSE21257 = expr_combat[,cohort$sample[cohort$projectid == "GSE21257"]]



saveRDS(tcga,"./Rdata/exprSet_correct_TARGET_tpm.rds")

saveRDS(GSE16091,"./Rdata/exprSet_correct_GSE16091_tpm.rds")
saveRDS(GSE21257,"./Rdata/exprSet_correct_GSE21257_tpm.rds")


