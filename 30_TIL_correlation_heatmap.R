


rm(list = ls())


project_dir = getwd()
result_dir = paste0(project_dir, "/results/ImmuneInfiltration/")
if(!dir.exists(result_dir)){
  dir.create(result_dir)
}



exprSet = readRDS("./Rdata/exprSet_TARGET_OS_tpm.rds")
features = readRDS("./Rdata/Features_selected_LASSO.rds")
exp = t(exprSet[features,])

inputData = readRDS("./Rdata/cibersort_result.rds")
inputData = inputData[,1:22]
inputData = inputData[,-21]

nc = cbind(inputData, exp)
nc[1:5,1:5]

library(Hmisc)
res_coef = rcorr(nc)$r[1:ncol(inputData),(ncol(nc)-length(features)+1):ncol(nc)]
res_coef[1:4,1:4]
res_coef = res_coef[,order(colnames(res_coef))]

res_p =  rcorr(nc)$P[1:ncol(inputData),(ncol(nc)-length(features)+1):ncol(nc)]
res_p[1:4,1:4]
res_p = res_p[,order(colnames(res_p))]

library(dplyr)
sig = res_p
for (i in 1:ncol(sig)) {
  sig[,i] = case_when(sig[,i] < 0.001 ~ "***",
                      sig[,i] < 0.01 ~ "**",
                      sig[,i] < 0.05 ~ "*",
                      .default =  "")
}
source("./scripts/modified_pheatmap.R")
p = pheatmap(res_coef,
         display_numbers =sig,
         angle_col =45,cluster_cols = F,
         color = colorRampPalette(c("#92b7d1", "white", "#d71e22"))(100),
         border_color = "white",
         treeheight_col = 0,
         treeheight_row = 0)
export::graph2pdf(file = paste0(result_dir, "Cor_ModelGenes_CIBERSORT.pdf"), width = 9, height = 7)
export::graph2png(file = paste0(result_dir, "Cor_ModelGenes_CIBERSORT.png"), width = 9, height = 7)


# rm(list = ls())

exprSet = readRDS("./Rdata/exprSet_TARGET_OS_tpm.rds")
features = readRDS("./Rdata/Features_selected_LASSO.rds")
exp = t(exprSet[features,])

inputData = read.csv("./results/ImmuneInfiltration/CIBERSORTx_Job4_Results.csv", row.names = 1)
inputData = inputData[,1:11]

nc = cbind(inputData, exp)
nc[1:5,1:5]

library(Hmisc)
res_coef = rcorr(nc)$r[1:ncol(inputData),(ncol(nc)-length(features)+1):ncol(nc)]
res_coef[1:4,1:4]
res_coef = res_coef[,order(colnames(res_coef))]

res_p =  rcorr(nc)$P[1:ncol(inputData),(ncol(nc)-length(features)+1):ncol(nc)]
res_p[1:4,1:4]
res_p = res_p[,order(colnames(res_p))]

library(dplyr)
sig = res_p
for (i in 1:ncol(sig)) {
  sig[,i] = case_when(sig[,i] < 0.001 ~ "***",
                      sig[,i] < 0.01 ~ "**",
                      sig[,i] < 0.05 ~ "*",
                      .default =  "")
}
source("./scripts/modified_pheatmap.R")
p = pheatmap(res_coef,
             display_numbers =sig,
             angle_col =45,cluster_cols = F,
             color = colorRampPalette(c("#92b7d1", "white", "#d71e22"))(100),
             border_color = "white",
             treeheight_col = 0,
             treeheight_row = 0)
export::graph2pdf(file = paste0(result_dir, "Cor_ModelGenes_CIBERSORTx.pdf"), width = 9, height = 6)
export::graph2png(file = paste0(result_dir, "Cor_ModelGenes_CIBERSORTx.png"), width = 9, height = 6)

