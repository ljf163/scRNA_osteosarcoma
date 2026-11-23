rm(list = ls())

library(dplyr)

exprSet = readRDS("./Rdata/exprSet_correct_TARGET_tpm.rds")
# exprSet = readRDS("./Rdata/exprSet_TCGA_TPM_tumor.rds")
degs = readRDS("./Rdata/os_realated_genes.rds")
surv = readRDS("./Rdata/survData_TARGET_OS.rds")
surv = tibble::remove_rownames(surv)
surv = tibble::column_to_rownames(surv, "sampleID")

exp = exprSet[rownames(exprSet) %in% degs,]
exp = na.omit(exp)
exp = exp %>% t() %>% as.data.frame()

identical(rownames(exp), rownames(surv))

rt = cbind(surv, exp)


pFilter=0.05;
train = rt

library(survival)
library(survminer)

project_dir = getwd()
result_dir = paste0(project_dir, "/results/MachineLearning/")
if(!dir.exists(result_dir)){
  dir.create(result_dir)
}

#############uniCOX#############
outTab = data.frame()
pFilter=0.05
sigGenes = c()
for (i in colnames(train[, 3:ncol(train)])) {
  if (sd(train[, i]) < 0.000001) {
    next
  }
  
  
  cox <- coxph(Surv(OS.time, OS) ~ train[, i], data = train)
  coxSummary = summary(cox)
  coxP = coxSummary$coefficients[, "Pr(>|z|)"]
  outTab = rbind(
    outTab,
    cbind(
      id = i,
      # KM = pValue,
      B = coxSummary$coefficients[, "coef"],
      SE = coxSummary$coefficients[, "se(coef)"],
      HR = coxSummary$conf.int[, "exp(coef)"],
      HR.95L = coxSummary$conf.int[, "lower .95"],
      HR.95H = coxSummary$conf.int[, "upper .95"],
      pvalue = coxP
    )
  )
  if ((coxP < pFilter)) {
    sigGenes = c(sigGenes, i)
  }
}

saveRDS(sigGenes, file = "./Rdata/Unicox_sig_genes.rds")
write.csv(outTab, file = paste0(result_dir, "Unicox_result.csv"),row.names = F, quote = F)
write.table(sigGenes, file = paste0(result_dir, "Unicox_SigGenes.txt"),row.names = F, quote = F,col.names = F)
