library(GEOquery)
library(stringr)
setwd("./datas/")

#---------------------- GSE234998---------------------
# RNA-seq
rm(list=ls())
gset = getGEO('GSE234998', destdir=".",getGPL = F)
gset=gset[[1]]
gset@annotation # GPL25318

exprSet = read.table("./GSE234998_rpkm_rna.txt.gz", header = T, sep = "\t")
exprSet = exprSet[,c(15,1:13)]
exprSet=t(sapply(split(exprSet,exprSet$symbol),function(x) colMeans(x[,2:(ncol(x))])))
exprSet = as.data.frame(exprSet)
colnames(exprSet) = gsub("\\.", "-", colnames(exprSet))

exp_fpkm=exprSet[rowMeans(exprSet)>0,]
fpkmToTpm <- function(fpkm)
{ 
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
exp_tpm <- apply(exp_fpkm, 2, fpkmToTpm)
exp_tpm = as.data.frame(exp_tpm)
colSums(exp_tpm)



pd = pData(gset)
pd = pd[, c(1,43)]
colnames(pd) = c("sampleID", "state")
pd$sampleID = str_replace(pd$sampleID, ".*_TCCC", "TCCC")
pd$sampleID = str_replace(pd$sampleID, "_.*", "")

exprSet = exprSet[, match(pd$sampleID, colnames(exprSet))]
identical(colnames(exprSet), pd$sampleID)

exp_tpm = exp_tpm[, match(pd$sampleID, colnames(exp_tpm))]
identical(colnames(exp_tpm), pd$sampleID)

saveRDS(exprSet, file = "../Rdata/exprSet_GSE234998_rpkm.rds")
saveRDS(exp_tpm, file = "../Rdata/exprSet_GSE234998_tpm.rds")
saveRDS(pd, file = "../Rdata/phenoData_GSE234998.rds")




#---------------------- GSE87624 ---------------------
# RNA-seq
rm(list=ls())
gset = getGEO('GSE87624', destdir=".",getGPL = F)
gset=gset[[1]]
gset@annotation # GPL570

exprSet = read.table("./GSE87624_Human_masked.txt.gz", header = T, sep = "\t")
exprSet=t(sapply(split(exprSet,exprSet$X),function(x) colMeans(x[,2:(ncol(x))])))
exprSet = as.data.frame(exprSet)
boxplot(exprSet)

exp_fpkm=exprSet[rowMeans(exprSet)>0,]
fpkmToTpm <- function(fpkm)
{ 
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
exp_tpm <- apply(exp_fpkm, 2, fpkmToTpm)
exp_tpm = as.data.frame(exp_tpm)
colSums(exp_tpm)

pd = pData(gset)
pd = pd[, c(1,40)]
colnames(pd) = c("sampleID", "state")
pd = pd[pd$state %in% c("metastasis", "normal bone", "primary"),]
pd = pd[order(pd$state),]

exprSet = exprSet[, match(pd$sampleID, colnames(exprSet))]
identical(colnames(exprSet), pd$sampleID)

exp_tpm = exp_tpm[, match(pd$sampleID, colnames(exp_tpm))]
identical(colnames(exp_tpm), pd$sampleID)

saveRDS(exprSet, file = "../Rdata/exprSet_GSE87624_fpkm.rds")
saveRDS(exp_tpm, file = "../Rdata/exprSet_GSE87624_tpm.rds")
saveRDS(pd, file = "../Rdata/phenoData_GSE87624.rds")


#---------------------- GSE33382 ---------------------
# RNA-seq
rm(list=ls())
gset = getGEO('GSE33382', destdir=".",getGPL = F)
gset=gset[[1]]
gset@annotation # GPL10295

exprSet = exprs(gset)
# exprSet = log2(exprSet + 1)
boxplot(exprSet)

gpl <- getGEO('GPL10295')
ids=Table(gpl)
colnames(ids)
ids = ids[,c(1,6)]
colnames(ids) = c("probe_id", "symbol")
ids$symbol = str_replace(ids$symbol, " .*","")
ids[ids == ""] = NA
ids = na.omit(ids)

exp = data.frame(probe_id = rownames(exprSet), exprSet)
exp = merge(ids, exp, by = "probe_id")
exprSet = exp[,-1]
exprSet=t(sapply(split(exprSet,exprSet$symbol),function(x) colMeans(x[,2:(ncol(x))])))
exprSet = limma::normalizeBetweenArrays(exprSet)
exprSet = as.data.frame(exprSet)
boxplot(exprSet)

pd = pData(gset)
identical(colnames(exprSet), pd$geo_accession)

pd = pd[,c(2,42,44)]

colnames(pd) = c("sampleID", "Grade", "Metastasis")


saveRDS(exprSet, file = "../Rdata/exprSet_GSE33382.rds")
saveRDS(pd, file = "../Rdata/phenoData_GSE33382.rds")


#---------------------- GSE21257  ---------------------
# RNA-seq
rm(list=ls())
gset = getGEO('GSE21257', destdir=".",getGPL = F)
gset=gset[[1]]
gset@annotation # GPL10295

exprSet = exprs(gset)
# exprSet = log2(exprSet + 1)
boxplot(exprSet)

gpl <- getGEO('GPL10295')
ids=Table(gpl)
colnames(ids)
ids = ids[,c(1,6)]
colnames(ids) = c("probe_id", "symbol")
ids$symbol = str_replace(ids$symbol, " .*","")
ids[ids == ""] = NA
ids = na.omit(ids)

exp = data.frame(probe_id = rownames(exprSet), exprSet)
exp = merge(ids, exp, by = "probe_id")
exprSet = exp[,-1]
exprSet=t(sapply(split(exprSet,exprSet$symbol),function(x) colMeans(x[,2:(ncol(x))])))
exprSet = limma::normalizeBetweenArrays(exprSet)
exprSet = as.data.frame(exprSet)
boxplot(exprSet)

pd = pData(gset)
identical(colnames(exprSet), pd$geo_accession)
# table(pd$source_name_ch1)


surv = pd[,c(2,42)]
surv$OS = str_replace(surv$`status:ch1`, " .*", "")
surv$OS = ifelse(surv$OS == "Alive", 0, 1)
surv$OS.time = str_replace(surv$`status:ch1`, "[A-Za-z]* [a-z]* ", "")
surv$OS.time = str_replace(surv$OS.time, " .*", "")
surv$OS.time = as.numeric(surv$OS.time) * 30

surv = surv[,3:4]

identical(colnames(exprSet), rownames(surv))

saveRDS(exprSet, file = "../Rdata/exprSet_GSE21257.rds")
saveRDS(surv, file = "../Rdata/survData_GSE21257.rds")
saveRDS(pd, file = "../Rdata/phenoData_GSE21257.rds")




#---------------------- GSE16091  ---------------------
# RNA-seq
rm(list=ls())
gset = getGEO('GSE16091', destdir=".",getGPL = F)
gset=gset[[1]]
gset@annotation # GPL96

exprSet = exprs(gset)
# exprSet = log2(exprSet + 1)
boxplot(exprSet)

library(AnnoProbe)
ids=idmap('GPL96',type = 'soft')
colnames(ids) = c("probe_id", "symbol")
ids$symbol = str_replace(ids$symbol, " .*","")
ids[ids == ""] = NA
ids = na.omit(ids)

exp = data.frame(probe_id = rownames(exprSet), exprSet)
exp = merge(ids, exp, by = "probe_id")
exprSet = exp[,-1]
exprSet=t(sapply(split(exprSet,exprSet$symbol),function(x) colMeans(x[,2:(ncol(x))])))
exprSet = limma::normalizeBetweenArrays(exprSet)
exprSet = as.data.frame(exprSet)
boxplot(exprSet)

pd = pData(gset)
identical(colnames(exprSet), pd$geo_accession)
# table(pd$source_name_ch1)


surv = pd[,c(34,35)]
colnames(surv) = c("OS", "OS.time")
surv$OS = ifelse(surv$OS == "D", 1, 0)
surv$OS.time = as.numeric(surv$OS.time)

identical(colnames(exprSet), rownames(surv))


saveRDS(exprSet, file = "../Rdata/exprSet_GSE16091.rds")
saveRDS(surv, file = "../Rdata/survData_GSE16091.rds")
saveRDS(pd, file = "../Rdata/phenoData_GSE16091.rds")


#---------------------- GSE200845  ---------------------
# RNA-seq
rm(list=ls())
gset = getGEO('GSE200845', destdir=".",getGPL = F)
gset=gset[[1]]
gset@annotation # GPL16791

pd = pData(gset)

pd = pd[pd$`chemical name:ch1` == "Etoposide",]
pd = na.omit(pd)

pd = pd[,c(49,50,52,54,55,56)]
pd$sampleID = paste0(pd$`assay plate:ch1`, "_", pd$`assay plate well:ch1`)
colnames(pd) = gsub("\\:.*","", colnames(pd))


exprSet = read.csv("./GSE200845_httr_u2os_pilot_GEO_count_data.csv.gz")

exprSet = exprSet[exprSet$sample_id %in% pd$sampleID,]
exprSet = exprSet[,-c(2:6)]
exprSet = exprSet %>%  tibble::remove_rownames() %>%  tibble::column_to_rownames( "sample_id") %>%
  t() %>% as.data.frame() %>% tibble::rownames_to_column("probe_id")

ids = read.csv("./GSE200845_hWTv2_u2os_probes.csv.gz")
ids = ids[,c(2,7)]
colnames(ids) = c("probe_id", "symbol")

tmp = merge(ids, exprSet, by = "probe_id")
tmp = tmp[,-1]

exprSet=t(sapply(split(tmp,tmp$symbol),function(x) colMeans(x[,2:(ncol(x))])))
exprSet = as.data.frame(round(exprSet))

pd = pd[order(pd$`dose level`,pd$sampleID),]
exprSet = exprSet[, match(pd$sampleID, colnames(exprSet))]

saveRDS(exprSet, file = "../Rdata/exprSet_GSE200845.rds")
saveRDS(pd, file = "../Rdata/phenoData_GSE200845.rds")


exp = data.frame(Name = rownames(exprSet), Description = "", exprSet)
write.table(exp,  file = "./exp_2_GSEA_GSE200845.txt", row.names = F, quote = F, sep = "\t")
