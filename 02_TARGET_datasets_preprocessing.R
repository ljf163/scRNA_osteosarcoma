rm(list = ls())
library(dplyr)

### count
exprSet = read.table("./datas/TARGET-OS.star_counts.tsv.gz", sep = "\t", header = T)
exprSet = exprSet[-c(1:4),]
colnames(exprSet)[1] = "gene_ID"

ensembl2symbol = read.table("./datas/gencode.v22.annotation.gene.probeMap", sep = "\t", header = T)

identical(exprSet$gene_ID, ensembl2symbol$id)

exprSet$gene_ID = ensembl2symbol$gene


exprSet=t(sapply(split(exprSet,exprSet$gene_ID),function(x) colMeans(x[,2:(ncol(x))])))
exprSet = as.data.frame(exprSet)
exprSet = round(2^exprSet - 1)

colnames(exprSet) = gsub("\\.", "-", colnames(exprSet))

saveRDS(exprSet, file = "./Rdata/exprSet_TARGET_OS_count.rds")


#------------FPKM

exprSet = read.table("./datas/TARGET-OS.star_counts.tsv.gz", sep = "\t", header = T)
exprSet = exprSet[-c(1:4),]

exprSet = exprSet %>% tibble::remove_rownames() %>% tibble::column_to_rownames("Ensembl_ID")
exprSet = 2^exprSet - 1
rownames(exprSet) = gsub("\\..*", "", rownames(exprSet))


geneid_length = readRDS("d:/Datas/TCGA_RNA_data_STAR/geneid_length.rds")
geneid_length$geneid = gsub("\\..*", "", geneid_length$geneid)

exp_count = exprSet[rownames(exprSet) %in% geneid_length$geneid,]
geneid_length = geneid_length[match( rownames(exp_count), geneid_length$geneid),]

identical(geneid_length$geneid, rownames(exp_count))
exp_fpkm <- apply(X = exp_count, 
                  MARGIN = 2, 
                  FUN = function(x) {
                    10^9 * x / geneid_length$efflen / sum(as.numeric(x))
                  })
exp_fpkm = as.data.frame(round(exp_fpkm, 4))

ensembl2symbol = read.table("./datas/gencode.v22.annotation.gene.probeMap", sep = "\t", header = T)
ensembl2symbol$id = gsub("\\..*", "", ensembl2symbol$id)

exp_fpkm = exp_fpkm  %>% tibble::rownames_to_column("id")

exp_fpkm = merge(ensembl2symbol, exp_fpkm, by = "id")
exp_fpkm = exp_fpkm[,-c(1,3:6)]
exp_fpkm=t(sapply(split(exp_fpkm,exp_fpkm$gene),function(x) colMeans(x[,2:(ncol(x))])))
exp_fpkm = as.data.frame(exp_fpkm)

colnames(exp_fpkm) = gsub("\\.", "-", colnames(exp_fpkm))

saveRDS(exp_fpkm, file = "./Rdata/exprSet_TARGET_OS_fpkm.rds")



exp_fpkm=exp_fpkm[rowMeans(exp_fpkm)>0,]
fpkmToTpm <- function(fpkm)
{ 
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
exp_tpm <- apply(exp_fpkm, 2, fpkmToTpm)
colSums(exp_tpm)
exp_tpm = as.data.frame(exp_tpm)
saveRDS(exp_tpm, file = "./Rdata/exprSet_TARGET_OS_tpm.rds")


#---------------------------surv
rm(list = ls())
surv = read.table("./datas/TARGET-OS.survival.tsv", sep = "\t", header = T)
surv = surv[,-3]
colnames(surv) = c("sampleID", "OS", "OS.time")

exp_count = readRDS("./Rdata/exprSet_TARGET_OS_count.rds")
exp_count = exp_count[,colnames(exp_count) %in% surv$sampleID]

surv = surv[match(colnames(exp_count), surv$sampleID),]

exp_fpkm = readRDS("./Rdata/exprSet_TARGET_OS_fpkm.rds")
exp_fpkm = exp_fpkm[,match(surv$sampleID, colnames(exp_fpkm))]

exp_tpm = readRDS("./Rdata/exprSet_TARGET_OS_tpm.rds")
exp_tpm = exp_tpm[,match(surv$sampleID, colnames(exp_tpm))]

identical(colnames(exp_count), surv$sampleID)
identical(colnames(exp_fpkm), surv$sampleID)
identical(colnames(exp_tpm), surv$sampleID)


saveRDS(exp_count, file = "./Rdata/exprSet_TARGET_OS_count.rds")
saveRDS(exp_fpkm, file = "./Rdata/exprSet_TARGET_OS_fpkm.rds")
saveRDS(exp_tpm, file = "./Rdata/exprSet_TARGET_OS_tpm.rds")
saveRDS(surv, file = "./Rdata/survData_TARGET_OS.rds")


#----------------------------------
rm(list = ls())

pd = read.table("./datas/TARGET-OS.clinical.tsv.gz", header = T , sep = "\t")
surv = readRDS("./Rdata/survData_TARGET_OS.rds")

pd = pd[pd$sample_id %in% surv$sampleID,]
pd = pd[,c(1,3,6,14,17)]
colnames(pd) = c("sampleID", "Gender", "Age", "Metastasis", "TumorSite")

pd$Age = round(as.numeric(pd$Age) / 365)
pd$Metastasis = str_replace(pd$Metastasis, " .*", "")
pd$Gender = str_replace(pd$Gender, " .*", "")
pd$TumorSite = str_replace(pd$TumorSite, " .*", "")


saveRDS(pd, file = "./Rdata/phenoData_TARGET_OS.rds")


