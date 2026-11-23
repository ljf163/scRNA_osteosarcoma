# HALLMARK GSVA

rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(clusterProfiler)
library(msigdbr)  #install.packages("msigdbr")
library(GSVA) 
library(GSEABase)
library(pheatmap)
library(limma)
library(BiocParallel)

project_dir = getwd()
result_dir = paste0(project_dir, "/results/GSEA/")
if(!dir.exists(result_dir)){
  dir.create(result_dir)
}

tb = read.table("./results/OS.PD.GSEA.Gsea.1713501712648/OSTEOSARCOMA DISEASE PROGRESSION GENE SET.tsv", header = T, sep = "\t")
core.genes = tb$SYMBOL[tb$CORE.ENRICHMENT == "Yes"]


# gmt_files = "D:/ProjectR/msigdb/h.all.v2023.2.Hs.symbols.gmt"
# genelist <- getGmt(gmt_files)
genelist = list(core.genes)
names(genelist) = "Osteosarcoma Disease Progression Gene Set"
dat = readRDS("./Rdata/exprSet_GSE200845.rds")
dat = as.matrix(dat)

gsva_mat <- gsva(expr=dat, 
                 gset.idx.list=genelist, 
                 kcdf="Poisson" ,#"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
                 verbose=T, 
                 parallel.sz = parallel::detectCores())
gsva_mat = t(gsva_mat)
saveRDS(gsva_mat, file = "./Rdata/gsva_result.rds")



pd = readRDS("./Rdata/phenoData_GSE200845.rds")
rt =cbind(pd, gsva_mat)

write.csv(rt, file =  paste0(result_dir, "GSVA_score_OS.PD.geneset.csv"), quote =F)


summarySE<- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                     conf.interval=.95,.drop=TRUE) {
  library(plyr)
  
  # New version of length which can handleNA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else      length(x)
  }
  
  # This does the summary. For each group'sdata frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars,.drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   =sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column   
  datac <- rename(datac, c("mean"= measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standarderror
  # Calculate t-statistic for confidenceinterval:
  # e.g., if conf.interval is .95, use .975(above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5,datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


colnames(rt)[c(3,8)] = c("concentration","GSVA_Score")
carss<-summarySE(rt, measurevar="GSVA_Score",groupvars=c("concentration"))

carss$concentration = factor(carss$concentration, levels = c("0.01 uM","0.03 uM","0.1 uM","0.3 uM" ,"1 uM" ,  "3 uM","10 uM"  ))
unique(carss$concentration)


mytheme = theme_bw() +
  theme(axis.text.x = element_text(size = 9),axis.text.y = element_text(size = 9))+
  theme(axis.title.y= element_text(size=12))+theme(axis.title.x = element_text(size = 12))+
  theme(legend.title=element_text(size=5),legend.text=element_text(size=5))+theme(legend.position = "none",panel.grid = element_blank())


ggplot(carss,aes(x=concentration, y=GSVA_Score,colour = "#EA5D2D")) +
  geom_errorbar(aes(ymin=GSVA_Score-se, ymax=GSVA_Score+se),width=.1,lwd =.8) +
  geom_line(group = 1, lwd = .8) +
  geom_point(size =2) + mytheme + ylab("GSVA Score") + xlab("Concentration of Etoposide")

ggsave(paste0(result_dir, "GSVA_Score_curv.pdf"), width = 6, height = 5)
ggsave(paste0(result_dir, "GSVA_Score_curv.png"), width = 6, height = 5)


