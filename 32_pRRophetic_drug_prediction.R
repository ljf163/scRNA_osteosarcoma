options("repos"=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

if(!requireNamespace("car",quietly = TRUE)) install.packages("car",update = F,ask = F)
if(!requireNamespace("ridge",quietly = TRUE)) install.packages("ridge",update = F,ask = F)
if(!requireNamespace("preprocessCore",quietly = TRUE)) BiocManager::install("preprocessCore",update = F,ask = F)
if(!requireNamespace("genefilter",quietly = TRUE)) BiocManager::install("genefilter",update = F,ask = F)
if(!requireNamespace("sva",quietly = TRUE)) BiocManager::install("sva",update = F,ask = F)


library(pRRophetic)
library(ggplot2)
library(reshape2)
library(ggpubr)
set.seed(12345)

rm(list = ls())
setwd("D:/Orders/OS_PD/bulkRNAseq")

project_dir = getwd()
result_dir = paste0(project_dir, "/results/pRRophetic_result/")
if(!dir.exists(result_dir)){
  dir.create(result_dir)
}


exprSet =readRDS("./Rdata/exprSet_TARGET_OS_tpm.rds")
exprSet = as.matrix(exprSet)

rs = readRDS("./Rdata/RS_TARGET.rds")
rs = rs %>% as.data.frame %>% tibble::rownames_to_column("sampleID")
rs = rs[match(colnames(exprSet), rs$sampleID),]




drugs = c(
  "A.443654", "A.770041", "ABT.263", "ABT.888", "AICAR",
  "AKT.inhibitor.VIII", "AMG.706", "AP.24534", "AS601245", "ATRA",
  "AUY922", "Axitinib", "AZ628", "AZD.0530", "AZD.2281", "AZD6244",
  "AZD6482", "AZD7762", "AZD8055", "BAY.61.3606", "Bexarotene",
  "BI.2536", "BIBW2992", "Bicalutamide", "BI.D1870", "BIRB.0796",
  "Bleomycin", "BMS.509744", "BMS.536924",
  "BMS.754807","Bortezomib", "Bosutinib", "Bryostatin.1", "BX.795", "Camptothecin",
    "CEP.701", "CGP.082996", "CGP.60474",
  "CHIR.99021", "CI.1040", "Cisplatin", "CMK", "Cyclopamine",
  "Cytarabine", "Dasatinib", "DMOG", "Docetaxel", "Doxorubicin",
  "Elesclomol", "Embelin", "Epothilone.B", "Erlotinib",
  "Etoposide", "FH535", "FTI.277", "GDC.0449", "GDC0941", "Gefitinib",
  "Gemcitabine", "GNF.2", "GSK269962A", "GSK.650394", "GW.441756",
  "GW843682X", "Imatinib", "IPA.3",
  "JNK.9L",
  "JNK.Inhibitor.VIII", "JW.7.52.1", "KIN001.135", "KU.55933",
  "Lapatinib", "Lenalidomide", "LFM.A13", "Metformin", "Methotrexate",
  "MG.132", "Midostaurin", "Mitomycin.C", "MK.2206", "MS.275", "Nilotinib",
  "NSC.87877", "NU.7441", "Nutlin.3a", "NVP.BEZ235", "NVP.TAE684",
  "Obatoclax.Mesylate", "OSI.906", "PAC.1", "Paclitaxel", "Parthenolide",
  "Pazopanib", "PD.0325901", "PD.0332991", "PD.173074", "PF.02341066",
   "PF.562271", "PHA.665752", "PLX4720", "Pyrimethamine",
  "QS11", "Rapamycin", "RDEA119", "RO.3306", "Roscovitine", "Salubrinal",
  "SB.216763", "SB590885", "Shikonin", "SL.0101.1", "Sorafenib",
  "S.Trityl.L.cysteine", "Sunitinib", "Temsirolimus", "Thapsigargin",
  "Tipifarnib",
  "Vinblastine", "Vinorelbine", "Vorinostat",
  "VX.680", "VX.702", "WH.4.023",
  "WZ.1.84", "X17.AAG",
  "X681640", "XMD8.85", "Z.LLNle.CHO", "ZM.447439","AG.014699","BMS.708163","CCT007093","CCT018159","EHT.1864", "JNJ.26854165","PF.4708671", "TW.37", "WO2009093972"
  )
out_tab = rs[,c(1,5)]
out_p = data.frame()
for(drugName in drugs) {
  print(paste("-----------", drugName, "---------------"))
  scores <- pRRopheticPredict(exprSet, drugName, "allSolidTumors", selection=1)
  out_tab = data.frame(out_tab,scores)
  result = data.frame(sampleID = names(scores), IC50 = scores)
  result$IC50[result$IC50 == "NaN"] = NA
  result = na.omit(result)
  result = merge(rs[,c(1,ncol(rs))], result, by="sampleID")
  
  score_high = result$IC50[result$group == "High"]
  score_low = result$IC50[result$group == "Low"]
  wt = wilcox.test(score_high, score_low)
  pvalue = wt$p.value
  out_p = rbind(out_p, cbind(drugName, pvalue))
  if(pvalue >= 0.05){next}

  inputData = melt(result)
  inputData = inputData[,-3]
  colnames(inputData) = c("sampleID", "Risk", "IC50")

  my_comparisons = list(c("High", "Low"))
  p<-ggplot(inputData,aes(x=Risk,y=IC50))+
    geom_boxplot(aes(fill=Risk)) +
    ggsci::scale_fill_npg()+
    ggtitle(drugName) + theme_classic()
  p2 = p + stat_compare_means(comparisons = my_comparisons, label = "p.signif",method  = "wilcox.test") +
    stat_compare_means(label.x = 0.5,method = "wilcox.test")

  plotPath = paste("./results/pRRophetic_result/", drugName, ".pdf", sep="")
  pdf(plotPath)
  print(p2)
  dev.off()
}
colnames(out_tab)[3:ncol(out_tab)] = drugs
write.csv(out_tab, file = "./results/pRRophetic_result/result.csv",quote =F)

colnames(out_p) = c("drugName", "pvalue")
write.csv(out_p, file = "./results/pRRophetic_result/pvalue.csv",quote =F)


