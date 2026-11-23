rm(list = ls())


project_dir = getwd()
result_dir = paste0(project_dir, "/results/ImmuneTherapy/")
if(!dir.exists(result_dir)){
  dir.create(result_dir)
}

library(stringr)
library(dplyr)

exprSet = readRDS("./Rdata/exprSet_TARGET_OS_tpm.rds")

exp_tide = t(apply(exprSet, 1, function(x){x-(mean(x))})) %>% as.data.frame()
write.table(exp_tide,file = paste0(result_dir, "tide_input.txt") ,sep = '\t',quote = F)


tideScore = read.csv("./results/ImmuneTherapy/TIDE_output.csv",header = T)
colnames(tideScore)[1] = 'sampleID'

risk_info =readRDS("./Rdata/RS_TARGET.rds")
risk_info = tibble::rownames_to_column(risk_info, "sampleID")
rt = merge(risk_info, tideScore,  by = "sampleID")
rt = rt[order(rt$group),]

write.csv(rt ,file =  paste0(result_dir, "TIDE_score.csv"), row.names = F, quote = F)


library(ggsci)
library(smplot2)
library(ggpubr)

my_comparisons = list(c("High",'Low'))

p = ggviolin(rt, x = "group", y = "TIDE",xlab = "Risk Group",legend.title = "",ylab = "TIDE",
             fill = "group", palette = "npg",
             add = "boxplot", add.params = list(fill = "white")) 
p + stat_compare_means(comparisons = my_comparisons, label = "p.signif", label.y = 4.2) +
  stat_compare_means(label.x = 0.5,label.y = 4.2)+ theme_bw() + theme(panel.grid=element_blank(),legend.position = "none")

ggsave(paste0(result_dir, "TIDE_violinPlot.pdf"),width = 7,height = 6)
ggsave(paste0(result_dir, "TIDE_violinPlot.png"),width = 7,height = 6)


# Dysfunction
p = ggviolin(rt, x = "risk", y = "Dysfunction",legend.title = "Risk",
             fill = "risk", palette =c("#DF8F44FF","#374E55FF"),
             add = "boxplot", add.params = list(fill = "white"))
p + stat_compare_means(comparisons = my_comparisons, label = "p.signif",label.y = 3.9) +
  stat_compare_means(label.y = 3.9,label.x = 0.5)

ggsave('./results/tide/Dysfunction_violinPlot.pdf',width = 8,height = 6)

# Exclusion
p = ggviolin(rt, x = "risk", y = "Exclusion",legend.title = "Risk",
             fill = "risk", palette =c("#DF8F44FF","#374E55FF"),
             add = "boxplot", add.params = list(fill = "white"))
p + stat_compare_means(comparisons = my_comparisons, label = "p.signif",label.y = 3.6) +
  stat_compare_means(label.y = 3.6,label.x = 0.5)

ggsave('./results/tide/Exclusion_violinPlot.pdf',width = 8,height = 6)



# 
p = ggboxplot(rt, x = "IRGPI group", y = "MSI.Expr.Sig",
              color = "IRGPI group", palette =c( "#FC4E07","blue"),
              add = "jitter", shape = "IRGPI group")
p + stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  stat_compare_means(label.y = 1,label.x = 0.5)





