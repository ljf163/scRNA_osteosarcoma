# CIBERSORT

rm(list = ls())
setwd("D:/Orders/OS_PD/bulkRNAseq")
exp_file =readRDS("./Rdata/exprSet_TARGET_OS_tpm.rds")

project_dir = getwd()
result_dir = paste0(project_dir, "/results/ImmuneInfiltration/")
if(!dir.exists(result_dir)){
  dir.create(result_dir)
}



write.table(exp_file,file = "./scripts/cibersort/exp_file.txt",row.names = T,col.names = T,quote = F,sep = "\t")

source("./scripts/cibersort/CIBERSORT.R")

# Define LM22 file
LM22_file <- "./scripts/cibersort/LM22.txt"
exp_file <- "./scripts/cibersort/exp_file.txt"

cibersort_output = CIBERSORT(LM22_file, exp_file, perm = 1000, QN = F)
saveRDS(cibersort_output,file = "./Rdata/cibersort_result.rds")
write.csv(cibersort_output, file = paste0(result_dir, "cibersort_results.csv"),quote =F)
print("CIBERSOTR analysis completed!")


print("Start ploting!")

# rm(list = ls())
library(reshape2)
library(ggpubr)
cibersort_output = readRDS("./Rdata/cibersort_result.rds")
cbsOut = cibersort_output[,1:22]

rs = readRDS("./Rdata/RS_TARGET.rds")

rs = rs[match(rownames(cibersort_output), rownames(rs)),]
identical(rownames(rs), rownames(cibersort_output))
group = rs$group

data = as.data.frame(cbsOut)

library(tidyverse)
data2 <- data %>% 
  add_column(class = group, .before = 1) %>% 
  gather("cell_types","percentages", -class) %>% 
  mutate(cell_types = factor(cell_types, levels = colnames(data))) %>% 
  mutate(class = factor(class))
data3 <- data2 %>% 
  group_by(class, cell_types) %>% 
  summarise("Mx" = max(percentages)) %>% 
  group_by(cell_types) %>% 
  summarise(lable = max(Mx))

x = c()
xend = c()
for(i in 1:length(unique(data2$cell_types))){
  x[i] = i - 0.2
  xend[i] = i + 0.2}
y = data3$lable +0.01
yend = data3$lable +0.01
line <-  data.frame(x, xend, y, yend)

library(ggthemes)

ggplot()+
  geom_violin(data = data2, aes(x = cell_types, y = percentages,color = class),
              position = position_dodge(1),scale = "width",show.legend = F)+
  geom_boxplot(data = data2, aes(cell_types, percentages,fill = class),
               position = position_dodge(1),width=.4, outlier.shape = NA)+
  stat_summary(data = data2, aes(cell_types, percentages,fill = class),
               fun.y = "median",geom = "point",shape = 16, 
               size = 2, color = "black",position = position_dodge(1))+
  ggsci::scale_fill_npg()+
  ggsci::scale_color_npg()+
  # scale_color_manual(values = c("#c62d17", "#266b69"))+
  # scale_fill_manual(values = c("#c62d17", "#266b69"))+
  ggpubr::stat_compare_means(data = data2,aes(cell_types, percentages,group = class),
                             label = "p.signif",
                             method = "wilcox.test",label.y = data3$lable + 0.03)+
  geom_segment(data = line, aes(x=x, y=y, xend= xend, yend=yend))+
  labs(x = "", fill = "",color = "", y = "Percentage")+
  theme_few(base_size = 14)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, colour = "black"),
        legend.position = "top")


ggsave(paste0(result_dir, "CIBERSORT_violinplot.pdf"),width = 10,height = 6)
ggsave(paste0(result_dir, "CIBERSORT_violinplot.png"),width = 10,height = 6)
