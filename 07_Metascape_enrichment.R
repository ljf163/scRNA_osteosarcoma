#setwd("./metascape_dotplot/")
options(stringsAsFactors = F)
if (!require(tidyverse)) {
  install.packages("tidyverse")
  require(tidyverse)
}
if (!require(stringr)) {
  install.packages("stringr")
  require(stringr)
}
if (!require(patchwork)) {
  install.packages("patchwork")
  require(patchwork)
}
if (!require(ggrepel)) {
  install.packages("ggrepel")
  require(ggrepel)
}
if (!require(RColorBrewer)) {
  install.packages("RColorBrewer")
  require(RColorBrewer)
}
if (!require(BuenColors)) {
   devtools::install_github("caleblareau/BuenColors")
   require(BuenColors)
}

rm(list = ls())

project_dir = getwd()
result_dir = paste0(project_dir, "/results/GSEA/")
if(!dir.exists(result_dir)){
  dir.create(result_dir)
}


metascape_dir <- "./results/GSEA/"
################################################################################
data <- read.csv(paste(metascape_dir,"/GO_AllLists.csv",sep = ""))
## plot1
data_p1 <- data[,c("Category","Description","Log.q.value.","Enrichment")]
colnames(data_p1) = c("Classification", "Pathways", "logP.value", "Enrichment")
dt = data_p1
dt$Classification[str_detect(dt$Classification,"KEGG")] <- "KEGG"
dt$Classification[str_detect(dt$Classification,"GO Biological Processes")] <- "GO:BP"
dt$Classification[str_detect(dt$Classification,"GO Molecular Functions")] <- "GO:MF"
dt$Classification[str_detect(dt$Classification,"GO Cellular Components")] <- "GO:CC"
dt$Classification[str_detect(dt$Classification,"Hallmark")] <- "Hallmark"
dt$Classification[str_detect(dt$Classification,"WikiPathways")] <- "Wiki"
dt$Classification[str_detect(dt$Classification,"Reactome")] <- "Reactome"
dt$Classification[str_detect(dt$Classification,"Oncogenic")] <- "Onco-Sig"
dt$Classification[str_detect(dt$Classification,"Transcription Factor Targets")] <- "TF"
dt$Classification[str_detect(dt$Classification,"DisGeNET")] <- "DisGeNET"
print(unique(dt$Classification))

dt = dt[order(dt$Classification,dt$logP.value),]
dt = dt[!(duplicated(dt$Pathways)),]
dt$Pathways = factor(dt$Pathways, levels =  dt$Pathways)

dt_s <- dt %>% group_by(Classification) %>% top_n(2)
p1 <- ggplot(dt_s, aes( x = Pathways,
                      y = logP.value, 
                      fill = Classification))+
  scale_fill_manual(values = c('#ee4c58','#56c1ab','#80c5d8','#437eb8','#c4a2c5','#e4e932','#00aeff','#00ad45'))+ 
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 45) )+
  geom_bar(stat = 'identity')+
  coord_flip()+
  ylab('log(q.value)')+
  xlab('') + theme(panel.grid = element_blank(),panel.background = element_blank())
p1
p2 <- ggplot(dt_s, aes( x = Pathways,
                      y = Enrichment,fill = Classification)) +
  scale_fill_manual(values = c('#ee4c58','#56c1ab','#80c5d8','#437eb8','#c4a2c5','#e4e932','#00aeff','#00ad45')) +
  scale_y_continuous(expand = c(0,0)) +
  geom_bar(stat = 'identity')+
  coord_flip()+
  ylab('Enrichment Score')+
  xlab('') +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank())
p2

p1 + p2 + plot_layout(guides = 'collect') + plot_layout(widths = c(2,2))
ggsave(paste0(result_dir, "metascape_barplot.pdf"),width = 13,height = 11)
ggsave(paste0(result_dir, "metascape_barplot.png"),width = 13,height = 11)

###############################################################################
col = c("#E41A1C","#377EB8","#4DAF4A", "#984EA3","#F29403", "#F781BF", "#BC9DCC", "#A65628")
dt_s <- dt
label_data <- dt_s
label_data$label <- NA

# label_data_t <- label_data %>% group_by(Classification) %>% top_n(n = 3, wt = -logP.value)
label_data_t <- label_data[rownames(label_data) %in% c(24,50, 31,77, 8,28,45,15,54,33),]
label_data_t$label <- label_data_t$Pathways

# index = which(label_data$Pathways %in% label_data_t$Pathways)
# label_data$label[index] <- as.character(label_data_t$Pathways)

table(dt_s$Classification)
#dt_title <- data.frame(Classification = unique(dt$Classification),num = as.vector(table(dt_s$Classification)))
dt_s$Classification <- factor(dt_s$Classification,levels = unique(dt_s$Classification))

ggplot(dt_s) +  geom_jitter(aes(Classification, -logP.value, color = Classification,size = -logP.value,alpha = -logP.value), width = 0.45) + pretty_plot() +
  geom_tile(aes(Classification, 0, fill = Classification),
            height=1,
            color = "white",
            alpha = 0.5,
            show.legend = F,
            width= 1) +
  geom_text(data = dt_s[!duplicated(dt_s$Classification), ],
            aes(Classification, 0, label = Classification),
            size = 0,
            color ="black",vjust = 0.5,nudge_y = 100)+
  geom_text_repel(
    data = label_data_t,
    aes(Classification, -logP.value, label = label,color = Classification),
    size = 3.5, max.overlaps = 500#,col = "#eb5424"
  ) +
  xlab("")+
  ylab("-log10(q value)")  +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(labels = unique(dt_s$Classification)) +
  scale_fill_manual(values = col) +
  scale_color_manual(values = col) +
  theme_classic() + ylim(c(0,12)) + theme_bw() +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5),legend.position = "none",plot.background = element_blank(),panel.grid = element_blank())
ggsave(paste0(result_dir, "metascape_mantleplot.pdf"),width = 12,height = 8)
ggsave(paste0(result_dir, "metascape_mantleplot.png"),width = 12,height = 8)


