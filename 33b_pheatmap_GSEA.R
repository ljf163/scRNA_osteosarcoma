rm(list = ls())
project_dir = getwd()
result_dir = paste0(project_dir, "/results/GSEA/")
if(!dir.exists(result_dir)){
  dir.create(result_dir)
}


tb = read.table("./results/OS.PD.GSEA.Gsea.1713501712648/OSTEOSARCOMA DISEASE PROGRESSION GENE SET.tsv", header = T, sep = "\t")
tb = tb[order(tb$RANK.METRIC.SCORE, decreasing = F),]

core.genes = tb$SYMBOL[tb$CORE.ENRICHMENT == "Yes"][1:50]


exprSet = readRDS("./Rdata/exprSet_GSE200845.rds")
pd = readRDS("./Rdata/phenoData_GSE200845.rds")


print("Heatmap!")
group_sample = data.frame(
  Etoposide = pd$`chemical concentration`,
  sampleID= pd$sampleID
)
# group_sample = group_sample[order(group_sample$group),]

annotation_col <-  data.frame(
  Etoposide = group_sample$Etoposide
)
rownames(annotation_col) = group_sample$sampleID


library(RColorBrewer)
display.brewer.all()


cols = brewer.pal(7,"PuRd")
names(cols) = unique(group_sample$Etoposide)
ann_colors = list(Etoposide = cols)


dat = t(scale(t(log2(exprSet[rownames(exprSet) %in% core.genes,] + 1 ))))
dat = na.omit(dat)
dat[dat > 2] = 2
dat[dat < -2] = -2
dat = dat[,match(rownames(annotation_col), colnames(dat))]
library(pheatmap)
pheatmap(dat, annotation_col = annotation_col, cluster_cols = F,show_colnames = F,
         annotation_colors = ann_colors,fontsize_row = 8,#cellwidth = 25,cellheight = 10,
         border_color = "white")
export::graph2pdf(file = paste0(result_dir, "OS.PD.geneset.top50.heatmap.pdf"), width = 10, height = 6)
export::graph2png(file = paste0(result_dir, "OS.PD.geneset.top50.heatmap.png"), width = 10, height = 6)


#--------------GSEA
library(ggplot2)

p1 <- ggplot(data = tb,aes(x=RANK.IN.GENE.LIST,y=RUNNING.ES,color=RANK.METRIC.SCORE)) +
  scale_color_gradient2(high = 'red',mid = 'pink',low = 'blue',
                        guide = guide_colorbar(barwidth = 8, 
                                               barheight = 1 , 
                                               ticks.colour = 'black', 
                                               ticks.linewidth = 2,
                                               ticks = F, 
                                               direction = 'horizontal', 
                                               title = 'Rank Scores', 
                                               title.position = 'top') ,
  ) +
  geom_point(size=2) + 
  geom_line(size=1) + 
  geom_segment(aes(xend=RANK.IN.GENE.LIST,yend=0)) + 
  theme_classic() +
  ylab('Enrichment score (ES)') +
  xlab('Gene ranked positions') +
  theme(axis.title = element_text(size = 12,face = 'bold'), 
        legend.position = c(0.85,0.95),
        legend.background = element_blank(),
        axis.text = element_text(size = 8,face = 'bold'),
        axis.line = element_line(size = 1),
        axis.ticks.length = unit(0.25,'cm'),
        axis.ticks = element_line(size = 1),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9,face = 'bold'),
        plot.title = element_text(size = 12,face = 'bold')
  ) +
  geom_hline(yintercept = 0,size=1) +
  labs(title = 'GSEA Enrichment plot: Osteosarcoma Disease Progression Gene Set')

p1



features = readRDS("./Rdata/Features_selected_LASSO.rds")
# my_gene <- subset(tb,CORE.ENRICHMENT == 'Yes')
my_gene = tb[tb$SYMBOL %in% features,]

library(ggrepel)


p1 + geom_text_repel(data = my_gene,aes(x=RANK.IN.GENE.LIST,y=RUNNING.ES,label=SYMBOL),
                     color='black', 
                     force = 0.05 ,
                     max.overlaps = 35 , 
                     hjust= -1.5 , 
                     size = 2 ,  
                     point.padding = 0.5 ,
                     fontface='italic', 
                     arrow = arrow(length = unit(0.005, "npc"), 
                                   ends = "last", type = "open"), 
                     segment.color= 'black', 
                     segment.size=0.5 , 
                     segment.alpha= 0.8) 
export::graph2pdf(file = paste0(result_dir, "OS.PD.geneset.GSEA.pdf"), width = 9, height = 6)
export::graph2png(file = paste0(result_dir, "OS.PD.geneset.GSEA.png"), width = 9, height = 6)
