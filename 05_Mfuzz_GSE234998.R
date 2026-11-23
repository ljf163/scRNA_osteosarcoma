
#------------------------GSE234998
rm(list = ls())
library(tidyverse)
library(Mfuzz)

project_dir = getwd()
result_dir = paste0(project_dir, "/results/Mfuzz/")
if(!dir.exists(result_dir)){
  dir.create(result_dir)
}

set.seed(123)

exprSet = readRDS("./Rdata/exprSet_GSE234998_tpm.rds")
exprSet = exprSet[-c(1:25),]
pd = readRDS("./Rdata/phenoData_GSE234998.rds")
table(pd$state)
colnames(exprSet) = pd$state

library(limma)
avereps_df<- t(limma::avereps(t(exprSet), ID = colnames(exprSet)))
View(avereps_df)
colnames(avereps_df)
avereps_df = avereps_df[,colnames(avereps_df)[c(2,4,3)]]


gene_tpm <- data.matrix(avereps_df)
eset <- new("ExpressionSet",exprs = gene_tpm)
dim(eset) 

gene.r <- filter.NA(eset, thres=0.25)

gene.f <- fill.NA(gene.r,mode="mean")
#gene.f <- fill.NA(gene.r,mode="knn")
#gene.f <- fill.NA(gene.r,mode="wknn")

tmp <- filter.std(gene.f,min.std=1) 

gene.s <- standardise(tmp)


c <- 9

m <- mestimate(gene.s) 
cl <- mfuzz(gene.s, c = c, m = m) 

m <- mestimate(gene.s) 
cluster.range <- 2:15
md  <- Dmin(gene.s,m=m,crange= cluster.range,repeats=3,visu=TRUE)

df <- data.frame(Cluster = cluster.range, Dmin = md)

ggplot(df, aes(x = Cluster, y = Dmin)) +
  geom_line(color = "black") +
  geom_point(color = "red", size = 3) +
  labs(title = NULL, x = "Cluster number", y = "Min. centroid distance") +
  theme_bw() + theme(panel.grid = element_blank())

export::graph2pdf(file = paste(result_dir, 'Dmin_clusters.pdf', sep = "/"),height = 4,width = 6)
export::graph2png(file = paste(result_dir, 'Dmin_clusters.png', sep = "/"),height = 4,width = 6)


library(vegan)
vegan <- cascadeKM(gene_tpm, inf.gr = 2, sup.gr = 10, iter = 1000)
vegan$results

plot(vegan, sortg=TRUE,grpmts.plot = TRUE)

library(ClusterGVis)


getClusters(exp = gene_tpm)

cl$size 
cl$cluster[cl$cluster == 1] 

## cluster cores
# membership values can also indicate the similarity of vectors to each other.
gene.s

cl.thres <- acore(gene.s,cl,min.acore=0.5)  ## a posteriori
head(cl.thres[[1]])
table(cl$cluster)
unlist(lapply(cl.thres, nrow))

library(RColorBrewer)
color.2 <- colorRampPalette(rev(c("#ff0000", "Yellow", "OliveDrab1")))(1000)

mfuzz.plot(gene.s,cl,mfrow=c(3,3),
           new.window= FALSE,
           time.labels= colnames(gene.s),
           colo = color.2)
export::graph2pdf(file = paste0(result_dir, 'mfuzz_clusters_GSE234998.pdf'),height = 7,width = 12)
export::graph2png(file = paste0(result_dir, 'mfuzz_clusters_GSE234998.png'),height = 7,width = 12)


genelist = names(cl$cluster[cl$cluster %in% c(2:6)])
saveRDS(genelist, file = "./Rdata/cluster_mfuzz_GSE234998.rds")
write.table(genelist, file = paste0(result_dir, "cluster_mfuzz_GSE234998.txt"), row.names = F, col.names = F, quote = F)






time_points <- 1:ncol(gene_tpm)  

cluster_centers <- cl$centers  

monotonic_clusters <- list()
for(i in 1:nrow(cluster_centers)) {
  rho <- cor(time_points, cluster_centers[i,], method="spearman")
  if(abs(rho) > 0.8) { 
    direction <- ifelse(rho > 0, "increasing", "decreasing")
    monotonic_clusters[[paste0("Cluster_", i)]] <- list(
      rho = rho,
      direction = direction
    )
  }
}

monotonic_clusters



