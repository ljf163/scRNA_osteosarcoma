rm(list = ls())



clusters = readRDS("./Rdata/cluster_mfuzz_GSE234998.rds")


degs = readRDS("./Rdata/DEG_metastasis-primary.rds")

degs = rownames(degs)[degs$change != "NOT"]


os.genes = intersect(degs, clusters)
saveRDS(os.genes, file = "./Rdata/os_realated_genes.rds")

write.table(os.genes, file = "./results/os.PD.genes.txt", row.names = F, quote = F, col.names = F)


project_dir = getwd()
result_dir = paste0(project_dir, "/results/DEG/")
if(!dir.exists(result_dir)){
  dir.create(result_dir)
}


library(ggvenn)
xx <- list(Mfuzz_Clusters = clusters,DEGs = degs)
p1 <- ggvenn(xx,show_percentage = T,show_elements = F,label_sep = ",",
             digits = 1,stroke_color = "white",
             text_color = "black",text_size = 4,set_name_size = 5,
             fill_color = c("dodgerblue", "goldenrod1"),
             set_name_color = c("dodgerblue", "goldenrod1"))
p1
ggsave("./results/Venn.pdf", width = 7, height = 7)
ggsave("./results/Venn.png", width = 7, height = 7)
