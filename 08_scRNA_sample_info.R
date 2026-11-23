

rm(list = ls())
gc()
library(stringr)
library(dplyr)

library(GEOquery)
library(stringr)
library(Seurat)


#---------------------- GSE152048---------------------
setwd("./00_Rawdata/")
gset = getGEO('GSE162454', destdir=".",getGPL = F)
gset=gset[[1]]
pd1 = pData(gset)
gset = getGEO('GSE152048', destdir=".",getGPL = F)
gset=gset[[1]]
pd2 = pData(gset)

pd1 <- pd1[,c(1,8,9,41)]
pd1 <- pd1 %>% mutate(sample = paste0("Patient", 1:nrow(pd1)),patient = paste0("Patient", 1:nrow(pd1)), dataset = "GSE162454", 
                      tumor = "Osteosarcoma", tumor_type = NA, organism = organism_ch1)
pd1 <- pd1[,5:10]

pd2 <- pd2[,c(1,9, 38)]
pd2 <- pd2 %>% mutate(sample = gsub(" .*", "",title), patient = gsub(" .*", "",title), dataset = "GSE152048", 
                      tumor = "Osteosarcoma", tumor_type = `tumor type:ch1`, organism = organism_ch1)
pd2 <- pd2[,4:9] %>% distinct()

sample_info <- rbind(pd1, pd2)
write.csv(sample_info, file = "./sample_info.csv", quote = F)
saveRDS(sample_info, file = "./00_Rawdatas/sample_info.rds")
setwd("..")

#---------------------- GSE152048---------------------
dir='./00_Rawdata/GSE152048/'     
samples=list.files( dir)
samples 

library(data.table)
sceList = lapply(samples,function(pro){ 
  # pro=samples[1] 
  print(pro) 
  sce=CreateSeuratObject(counts =  Read10X(file.path(dir,pro) ) ,
                         project =  pro ,
                         min.cells = 5,
                         min.features = 300,)
  sce$sample = pro
  sce$Batch = "GSE152048"
  return(sce)
})


names(sceList)  

library(stringr)
samples
names(sceList) =  samples

sce.all <- merge(sceList[[1]], y= sceList[ -1 ] ,
                 add.cell.ids =  samples)
qs::qsave(sce.all ,file = "./00_Rawdata/sce_raw_GSE152048.qs")



#---------------------- GSE162454---------------------
dir='./00_Rawdata/GSE162454/'     
samples=list.files( dir)
samples 

library(data.table)
sceList = lapply(samples,function(pro){ 
  # pro=samples[1] 
  print(pro) 
  sce=CreateSeuratObject(counts =  Read10X(file.path(dir,pro) ) ,
                         project =  pro ,
                         min.cells = 5,
                         min.features = 300,)
  sce$sample = pro
  sce$Batch = "GSE162454"
  return(sce)
})


names(sceList)  

library(stringr)
samples
names(sceList) =  samples

sce.all <- merge(sceList[[1]], y= sceList[ -1 ] ,
                 add.cell.ids =  samples)
qs::qsave(sce.all ,file = "./00_Rawdata/sce_raw_GSE162454.qs")

