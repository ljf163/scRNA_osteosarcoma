library(dplyr)
library(stringr)
library(purrr)

fs1 <- list.files("Rdata/", pattern = "phenoData", full.names = T)


pd.list <- lapply(fs1, function(x) readRDS(x))
names(pd.list) <- basename(fs1) %>% str_replace(".rds", "") %>% str_replace("phenoData_", "")



fs2 <- list.files("Rdata/", pattern = "surv", full.names = T)

surv.list <- lapply(fs2, function(x) {
  tmp <- readRDS(x)
  if(!"sampleID" %in%  colnames(tmp)){
    tmp <- tmp %>% tibble::rownames_to_column("sampleID")
  }
  tmp$dataset <- basename(x) %>% str_replace(".rds", "") %>% str_replace("survData_", "")
  return(tmp)
})



surv <- do.call(rbind, surv.list)
surv <- surv

GSE16091 <- pd.list[["GSE16091"]]
GSE16091 <- data.frame(sampleID = GSE16091$geo_accession, dataset = "GSE16091", tissue = "Human Osteosarcoma")

GSE234998 <- pd.list[["GSE234998"]]
GSE234998 <- GSE234998 %>% mutate(geo_accession = rownames(GSE234998), dataset= "GSE234998")

GSE234998 <- GSE234998[!str_detect(GSE234998$state, "control"),]


GSE87624 <- pd.list[["GSE87624"]]
GSE87624 <- GSE87624 %>% mutate(geo_accession = rownames(GSE87624), dataset= "GSE87624")
GSE87624 <- GSE87624[!str_detect(GSE87624$state, "normal"),]


GSE21257 <- pd.list[["GSE21257"]]
GSE21257 <- data.frame(sampleID = GSE21257$geo_accession, tissue = GSE21257$source_name_ch1,
                       geo_accession = GSE21257$geo_accession,dataset = "GSE21257",
                       Age = round((gsub(" .*", "",GSE21257$`age:ch1`) %>% as.numeric()) / 12, 0) ,
                       Gender = GSE21257$`gender:ch1`,
                       Metastasis = GSE21257$`group:ch1`,
                       TumorSite = GSE21257$`tumor location:ch1`
                       )


TARGET_OS <- pd.list[["TARGET_OS"]]
TARGET_OS <- TARGET_OS %>% mutate(dataset= "TARGET_OS")
colnames(TARGET_OS)

GSE200845 <- pd.list[["GSE200845"]]
GSE200845 <- GSE200845 %>% mutate(dataset = "GSE200845")


sample_info <- data.frame(sampleID = c(
  TARGET_OS$sampleID,GSE16091$sampleID, GSE21257$sampleID, GSE234998$sampleID, GSE87624$sampleID,GSE200845$sampleID
))


sample_info <- sample_info %>% 
  left_join(surv[1:3], by = "sampleID")%>%
  left_join(TARGET_OS, by = "sampleID")%>% 
  left_join(GSE16091, by = "sampleID")%>% 
  left_join(GSE21257, by = "sampleID")%>% 
  left_join(GSE234998, by = "sampleID")%>% 
  left_join(GSE87624, by = "sampleID")%>% 
  left_join(GSE200845, by = "sampleID")

openxlsx::write.xlsx(sample_info, file = "results/sample_info_all.xlsx")


