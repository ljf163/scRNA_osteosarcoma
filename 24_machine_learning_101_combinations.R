library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(CoxBoost)
library(survivalsvm)
library(dplyr)
library(tibble)
library(BART)
library(miscTools)
library(compareC)
library(ggplot2)
library(ggsci)
library(tidyr)
library(ggbreak)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 




#-------------
rm(list = ls())
gc()


tcga = readRDS("./Rdata/exprSet_correct_TARGET_tpm.rds")
GSE16091 = readRDS("./Rdata/exprSet_correct_GSE16091_tpm.rds")
GSE21257 = readRDS("./Rdata/exprSet_correct_GSE21257_tpm.rds")

surv_tcga = readRDS("./Rdata/survData_TARGET_OS.rds")
surv_tcga = tibble::remove_rownames(surv_tcga)
surv_tcga = tibble::column_to_rownames(surv_tcga, "sampleID")
surv_GSE16091 = readRDS("./Rdata/survData_GSE16091.rds")
surv_GSE21257 = readRDS("./Rdata/survData_GSE21257.rds")

colnames(surv_tcga) =  c("OS", "OS.time")
colnames(surv_GSE16091) =  c("OS", "OS.time")
colnames(surv_GSE21257) =  c("OS", "OS.time")


features = readRDS("./Rdata/Unicox_sig_genes.rds")
features = features[features %in% rownames(tcga) & 
                      features %in% rownames(GSE16091) &
                      features %in% rownames(GSE21257) 
]

tcga = tcga[features,] %>% t() %>% as.data.frame()
GSE16091 = GSE16091[features,] %>% t() %>% as.data.frame()
GSE21257 = GSE21257[features,] %>% t() %>% as.data.frame()


identical(rownames(tcga), rownames(surv_tcga))
identical(rownames(GSE16091), rownames(surv_GSE16091))
identical(rownames(GSE21257), rownames(surv_GSE21257))


tcga = data.frame(sample = rownames(tcga), surv_tcga, tcga)
GSE16091 = data.frame(sample = rownames(GSE16091), surv_GSE16091, GSE16091)
GSE21257 = data.frame(sample = rownames(GSE21257), surv_GSE21257, GSE21257)

tcga$OS.time = tcga$OS.time / 365
GSE16091$OS.time = GSE16091$OS.time / 365
GSE21257$OS.time = GSE21257$OS.time / 365



mm <- list(TCGA = tcga,
           GSE21257 = GSE21257, GSE16091 = GSE16091)

mm <- lapply(mm,function(x){
  x[,-c(1:3)] <- scale(x[,-c(1:3)])
  return(x)})

result <- data.frame()
est_data <- mm$TCGA
val_data_list <- mm
pre_var <- colnames(est_data)[-c(1:3)]
est_dd <- est_data[, c('OS.time', 'OS', pre_var)]
val_dd_list <- lapply(val_data_list, function(x){x[, c('OS.time', 'OS', pre_var)]})

rf_nodesize <- 5
seed <- 222222



library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
dd2 <- dd2[order(dd2$All, decreasing = T),]
dt <- dd2[, c( "TARGET","GSE16091", "GSE21257")]
rownames(dt) <- dd2$Model
cls = ggsci::pal_nejm()(3)
names(cls) = c( "TARGET","GSE16091", "GSE21257")
col_ha <- HeatmapAnnotation(which = "col", Cohort = c("TARGET","GSE16091", "GSE21257"),
                            annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
                            annotation_name_side = "left",
                            col = list(Cohort=cls),
                            annotation_legend_param = list(Cohort=list(title = "Cohort",
                                                                       title_position = "topleft",
                                                                       title_gp = gpar(fontsize = 12, fontface = "bold"),
                                                                       labels_rot = 0,
                                                                       legend_height = unit(1,"cm"),
                                                                       legend_width = unit(5,"mm"),
                                                                       labels_gp = gpar(fontsize = 9,
                                                                                        fontface = "bold"))
                            )
)
row_ha <- rowAnnotation('Mean Cindex' = anno_barplot(round(rowMeans(dt[,1:3]), 3), bar_width = 1, add_numbers = T,
                                                     labels = c("Mean Cindex"), height = unit(1, "mm"),
                                                     gp = gpar(col = "white", fill = "skyblue1"), numbers_gp = gpar(fontsize = 8),
                                                     axis_param = list(at = c(0, 0.5, 1),
                                                                       labels = c("0", "0.5", "1")),
                                                     width = unit(2.5, "cm")),
                        annotation_name_side = "bottom",
                        annotation_name_gp = gpar(fontsize = 9, fontface = "bold", angle = 90))

cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(
    round(dt[i, j], 2), 
    x, y,
    gp = gpar(
      fontsize = 8
    ))
}


heatmap <- Heatmap(dt,name = " ", 
                   heatmap_legend_param = list(title="",title_position = "topleft", labels_rot = 0,
                                               legend_height = unit(8,"cm"),
                                               legend_width = unit(5,"mm"),
                                               labels_gp = gpar(fontsize = 15, fontface = "bold")),
                   border = TRUE,
                   column_split = factor(colnames(dt), levels = c("TARGET","GSE16091", "GSE21257")),
                   column_gap = unit(3, "mm"),
                   show_column_names = F,
                   show_row_names = T,
                   col = colorRamp2(c(0.4,0.65,0.9), c("#4DBBD5B2", "white", "#E64B35B2")), 
                   column_title ="", 
                   #row_title ="Intersect Gene",
                   column_title_side = "top", 
                   row_title_side = "left",
                   row_title_rot = 90, 
                   column_title_gp = gpar(fontsize = 12, fontface = "bold",col = "black"), 
                   #row_title_gp = gpar(fontsize = 15, fontface = "bold",col = "black"),
                   cluster_columns =F,
                   cluster_rows = F,
                   column_order = c("TARGET","GSE16091", "GSE21257"),
                   show_row_dend = F, 
                   cell_fun = cell_fun,
                   top_annotation = col_ha,
                   right_annotation = row_ha
)
draw(heatmap)
dev.copy2pdf(file = paste0(result_dir, "ComplexHeatmap.pdf"), width = 10, height = 18)

# LASSO + GBM

x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit <- glmnet(x1, x2, family = "cox",alpha = 1, maxit = 1000)
cvfit <- cv.glmnet(x1, x2, family = "cox",alpha = 1, maxit = 1000)
coef <- coef(fit, s = cvfit$lambda.min)
rid = row.names(coef)[which(coef != 0)]
est_dd2 <- est_data[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('OS.time', 'OS', rid)]})
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 8)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})

saveRDS(rid, file = "./Rdata/Features_selected_LASSO.rds")
saveRDS(fit, file = "./Rdata/model_fit_LASSO_GBM.rds")
saveRDS(rs, file = "./Rdata/RS_LASSO_GBM.rds")






#Lasso
library(broom)

x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit <- glmnet(x1, x2, family = "cox",alpha = 1, maxit = 1000)
cvfit <- cv.glmnet(x1, x2, family = "cox",alpha = 1, maxit = 1000)

tidy_df <- broom::tidy(fit)
head(tidy_df)

library(ggplot2)
library(RColorBrewer)


mypalette <- c(brewer.pal(11,"BrBG"),brewer.pal(11,"Spectral"),
               brewer.pal(5,"Accent"),brewer.pal(12,"Set3"),
               brewer.pal(9,"Pastel1"),brewer.pal(8,"Pastel2"))

mypalette = CellChat::scPalette(100)

p1 <- ggplot(tidy_df, aes(log(lambda), estimate, group = term,color=term)) +
  geom_line(size=1.2)+
  geom_hline(yintercept = 0)+
  xlab("Log Lambda")+
  ylab("Coefficients")+
  scale_color_manual(name="Gene",values = mypalette)+
  ggthemes::theme_few() + guides(
    color = guide_legend(
      ncol = 2, 
      byrow = TRUE,
      reverse = T)
  )
p1
ggsave(paste0(result_dir,"lasso_lambda.pdf"), width = 8, height = 6)
ggsave(paste0(result_dir,"lasso_lambda.png"), width = 8, height = 6)

plot(cvfit)

tidy_cvdf <- broom::tidy(cvfit)
tidy_cvdf$lambda = log(tidy_cvdf$lambda)

p2 <- ggplot()+
  geom_point(data=tidy_cvdf, aes(lambda,estimate),color = ggsci::pal_lancet()(2)[2])+
  geom_errorbar(data = tidy_cvdf, aes(x=lambda,ymin=conf.low,ymax=conf.high))+
  geom_vline(xintercept =log(cvfit$lambda.min),lty="dashed")+
  geom_vline(xintercept =log(cvfit$lambda.1se),lty="dashed")+
  xlab("Log Lambda")+
  ylab("C index")+
  ggthemes::theme_few()
for (i in 1:nrow(tidy_cvdf)) {
  p2 <- p2+ annotate("text", x=tidy_cvdf$lambda[i], y=tidy_cvdf$conf.high[i] + 0.08, label= tidy_cvdf$nzero[i],size = 2)
}

p2
ggsave(paste0(result_dir,"lasso_cvfit.pdf"), width = 7, height = 6)
ggsave(paste0(result_dir,"lasso_cvfit.png"), width = 7, height = 6)



# LASSO + GBM

x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit <- glmnet(x1, x2, family = "cox",alpha = 1, maxit = 1000)
cvfit <- cv.glmnet(x1, x2, family = "cox",alpha = 1, maxit = 1000)
coef <- coef(fit, s = cvfit$lambda.min)
rid = row.names(coef)[which(coef != 0)]
est_dd2 <- est_data[,c('OS.time', 'OS', rid)]
val_dd_list2 <- lapply(val_data_list, function(x){x[,c('OS.time', 'OS', rid)]})
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10, n.cores = 8)
rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})

cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])}))%>%
  rownames_to_column('ID')




# ------------------------------
# ------------------------------

set.seed(seed)
B <- 100  
boot_result <- data.frame() 

for (b in 1:B) {
  message("Bootstrap iteration: ", b)
  
  # ==== Step 1: Bootstrap resample training data ====
  boot_idx <- sample(1:nrow(est_dd), replace = TRUE)
  boot_data <- est_dd[boot_idx, ]
  
  # LASSO + GBM
  
  x1 <- as.matrix(boot_data[, pre_var])
  x2 <- as.matrix(Surv(boot_data$OS.time, boot_data$OS))
  set.seed(seed)
  fit <- glmnet(x1, x2, family = "cox",alpha = 1, maxit = 1000)
  cvfit <- cv.glmnet(x1, x2, family = "cox",alpha = 1, maxit = 1000)
  coef <- coef(fit, s = cvfit$lambda.min)
  rid = row.names(coef)[which(coef != 0)]
  if (length(rid) == 0) next
  est_dd2 <- est_data[,c('OS.time', 'OS', rid)]
  val_dd_list2 <- lapply(val_data_list, function(x){x[,c('OS.time', 'OS', rid)]})
  set.seed(seed)
  fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
             n.trees = 10000,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 10, n.cores = 6)
  # find index for number trees with minimum CV error
  best <- which.min(fit$cv.error)
  set.seed(seed)
  fit <- gbm(formula = Surv(OS.time,OS)~., data = est_dd2, distribution = 'coxph',
             n.trees = best,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 10, n.cores = 8)
  rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, x, n.trees = best, type = 'link')))})
  
  cc <- data.frame(Cindex = sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])}))%>%
    rownames_to_column('ID')
  cc <- cc %>% mutate(boot_number = b)
  boot_result <- rbind(boot_result, cc)

}


cls = ggsci::pal_nejm()(3)
names(cls) = c( "TARGET","GSE16091", "GSE21257")

bootstrap_results <- boot_result %>% na.omit()

bootstrap_results$ID <- str_replace(bootstrap_results$ID,"TCGA", "TARGET")
bootstrap_results$ID <- factor(bootstrap_results$ID, levels = c("TARGET","GSE16091", "GSE21257"))

colnames(bootstrap_results)[1] <- "Cohort"

plot_bootstrap_distribution <- function(bootstrap_results) {
  library(ggplot2)
  library(ggpubr)
  
  p1 <- ggplot(bootstrap_results, aes(x = Cohort, y = Cindex, fill = Cohort)) +
    geom_violin(alpha = 0.7, trim = FALSE) +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
    geom_jitter(aes(colour = Cohort),width = 0.05, alpha = 0.3) +
    # geom_hline(yintercept = original_performance, 
    #            linetype = "dashed", color = "red", size = 1) +
    # geom_point(aes(y = original_performance), 
    #            color = "red", size = 3, shape = 18) +
    labs(title = "Bootstrap Performance Distribution",
         subtitle = "Violin plot + Boxplot of C-index across bootstrap samples",
         x = "Dataset", y = "C-index") +
    theme_minimal() +
    theme(legend.position = "none") +
    scale_fill_manual(values = cls) +
    scale_colour_manual(values = cls)
  
  return(p1)
}

p1 <- plot_bootstrap_distribution(bootstrap_results)

plot_bootstrap_density <- function(bootstrap_results) {
  p2 <- ggplot(bootstrap_results, aes(x = Cindex, fill = Cohort)) +
    geom_density(alpha = 0.6, adjust = 1.5) +
    # geom_vline(data = original_performance, 
    #            aes(xintercept = Cindex, color = ID),
    #            linetype = "dashed", size = 1) +
    labs(title = "Bootstrap Performance Density",
         subtitle = "Probability density of C-index across bootstrap samples",
         x = "C-index", y = "Density") +
    theme_minimal() +
    scale_fill_manual(values = cls) +
    scale_color_manual(values = cls) +
    facet_wrap(~Cohort, ncol = 1)
  
  return(p2)
}
p2 <- plot_bootstrap_density(bootstrap_results)


plot_bootstrap_qq <- function(bootstrap_results) {
  p6 <- ggplot(bootstrap_results, aes(sample = Cindex)) +
    stat_qq(alpha = 0.7) +
    stat_qq_line(color = "red", linetype = "dashed") +
    labs(title = "Bootstrap Distribution Q-Q Plot",
         subtitle = "Normality assessment of bootstrap C-index distribution",
         x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_minimal() +
    facet_wrap(~Cohort, scales = "free")
  
  return(p6)
}
p3 <- plot_bootstrap_qq(bootstrap_results)


library(patchwork)

p1
ggsave(file = paste0(result_dir, "bootstrap_distribution.pdf"), width = 8, height = 6)
ggsave(file = paste0(result_dir, "bootstrap_distribution.png"), width = 8, height = 6)

p2
ggsave(file = paste0(result_dir, "bootstrap_density.pdf"), width = 8, height = 8)
ggsave(file = paste0(result_dir, "bootstrap_density.png"), width = 8, height = 8)

p3
ggsave(file = paste0(result_dir, "bootstrap_qq.pdf"), width = 12, height = 4)
ggsave(file = paste0(result_dir, "bootstrap_qq.png"), width = 12, height = 4)


openxlsx::write.xlsx(bootstrap_results, file = paste0(result_dir, "bootstrap_results.100.xlsx"))
