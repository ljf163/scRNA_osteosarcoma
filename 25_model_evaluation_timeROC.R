rm(list = ls())
library(survival)
library(survminer)
# library(survivalROC)
library(timeROC)
library(patchwork)

project_dir = getwd()
result_dir = paste0(project_dir, "/results/ModelEvaluation/")
if(!dir.exists(result_dir)){
  dir.create(result_dir)
}

rs = readRDS("./Rdata/RS_LASSO_GBM.rds")

#---------------- TARGET -------------------------
rt = rs$TCGA
rt = rt[order(rt$RS, decreasing = F),]


res.cut <- surv_cutpoint(rt, 
                         time = "OS.time", 
                         event = "OS", 
                         variables = c("RS") 
)
# rt$group = ifelse(rt$RS > as.numeric(res.cut$cutpoint[1]), "High", "Low")

rt$group = ifelse(rt$RS > median(rt$RS), "High", "Low")
saveRDS(rt, file = "./Rdata/RS_TARGET.rds")
write.csv(rt, file = paste0(result_dir, "RS_TARGET.csv"), quote = F)

fit <- survfit(Surv(OS.time, OS) ~ group, data = rt)

p = ggsurvplot(fit, data = rt,pval = TRUE,
           palette = "npg",
           legend.title = "RiskScore",legend.labs = c("High", "Low"),
           conf.int = TRUE,legend=c(0.9,0.15),
           # surv.median.line = "hv",
           risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(),
           ggtheme = ggthemes::theme_few()
)
p$plot <- p$plot + labs(title = "TARGET")+theme(plot.title =element_text(size=12,hjust=0.5))
p
export::graph2pdf(file = paste0(result_dir, "TARGET_survival_curve.pdf"),width = 7,height = 6)
export::graph2png(file = paste0(result_dir, "TARGET_survival_curve.png"),width = 7,height = 6)

#AUC
timeROC <- with(rt, 
                timeROC(T = OS.time,
                        delta = OS,
                        marker = RS, 
                        cause = 1, 
                        times = c(1,3,5), 
                        ROC = TRUE, 
                        iid = TRUE, 
                        weighting = "marginal")) 
print(timeROC)
df <- data.frame(FPR = as.numeric(timeROC$FP),
                 TPR = as.numeric(timeROC$TP),
                 time = rep(as.factor(c(1,3,5)), each = nrow(timeROC$TP)))
mytheme <- theme(axis.title = element_text(size = 13),
                 axis.text = element_text(size = 11),
                 plot.title = element_text(size = 14,
                                           hjust = 0.5,
                                           face = "bold"),
                 panel.grid = element_blank(),
                 legend.title = element_text(size = 13),
                 legend.text = element_text(size = 13),
                 legend.background = element_blank(),
                 legend.position = c(1, 0),
                 legend.justification = c(1, 0))
ggplot() +
  geom_line(data = df,aes(x = FPR, y = TPR, color = time), size = 1) +
  geom_line(aes(x = c(0,1),y = c(0,1)),color = "grey") +
  scale_color_manual(name = NULL, values = c("#edae11","#0f6657","#c74732"),
                     labels = c(paste0("1 Years(AUC = ", round(timeROC$AUC[1],3), ")") ,
                                paste0("3 Years(AUC = ", round(timeROC$AUC[2],3), ")"),
                                paste0("5 Years(AUC = ", round(timeROC$AUC[3],3), ")")))+
  theme_bw() +
  mytheme +
  labs(x = "1 - Specificity (FPR)",
       y = "Sensitivity (TPR)") + ggtitle("TARGET")
ggsave(file =paste0(result_dir, "TARGET_ROC_curve.pdf") ,width = 7,height = 6)
ggsave(file =paste0(result_dir, "TARGET_ROC_curve.png") ,width = 7,height = 6)



#---------------- GSE16091 -------------------------
rt = rs$GSE16091
rt = rt[order(rt$RS, decreasing = F),]


res.cut <- surv_cutpoint(rt, 
                         time = "OS.time", 
                         event = "OS", 
                         variables = c("RS") 
)
rt$group = ifelse(rt$RS > as.numeric(res.cut$cutpoint[1]), "High", "Low")
write.csv(rt, file = paste0(result_dir, "RS_GSE16091.csv"), quote = F)

fit <- survfit(Surv(OS.time, OS) ~ group, data = rt)

p = ggsurvplot(fit, data = rt,pval = TRUE,
               palette = "npg",
               legend.title = "RiskScore",legend.labs = c("High", "Low"),
               conf.int = TRUE,legend=c(0.9,0.15),
               # surv.median.line = "hv",
               risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(),
               ggtheme = ggthemes::theme_few()
)
p$plot <- p$plot + labs(title = "GSE16091")+theme(plot.title =element_text(size=12,hjust=0.5))
p
export::graph2pdf(file = paste0(result_dir, "GSE16091_survival_curve.pdf"),width = 7,height = 6)
export::graph2png(file = paste0(result_dir, "GSE16091_survival_curve.png"),width = 7,height = 6)

#AUC
timeROC <- with(rt, 
                timeROC(T = OS.time,
                        delta = OS,
                        marker = RS, 
                        cause = 1, 
                        times = c(1,3,5), 
                        ROC = TRUE, 
                        iid = TRUE, 
                        weighting = "marginal")) 
print(timeROC)
df <- data.frame(FPR = as.numeric(timeROC$FP),
                 TPR = as.numeric(timeROC$TP),
                 time = rep(as.factor(c(1,3,5)), each = nrow(timeROC$TP)))
mytheme <- theme(axis.title = element_text(size = 13),
                 axis.text = element_text(size = 11),
                 plot.title = element_text(size = 14,
                                           hjust = 0.5,
                                           face = "bold"),
                 panel.grid = element_blank(),
                 legend.title = element_text(size = 13),
                 legend.text = element_text(size = 13),
                 legend.background = element_blank(),
                 legend.position = c(1, 0),
                 legend.justification = c(1, 0))
ggplot() +
  geom_line(data = df,aes(x = FPR, y = TPR, color = time), size = 1) +
  geom_line(aes(x = c(0,1),y = c(0,1)),color = "grey") +
  scale_color_manual(name = NULL, values = c("#edae11","#0f6657","#c74732"),
                     labels = c(paste0("1 Years(AUC = ", round(timeROC$AUC[1],3), ")") ,
                                paste0("3 Years(AUC = ", round(timeROC$AUC[2],3), ")"),
                                paste0("5 Years(AUC = ", round(timeROC$AUC[3],3), ")")))+
  theme_bw() +
  mytheme +
  labs(x = "1 - Specificity (FPR)",
       y = "Sensitivity (TPR)") + ggtitle("GSE16091")
ggsave(file = paste0(result_dir, "GSE16091_ROC_curve.pdf"),width = 7,height = 6)
ggsave(file = paste0(result_dir, "GSE16091_ROC_curve.png"),width = 7,height = 6)

#----------------   GSE21257 -------------------------
rt = rs$GSE21257
rt = rt[order(rt$RS, decreasing = F),]

res.cut <- surv_cutpoint(rt, 
                         time = "OS.time", 
                         event = "OS", 
                         variables = c("RS") 
)
rt$group = ifelse(rt$RS > as.numeric(res.cut$cutpoint[1]), "High", "Low")
write.csv(rt, file = paste0(result_dir, "RS_GSE21257.csv"), quote = F)

fit <- survfit(Surv(OS.time, OS) ~ group, data = rt)

p = ggsurvplot(fit, data = rt,pval = TRUE,
               palette = "npg",
               legend.title = "RiskScore",legend.labs = c("High", "Low"),
               conf.int = TRUE,legend=c(0.9,0.15),
               # surv.median.line = "hv",
               risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(),
               ggtheme = ggthemes::theme_few()
)
p$plot <- p$plot + labs(title = "GSE21257")+theme(plot.title =element_text(size=12,hjust=0.5))
p
export::graph2pdf(file = paste0(result_dir, "GSE21257_survival_curve.pdf"),width = 7,height = 6)
export::graph2png(file = paste0(result_dir, "GSE21257_survival_curve.png"),width = 7,height = 6)

#AUC
timeROC <- with(rt, 
                timeROC(T = OS.time,
                        delta = OS,
                        marker = RS, 
                        cause = 1, 
                        times = c(2,3,5), 
                        ROC = TRUE, 
                        iid = TRUE, 
                        weighting = "marginal")) 
print(timeROC)
df <- data.frame(FPR = as.numeric(timeROC$FP),
                 TPR = as.numeric(timeROC$TP),
                 time = rep(as.factor(c(1,3,5)), each = nrow(timeROC$TP)))
mytheme <- theme(axis.title = element_text(size = 13),
                 axis.text = element_text(size = 11),
                 plot.title = element_text(size = 14,
                                           hjust = 0.5,
                                           face = "bold"),
                 panel.grid = element_blank(),
                 legend.title = element_text(size = 13),
                 legend.text = element_text(size = 13),
                 legend.background = element_blank(),
                 legend.position = c(1, 0),
                 legend.justification = c(1, 0))
ggplot() +
  geom_line(data = df,aes(x = FPR, y = TPR, color = time), size = 1) +
  geom_line(aes(x = c(0,1),y = c(0,1)),color = "grey") +
  scale_color_manual(name = NULL, values = c("#edae11","#0f6657","#c74732"),
                     labels = c(paste0("1 Years(AUC = ", round(timeROC$AUC[1],3), ")") ,
                                paste0("3 Years(AUC = ", round(timeROC$AUC[2],3), ")"),
                                paste0("5 Years(AUC = ", round(timeROC$AUC[3],3), ")")))+
  theme_bw() +
  mytheme +
  labs(x = "1 - Specificity (FPR)",
       y = "Sensitivity (TPR)") + ggtitle("GSE21257")
ggsave(file = paste0(result_dir, "GSE21257_ROC_curve.pdf"),width = 7,height = 6)
ggsave(file = paste0(result_dir, "GSE21257_ROC_curve.png"),width = 7,height = 6)
