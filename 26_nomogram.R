# Nomogram
rm(list = ls())

project_dir = getwd()
result_dir = paste0(project_dir, "/results/Nomogram/")
if(!dir.exists(result_dir)){
  dir.create(result_dir)
}

rs = readRDS("./Rdata/RS_LASSO_GBM.rds")
rs_tcga = rs$TCGA
rs_tcga$sampleID = rownames(rs_tcga)

pd = readRDS("./Rdata/phenoData_TARGET_OS.rds")

rt = merge(rs_tcga, pd, by = "sampleID")
str(rt)

rt = rt[,1:7]
rt$Gender = as.numeric(factor(rt$Gender))
rt$Metastasis = as.numeric(factor(rt$Metastasis, levels = c("Non-metastatic","Metastatic")))

# colnames(rt)[4] = "PCDRGI"

str(rt)



library(survival)
library(survminer)
library(rms)



str(rt)

y<- Surv(time = rt$OS.time,event = rt$OS == 1)
Uni_cox<- function(x)
{ FML<- as.formula(paste0 ("y~",x))
cox<- coxph(FML,data=rt)
cox1<-summary(cox)
HR     <- cox1$coefficients[,2]
PValue <- cox1$coefficients[,5]
CI95L  <- cox1$conf.int[,3]
CI95U  <- cox1$conf.int[,4]
CI95M  <- (CI95L +CI95U)/2
Uni_cox<- data.frame('Characteristics' = x,
                     'Hazard Ratio' = HR,
                     'P Value' = PValue,
                     "CI95L" = CI95L,
                     "CI95U" = CI95U,
                     "CI95M " = CI95M )
return(Uni_cox)}         

names(rt)

variable.names<- colnames(rt)[4:ncol(rt)] 

cox2 <- lapply(variable.names, Uni_cox)
library(plyr)
UniCox <- ldply(cox2,data.frame)
UniCox
write.csv(UniCox, file = paste0(result_dir, "Unicox_result.csv"),row.names = F)

unique(UniCox$Characteristics[UniCox$P.Value < 0.05])

paste(unique(UniCox$Characteristics[UniCox$P.Value < 0.05]), collapse = "+")


# rt = readRDS("./Rdata/nomogram_rt.rds")
str(rt)

ddist <- datadist(rt)
options(datadist="ddist")

cox_n2=coxph(Surv(OS.time,OS == 1)~RS+Metastasis, data=rt)
# cox_n2=step(cox_n1,direction = "both")
cox_nSum=summary(cox_n2)


outMultiTab = data.frame()
outMultiTab = cbind(
  coef = cox_nSum$coefficients[, "coef"],
  HR = cox_nSum$conf.int[, "exp(coef)"],
  HR.95L = cox_nSum$conf.int[, "lower .95"],
  HR.95H = cox_nSum$conf.int[, "upper .95"],
  pvalue = cox_nSum$coefficients[, "Pr(>|z|)"]
)
outMultiTab = cbind(id = row.names(outMultiTab), outMultiTab)
write.csv(outMultiTab, file = paste0(result_dir, "Multicox_result.csv"))

ggforest(
  cox_n2,data = rt,
  main = "Hazard ratio",
  cpositions = c(0.02, 0.22, 0.4),
  fontsize = 0.7,
  refLabel = "reference",
  noDigits = 2
)
ggsave(paste0(result_dir,"nomo_forest.pdf"), height = 6, width = 8)



cox <- cph(Surv(OS.time,OS == 1)~RS+Metastasis,
           x=T,y=T,
           data=rt,
           surv=T);
cox
summary(cox)

# rt$time = rt$time * 365


library(regplot)
ndat = rt
ndat$Metastasis = ifelse(ndat$Metastasis == 1, "Non-metastatic","Metastatic")
coxph <- coxph(Surv(OS.time, OS) ~ RS+Metastasis, data = ndat)


regplot(coxph,
        observation=ndat[33,], 
        points=TRUE,
        plots=c("density","no plot"),
        failtime = c(1,3,5), 
        odds=F,
        droplines=T,
        leftlabel=T,
        prfail = TRUE, 
        showP = T, 
        rank="range", 
        interval="confidence",
        title="Cox regression"
) 
export::graph2pdf(file = paste0(result_dir, "Nomogram.pdf"),height = 8, width = 10)
export::graph2png(file = paste0(result_dir, "Nomogram.png"),height = 8, width = 10)



surv <- Survival(cox)
surv1 <- function(x)surv(1,lp=x)
surv2 <- function(x)surv(3,lp=x)
surv3 <- function(x)surv(5,lp=x)
nomoPlot<-nomogram(cox,
                   fun=list(surv1,surv2,surv3)
)

library(nomogramFormula)
results <- formula_rd(nomogram = nomoPlot)
nomo_score <- points_cal(formula = results$formula, rd = rt)

nomo_level<-as.vector(ifelse(nomo_score > median(nomo_score),"High","Low"))
nomoScore = cbind(id=rt$sampleID,cbind(rt[,(ncol(rt)-1):ncol(rt)],nomo_score,nomo_level))
write.csv(nomoScore,paste0(result_dir, "nomoScore.csv"),quote=F)

nomo_output = cbind(rt,nomoScore[,4:5])




ddist <- datadist(rt)
options(datadist="ddist")


nomo1 <- cph(Surv(OS.time,OS == 1)~RS+Metastasis,
             x=T,y=T,
             data=rt,
             surv=T,
             time.inc = 1)
nomo3 <- cph(Surv(OS.time,OS == 1)~RS+Metastasis,
             x=T,y=T,
             data=rt,
             surv=T,
             time.inc = 3)
nomo5 <- cph(Surv(OS.time,OS == 1)~RS+Metastasis,
             x=T,y=T,
             data=rt,
             surv=T,
             time.inc = 5)

p1<- calibrate(nomo1,
               cmethod='KM',
               method='boot',
               u=1,
               m=floor(nrow(rt)/3), 
               B=1000)

p3<- calibrate(nomo3,
               cmethod='KM',
               method='boot',
               u=3,
               m=floor(nrow(rt)/3), 
               B=1000)

p5<- calibrate(nomo5,
               cmethod='KM',
               method='boot',
               u=5,
               m=floor(nrow(rt)/3), 
               B=1000)




par(mfrow=c(1,3))
plot(p1,
     add=F,
     conf.int=T,
     subtitles = F,
     cex.subtitles=0.8, 
     lwd=2,
     lty=1,
     errbar.col="#edae11",
     xlim=c(0,1),
     ylim=c(0,1),
     xlab="Predicted OS",
     ylab="Actual OS",
     col="#edae11")
legend(x = 0.75,y= 0.2,bty = "n", legend=c("1 years"), col=c("#edae11"), lwd=2)
abline(0,1,lty=3,lwd=1,col="grey")

plot(p3,
     add=F,
     conf.int=T,
     subtitles = F,
     cex.subtitles=0.8, 
     lwd=2,
     lty=1,
     errbar.col="#0f6657",
     xlim=c(0,1),
     ylim=c(0,1),
     xlab="Predicted OS",
     ylab="Actual OS",
     col="#0f6657")
legend(x = 0.75,y= 0.2,bty = "n", legend=c("3 years"), col=c("#0f6657"), lwd=2)
abline(0,1,lty=3,lwd=1,col="grey")

plot(p5,
     add=F,
     conf.int=T,
     subtitles = F,
     cex.subtitles=0.8, 
     lwd=2,
     lty=1,
     errbar.col="#c74732",
     xlim=c(0,1),
     ylim=c(0,1),
     xlab="Predicted OS",
     ylab="Actual OS",
     col="#c74732")

legend(x = 0.75,y= 0.2,bty = "n", legend=c("5 years"), col=c("#c74732"), lwd=2)
abline(0,1,lty=3,lwd=1,col="grey")

export::graph2pdf(file = paste0(result_dir, "calibrate.pdf"), height = 5, width = 12)
export::graph2png(file = paste0(result_dir, "calibrate.png"), height = 5, width = 12)


library(survivalROC)
train = nomo_output
timeROC <- with(train, 
                timeROC(T = OS.time,
                        delta = OS,
                        marker = nomo_score, 
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
       y = "Sensitivity (TPR)")
ggsave(file = paste0(result_dir, "nomo_ROC.pdf"),width = 7,height = 6)
ggsave(file = paste0(result_dir, "nomo_ROC.png"),width = 7,height = 6)

if (!require(ggDCA)) {
  devtools::install_github("yikeshu0611/ggDCA")
}
library(rmda)
library(dcurves)
library(foreign)
library(survival)
library(dplyr)
library(ggsci)

Nomogram = coxph(Surv(OS.time,OS==1) ~RS+Metastasis,data = rt)
# pathologic_T = coxph(Surv(OS.time,OS==1) ~pathologic_T,data = rt)
# Gleason_Score = coxph(Surv(OS.time,OS==1) ~Gleason_Score,data = rt)
# pathologic_M = coxph(Surv(OS.time,OS==1) ~pathologic_M,data = rt)
# Age = coxph(Surv(OS.time,OS==1) ~Age,data = rt)
Metastasis = coxph(Surv(OS.time,OS==1) ~Metastasis,data = rt)
RS = coxph(Surv(OS.time,OS==1) ~RS,data = rt)

df2 <- ggDCA::dca(Nomogram, RS,Metastasis,times = c(1, 3, 5))
df2$timecut = ifelse(df2$time == 1, "1 years", ifelse(df2$time == 3, "3 years", "5 years"))
df2$timecut = factor(df2$timecut, levels = c("1 years", "3 years", "5 years"))

ggplot(df2,linetype = F)+
  # scale_color_jama(name="Model Type")+
  theme_bw()+labs(legend = "")+theme(panel.grid = element_blank())+
  facet_wrap(~timecut)
ggsave(paste0(result_dir, "DCA_curves.pdf"), width = 14, height = 6)
ggsave(paste0(result_dir, "DCA_curves.png"), width = 14, height = 6)


# Train
library(survival)
library(survminer)
fit <- survfit(Surv(OS.time,OS) ~ nomo_level, data = nomo_output)
ggsurvplot(fit, data = nomo_output,pval = TRUE,
           palette = "jama",
           legend.title = "Nomo Score",legend.labs = c("High", "Low"),
           conf.int = TRUE,legend=c(0.9,0.25),
           surv.median.line = "hv",
           risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(),
           ggtheme = ggthemes::theme_few()
)
export::graph2pdf(file = paste0(result_dir, "nomo_Survival.pdf"), width=8, height=6)
export::graph2png(file = paste0(result_dir, "nomo_Survival.png"), width=8, height=6)

library(forestplot)
library(stringr)

# rm(list = ls())
cox2 = read.csv("./results/Nomogram/Unicox_result.csv", header = T)

cox2$HR = round(as.numeric(cox2$Hazard.Ratio), 3)
cox2$HR.95L = round(as.numeric(cox2$CI95L), 3)
cox2$HR.95H = round(as.numeric(cox2$CI95U), 3)
cox2$pvalue = round(as.numeric(cox2$P.Value), 3)
cox2$HR2 = paste0(cox2$HR, "(", cox2$HR.95L, "-", cox2$HR.95H, ")", sep= "")
cox2 = cox2[,-c(2:6)]

cox2$Pvalue = ifelse(cox2$pvalue < 0.001, "< 0.001 ***",ifelse(cox2$pvalue < 0.01, "< 0.01 **", ifelse(cox2$pvalue < 0.05, "< 0.05 *", cox2$pvalue)))
tiltle = c("Characteristics", NA,NA,NA,NA, "HR", "pvalue")
names(tiltle) = colnames(cox2)
tiltle = t(as.data.frame(tiltle))
cox2 = rbind(tiltle, cox2)
str(cox2)
cox2$HR = as.numeric(cox2$HR)
cox2$HR.95L = as.numeric(cox2$HR.95L)
cox2$HR.95H = as.numeric(cox2$HR.95H)
str(cox2)


forestplot(
  cox2[, c(1, 6, 7)],
  mean = cox2[, 2],
  lower = cox2[, 3],
  upper = cox2[, 4],
  zero = 1,
  boxsize = 0.1,
  col = fpColors(box = '#136ad5', lines = 'black', zero = 'grey'),
  lty.ci = "solid",
  graph.pos = 3,
  is.summary=c(T,rep(F,37)),
  align = "l",
  hrzl_lines = list(
    "1" = gpar(lty=1),
    "2" = gpar(lty=1),
    "6"= gpar(lty=1)),
  colgap = unit(5, 'mm')
)
export::graph2pdf(file = paste0(result_dir, "Unicox_forest.pdf"), width = 9 ,height = 5)
export::graph2png(file = paste0(result_dir, "Unicox_forest.png"), width = 9 ,height = 5)

cox2 = read.csv("./results/Nomogram/Multicox_result.csv", header = T)

cox2$HR = round(as.numeric(cox2$HR), 3)
cox2$HR.95L = round(as.numeric(cox2$HR.95L), 3)
cox2$HR.95H = round(as.numeric(cox2$HR.95H), 3)
cox2$pvalue = round(as.numeric(cox2$pvalue), 3)
cox2$HR2 = paste0(cox2$HR, "(", cox2$HR.95L, "-", cox2$HR.95H, ")", sep= "")
cox2 = cox2[,c(-2,-3)]

cox2$Pvalue = ifelse(cox2$pvalue < 0.001, "< 0.001 ***",ifelse(cox2$pvalue < 0.01, "< 0.01 **", ifelse(cox2$pvalue < 0.05, "< 0.05 *", cox2$pvalue)))
tiltle = c("Characteristics", NA,NA,NA,NA, "HR", "pvalue")
cox2[cox2 == ""] = NA
names(tiltle) = colnames(cox2)
tiltle = t(as.data.frame(tiltle))
cox2 = rbind(tiltle, cox2)
str(cox2)
cox2$HR = as.numeric(cox2$HR)
cox2$HR.95L = as.numeric(cox2$HR.95L)
cox2$HR.95H = as.numeric(cox2$HR.95H)
str(cox2)



forestplot(
  cox2[, c(1, 6, 7)],
  mean = cox2[, 2],
  lower = cox2[, 3],
  upper = cox2[, 4],
  zero = 1,
  boxsize = 0.1,
  col = fpColors(box = 'red3', lines = 'black', zero = 'grey'),
  lty.ci = "solid",
  graph.pos = 3,
  is.summary=c(T,rep(F,37)),
  align = "l",
  hrzl_lines = list(
    "1" = gpar(lty=1),
    "2" = gpar(lty=1),
    "4"= gpar(lty=1)),
  colgap = unit(5, 'mm')
)

export::graph2pdf(file = paste0(result_dir, "Multicox_forest.pdf"), width = 9 ,height = 3)
export::graph2png(file = paste0(result_dir, "Multicox_forest.png"), width = 9 ,height = 3)
