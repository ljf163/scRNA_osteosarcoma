supplementary_table <- data.frame(
  Model = character(),
  Feature_Selection = character(),
  Modeling_Method = character(),
  Key_Parameters = character(),
  Tuning_Method = character(),
  stringsAsFactors = FALSE
)

rsf_models <- data.frame(
  Model = c("RSF", "RSF + CoxBoost", "RSF + Enet", "RSF + GBM", 
            "RSF + Lasso", "RSF + plsRcox", "RSF + Ridge", 
            "RSF + StepCox", "RSF + SuperPC", "RSF + survival-SVM"),
  Feature_Selection = "RSF (Random Survival Forest)",
  Modeling_Method = c("RSF", "CoxBoost", "Elastic Net", "GBM", 
                      "Lasso", "plsRcox", "Ridge", "StepCox", 
                      "SuperPC", "survival-SVM"),
  Key_Parameters = c(
    "ntree=1000, nodesize=5, splitrule=logrank",
    "stepno=cv.optimal, penalty=optimized",
    "alpha=0.1-0.9, lambda=cv.min",
    "n.trees=cv.best, interaction.depth=3, shrinkage=0.001",
    "alpha=1, lambda=cv.min",
    "nt=cv.optimal",
    "alpha=0, lambda=cv.min",
    "direction=both/backward/forward",
    "n.threshold=20, n.components=1, s0.perc=0.5",
    "gamma.mu=1"
  ),
  Tuning_Method = c(
    "Fixed parameters",
    "10-fold CV for stepno, optimized penalty",
    "10-fold CV for lambda, grid search for alpha",
    "10-fold CV for n.trees",
    "10-fold CV for lambda",
    "CV for number of components",
    "10-fold CV for lambda",
    "Stepwise selection",
    "10-fold CV for threshold",
    "Fixed parameters"
  )
)

individual_models <- data.frame(
  Model = c("Enet", "StepCox", "CoxBoost", "plsRcox", 
            "SuperPC", "GBM", "survival-SVM", "Ridge", "Lasso"),
  Feature_Selection = "None",
  Modeling_Method = c("Elastic Net", "StepCox", "CoxBoost", "plsRcox",
                      "SuperPC", "GBM", "survival-SVM", "Ridge", "Lasso"),
  Key_Parameters = c(
    "alpha=0.1-0.9, lambda=cv.min",
    "direction=both/backward/forward",
    "stepno=cv.optimal, penalty=optimized",
    "nt=cv.optimal",
    "n.threshold=20, n.components=1, s0.perc=0.5",
    "n.trees=cv.best, interaction.depth=3, shrinkage=0.001",
    "gamma.mu=1",
    "alpha=0, lambda=cv.min",
    "alpha=1, lambda=cv.min"
  ),
  Tuning_Method = c(
    "10-fold CV for lambda, grid search for alpha",
    "Stepwise selection",
    "10-fold CV for stepno, optimized penalty",
    "CV for number of components",
    "10-fold CV for threshold",
    "10-fold CV for n.trees",
    "Fixed parameters",
    "10-fold CV for lambda",
    "10-fold CV for lambda"
  )
)

stepcox_combinations <- data.frame(
  Model = c("StepCox + CoxBoost", "StepCox + Enet", "StepCox + GBM",
            "StepCox + Lasso", "StepCox + plsRcox", "StepCox + Ridge",
            "StepCox + RSF", "StepCox + SuperPC", "StepCox + survival-SVM"),
  Feature_Selection = "StepCox",
  Modeling_Method = c("CoxBoost", "Elastic Net", "GBM", "Lasso", 
                      "plsRcox", "Ridge", "RSF", "SuperPC", "survival-SVM"),
  Key_Parameters = c(
    "stepno=cv.optimal, penalty=optimized",
    "alpha=0.1-0.9, lambda=cv.min",
    "n.trees=cv.best, interaction.depth=3, shrinkage=0.001",
    "alpha=1, lambda=cv.min",
    "nt=cv.optimal",
    "alpha=0, lambda=cv.min",
    "ntree=1000, nodesize=5, splitrule=logrank",
    "n.threshold=20, n.components=1, s0.perc=0.5",
    "gamma.mu=1"
  ),
  Tuning_Method = c(
    "10-fold CV for stepno, optimized penalty",
    "10-fold CV for lambda, grid search for alpha",
    "10-fold CV for n.trees",
    "10-fold CV for lambda",
    "CV for number of components",
    "10-fold CV for lambda",
    "Fixed parameters",
    "10-fold CV for threshold",
    "Fixed parameters"
  )
)

coxboost_combinations <- data.frame(
  Model = c("CoxBoost + Enet", "CoxBoost + GBM", "CoxBoost + Lasso",
            "CoxBoost + plsRcox", "CoxBoost + Ridge", "CoxBoost + StepCox",
            "CoxBoost + SuperPC", "CoxBoost + survival-SVM", "CoxBoost + RSF"),
  Feature_Selection = "CoxBoost",
  Modeling_Method = c("Elastic Net", "GBM", "Lasso", "plsRcox", "Ridge",
                      "StepCox", "SuperPC", "survival-SVM", "RSF"),
  Key_Parameters = c(
    "alpha=0.1-0.9, lambda=cv.min",
    "n.trees=cv.best, interaction.depth=3, shrinkage=0.001",
    "alpha=1, lambda=cv.min",
    "nt=cv.optimal",
    "alpha=0, lambda=cv.min",
    "direction=both/backward/forward",
    "n.threshold=20, n.components=1, s0.perc=0.5",
    "gamma.mu=1",
    "ntree=1000, nodesize=5, splitrule=logrank"
  ),
  Tuning_Method = c(
    "10-fold CV for lambda, grid search for alpha",
    "10-fold CV for n.trees",
    "10-fold CV for lambda",
    "CV for number of components",
    "10-fold CV for lambda",
    "Stepwise selection",
    "10-fold CV for threshold",
    "Fixed parameters",
    "Fixed parameters"
  )
)

lasso_combinations <- data.frame(
  Model = c("Lasso + CoxBoost", "Lasso + GBM", "Lasso + plsRcox",
            "Lasso + RSF", "Lasso + StepCox", "Lasso + SuperPC",
            "Lasso + survival-SVM"),
  Feature_Selection = "Lasso",
  Modeling_Method = c("CoxBoost", "GBM", "plsRcox", "RSF", "StepCox",
                      "SuperPC", "survival-SVM"),
  Key_Parameters = c(
    "stepno=cv.optimal, penalty=optimized",
    "n.trees=cv.best, interaction.depth=3, shrinkage=0.001",
    "nt=cv.optimal",
    "ntree=1000, nodesize=5, splitrule=logrank",
    "direction=both/backward/forward",
    "n.threshold=20, n.components=1, s0.perc=0.5",
    "gamma.mu=1"
  ),
  Tuning_Method = c(
    "10-fold CV for stepno, optimized penalty",
    "10-fold CV for n.trees",
    "CV for number of components",
    "Fixed parameters",
    "Stepwise selection",
    "10-fold CV for threshold",
    "Fixed parameters"
  )
)

supplementary_table <- rbind(
  rsf_models,
  individual_models,
  stepcox_combinations,
  coxboost_combinations,
  lasso_combinations
)

supplementary_table$Notes <- ""

supplementary_table$Notes[supplementary_table$Model == "RSF"] <- "nodesize tuned, seed=222222"
supplementary_table$Notes[supplementary_table$Model %in% c("CoxBoost", "RSF + CoxBoost")] <- "start.penalty=500, maxstepno=500"
supplementary_table$Notes[supplementary_table$Model %in% c("GBM", "RSF + GBM")] <- "n.minobsinnode=10, cv.folds=10"
supplementary_table$Notes[supplementary_table$Model == "SuperPC"] <- "min.features=5, max.features=all"

write.csv(supplementary_table, 
          file = paste0(result_dir, "Supplementary_Table_MachineLearning_Parameters.csv"),
          row.names = FALSE)

head(supplementary_table)

parameter_summary <- data.frame(
  Parameter_Category = c("Tree-based", "Regularization", "Boosting", "Dimension Reduction", "Stepwise"),
  Methods = c("RSF, GBM", "Lasso, Ridge, Elastic Net", "CoxBoost, GBM", "plsRcox, SuperPC", "StepCox"),
  Key_Parameters = c("ntree, nodesize, splitrule", "alpha, lambda", "stepno, penalty, n.trees", "nt, n.components", "direction"),
  Tuning_Approach = c("Fixed or CV", "Cross-validation", "Cross-validation", "Cross-validation", "Stepwise selection")
)

project_dir = getwd()
result_dir = paste0(project_dir, "/results/MachineLearning/")
if(!dir.exists(result_dir)){
  dir.create(result_dir)
}

write.csv(supplementary_table,
          file = paste0(result_dir, "Parameter_Summary.csv"),
          row.names = FALSE)

cat("Supplementary tables have been generated and saved to:", result_dir, "\n")
cat("Total number of models evaluated:", nrow(supplementary_table), "\n")
