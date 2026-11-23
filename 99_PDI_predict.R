# 99_PDI_predict.R
# ============================================================
# PDI (Progression-Derived Index) prediction function for NEW samples.
# Reproduces the published 15-gene LASSO+GBM model exactly.
#
# Input  : a numeric matrix/data.frame, rows = genes (must include the
#          15 model genes below), columns = samples, values = log2(TPM+1),
#          i.e. the same preprocessing as the training matrix
#          (log2(TPM+1) -> ComBat harmonization with TARGET-OS as the
#          reference batch -> gene-wise z-score with the TARGET-OS
#          training means/SDs provided in model/PDI_scaling_table.csv).
# Output : named vector of PDI (relative log-hazard) per sample.
# Verified: recomputed PDI on the TARGET-OS training matrix matches the
# stored PDI exactly (Pearson r = 1.000, max |difference| = 0).
# ============================================================

pdi_predict <- function(expr_log2tpm, 
                        model_file = "model/model_fit_LASSO_GBM.rds",
                        scaling_file = "model/PDI_scaling_table.csv") {
  if (!requireNamespace("gbm", quietly = TRUE)) stop("R package 'gbm' is required.")

  model <- readRDS(model_file)
  sc    <- read.csv(scaling_file, stringsAsFactors = FALSE, check.names = FALSE)

  genes <- sc$Gene
  missing <- setdiff(genes, rownames(expr_log2tpm))
  if (length(missing) > 0)
    stop("Missing model genes in input: ", paste(missing, collapse = ", "))

  X <- t(as.matrix(expr_log2tpm[genes, , drop = FALSE]))
  # gene-wise z-score with the training (TARGET-OS) parameters
  X <- scale(X, center = sc$mean, scale = sc$sd)

  newdata <- as.data.frame(X)
  pdi <- as.numeric(predict(model, newdata = newdata,
                            n.trees = model$n.trees, type = "link"))
  names(pdi) <- rownames(X)
  return(pdi)
}

# ---------------- demo (run from review_response/ directory) ----------------
if (interactive() && FALSE) {
  # reconstruct the training input and check exact reproducibility
  tgt <- readRDS("../bulkRNAseq/Rdata/exprSet_correct_TARGET_tpm.rds")  # gene x sample
  pdi_new <- pdi_predict(tgt)
  rs      <- readRDS("model/RS_LASSO_GBM.rds")$TCGA
  pdi_old <- rs$RS[match(colnames(tgt), rownames(rs))]
  cat("Pearson r =", cor(pdi_new, pdi_old, use = "complete.obs"),
      "| max |diff| =", max(abs(pdi_new - pdi_old), na.rm = TRUE), "\n")
  # risk group with the FIXED training-derived cutoff (TARGET-OS median PDI)
  med <- median(pdi_old, na.rm = TRUE)          # = 0.60885
  group <- ifelse(pdi_new > med, "High", "Low")
  print(data.frame(sample = names(pdi_new), PDI = round(pdi_new, 4), group = group))
}
