# scRNA_osteosarcoma
Integrative bulk and single-cell transcriptomic analyses identify progression-associated programs and a prognostic index in osteosarcoma



# Analysis code — "Integrative bulk and single-cell transcriptomic analyses identify progression-associated programs and a prognostic index in osteosarcoma"

## Execution order

| # | Script | Purpose | Main outputs |
|---|--------|---------|--------------|
| 01 | 01_GEO_datasets_preprocessing.R | Download & preprocess GEO bulk datasets | exprSet, phenoData |
| 02 | 02_TARGET_datasets_preprocessing.R | TARGET-OS count→FPKM→TPM, survival, clinical | exprSet_TARGET_OS, survData_TARGET_OS |
| 03 | 03_sample_inclusion_list.R | Merge per-sample metadata | sample_info_all.xlsx (Supplementary Table 1) |
| 04 | 04_DEG_GSE87624_limma.R | DEG: metastatic-lesion vs primary-lesion (limma) | DEG_metastasis-primary, volcano (Fig 1A) |
| 05 | 05_Mfuzz_GSE234998.R | Mfuzz soft clustering (pseudo-time) | cluster_mfuzz_GSE234998 (Fig 1C, S1A) |
| 06 | 06_progression_gene_intersection.R | Intersect DEG × Mfuzz modules | os_realated_genes (339 genes; Fig 1B)                        |
| 07 | 07_Metascape_enrichment.R | Metascape functional enrichment | enrichment plots (Fig 1D) |
| 08 | 08_scRNA_sample_info.R | scRNA sample metadata | sample_info.csv |
| 09 | 09_scRNA_QC_Harmony_clustering.R | QC, Harmony integration, clustering, clustree | seurat_object (Fig S1B–D) |
| 10 | 10_scRNA_annotation.R | SingleR + manual cell-type annotation | single cell atlas (Fig 1E) |
| 11 | 11_scRNA_cell_markers.R | Cell-type marker detection | cluster marker genes (Fig S1F–G) |
| 12 | 12_AUCell_scoring.R | AUCell scoring of 339-gene set | AUC scores (Fig 1F, S1I) |
| 13 | 13_AddModuleScore_scoring.R | AddModuleScore of 339-gene set | module scores (Fig 1F) |
| 14 | 14_Scissor_prognosis_cells.R | Scissor: prognosis-associated cells | Scissor labels (Fig 2A–B) |
| 15 | 15_inferCNV_malignant_cells.R | inferCNV + malignant classification | CNV scores (Fig 2C, S2E) |
| 16 | 16_CellChat_communication.R | CellChat ligand–receptor inference | cellchat results (Fig 2D, S2A–C) |
| 17 | 17_pySCENIC_change.py | CSV→loom conversion helper (python) | SCENIC input |
| 18 | 18_pySCENIC_pipeline.sh | pySCENIC GRNBoost2 → ctx → AUCell (bash) | SCENIC output |
| 19 | 19_SCENIC_regulons_visual.R | SCENIC regulon visualization (R part) | TF heatmap/network (Fig 2E, S2D) |
| 20 | 20_signature_feature_plots.R | Feature/density plots of model genes | Fig 2F–G, S2F, S3C |
| 21 | 21_feature_plot_density_helper.R | Helper function for density plots | — |
| 22 | 22_univariate_Cox.R | Univariate Cox in TARGET-OS | Unicox genes(40 genes) |
| 23 | 23_ComBat_harmonization_ML_cohorts.R | ComBat harmonization (ref = TARGET) of TARGET/GSE16091/GSE21257 | exprSet_correct |
| 24 | 24_machine_learning_101_combinations.R | 101 ML combinations, LASSO+GBM final model, bootstrap | output_C_index, model_fit_LASSO_GBM.rds, RS_LASSO_GBM.rds (Fig 3) |
| 25 | 25_model_evaluation_timeROC.R | KM + time-dependent ROC (1/3/5 y) | Fig 3D–I, RS_*.csv |
| 26 | 26_nomogram.R | Univariate/multivariate Cox, nomogram, calibration, DCA | Fig 4 |
| 27 | 27_ML_parameter_table.R | Full ML hyperparameter table | Supplementary Table 3 |
| 28–28c | 28_CIBERSORT_LM22.R (+ tool + LM22) | CIBERSORT (LM22) immune deconvolution | cibersort_result.rds (Fig 5A, 5C) |
| 29 | 29_CIBERSORTx_custom_reference.R | CIBERSORTx with custom scRNA reference | OS-specific deconvolution (Fig 5B, 5D) |
| 30–30b | 30_TIL_correlation_heatmap.R (+ helper) | Model-gene × immune-cell correlation | Fig 5C–D |
| 31 | 31_TIDE_visualization.R | TIDE immune-checkpoint prediction | Fig 6 (TIDE) |
| 32 | 32_pRRophetic_drug_prediction.R | pRRophetic IC50 (GDSC) for 139 drugs | Fig 6A–F |
| 33–33b | 33_GSVA_etoposide.R (+ helper) | GSEA (java GSEA) & GSVA on GSE200845 etoposide series | Fig 6G–H, S4 |
| 99 | 99_PDI_predict.R | **PDI prediction function for new samples** | PDI scores |

## Software versions

- R 4.2.1
- Key packages: limma 3.52.2 ; Mfuzz 2.58.0; Seurat 4.3.0 ; Scissor 2.0.0; inferCNV 1.14.2; CellChat 1.6.1; AUCell 1.20.2; clustree 0.5.1;
  rms 6.7-0; ggDCA 1.2; pRRophetic 0.2; survival 3.5.7; glmnet 4.1.8; gbm 2.1.8.1.
- Python (SCENIC): pySCENIC GRNBoost2 workflow (hg19-tss-centered-10kb-10species.mc9nr.genes_vs_motifs.rankings.feather,
  motifs-v9-nr.hgnc-m0.001-o0.0.tbl, hs_hgnc_tfs.txt).

## Random seeds

- ML pipeline : `seed <- 222222` (all models, CV, bootstrap loops);
  RF `nodesize = 5`.
- Mfuzz: `set.seed(123)`; Dmin criterion k = 2–15, repeats = 3.
- AUCell threshold exploration : `set.seed(12345)`.
- inferCNV k-means: `set.seed(123)`.
- pRRophetic: `set.seed(12345)`.
- Audit bootstrap CIs: `seed = 2024`, B = 1000.

## PDI for new samples

```r
source("99_PDI_predict.R")
pdi <- pdi_predict(expr)   # expr: gene x sample matrix, log2(TPM+1), ComBat-adjusted
group <- ifelse(pdi > 0.60885, "High", "Low")   # fixed TARGET-OS median cutoff
```
