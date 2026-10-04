# 构建 Supplementary Tables S1–S10（CSV, UTF-8, 首行表头）→ supplementary/tables/
# 原则：只整合已落盘的结果文件，不新造数字。
#   - S1/S2 的 meta SE 由 metafor Wald 关系 z = beta/SE 代数反推（beta_re、pval_re 为全精度存储），
#     列名标注 se_wald_derived，README 中说明。
#   - S5 的 95% CI 由文件中已有的 2x2 计数以 fisher.test 精确重算（属确定性推导，非新分析）。
set.seed(123)
suppressMessages(library(data.table))

proc  <- "data/processed/crc_organotropism"
revd  <- "results/revision"
outdir <- "supplementary/tables"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

get_de <- function(f, valcol, newname) {
  d <- fread(file.path(proc, "de", f))
  d <- unique(d[, .(symbol, v = get(valcol))], by = "symbol")
  setnames(d, c("gene", newname))
  d
}
se_wald <- function(beta, p) {
  z <- qnorm(1 - p / 2)
  ifelse(p > 0 & is.finite(z) & z != 0, abs(beta) / abs(z), NA_real_)
}

## ---------- S1 / S2 ----------
build_core <- function(core_file, het_file, de_list, retain190609_file,
                       extra_retain = NULL, out_name) {
  core <- fread(file.path(proc, "meta", core_file))
  het  <- fread(file.path(proc, "meta", het_file))
  dt <- het[gene %in% core$gene]
  stopifnot(nrow(dt) == nrow(core))
  dt[, se_wald_derived := se_wald(beta_re, pval_re)]
  for (x in de_list) dt <- merge(dt, get_de(x$file, x$col, x$name), by = "gene", all.x = TRUE)
  ret <- fread(file.path(proc, "meta", retain190609_file))
  dt[, retained_after_purity_adj_GSE190609 := gene %in% ret$gene]
  if (!is.null(extra_retain)) {
    er <- fread(extra_retain$file)
    dt <- merge(dt, er[, .(gene, v = get(extra_retain$col))], by = "gene", all.x = TRUE)
    setnames(dt, "v", extra_retain$name)
  }
  setcolorder(dt, c("gene", "beta_re", "se_wald_derived", "pval_re", "fdr_re",
                    "k", "sign_concord", "tau2", "I2", "Q", "Q_p",
                    "beta_fe", "pval_fe", "fdr_fe"))
  setnames(dt, c("beta_re", "pval_re", "fdr_re", "k", "sign_concord"),
               c("meta_beta_re", "meta_pval_re", "meta_fdr_re", "k_cohorts", "sign_concordant"))
  setorder(dt, meta_pval_re)
  fwrite(dt, file.path(outdir, out_name))
  cat(out_name, "rows:", nrow(dt), "\n")
  invisible(dt)
}

s1 <- build_core(
  "peritoneal_tropism_genes.tsv", "revision_heterogeneity_per_gene_PM.tsv.gz",
  list(list(file = "GSE190609_PM_vs_PT.tsv.gz",  col = "logFC", name = "log2FC_GSE190609_PM_vs_PT"),
       list(file = "GSE225182_PM_vs_PT.tsv.gz",  col = "logFC", name = "log2FC_GSE225182_PM_vs_PT"),
       list(file = "GSE190609_PM_vs_LM.tsv.gz",  col = "logFC", name = "log2FC_GSE190609_PM_vs_LM")),
  "PM_core_purity_retained.tsv", NULL, "TableS1_PM_core_49genes.csv")

s2 <- build_core(
  "liver_tropism_genes.tsv", "revision_heterogeneity_per_gene_LM.tsv.gz",
  list(list(file = "GSE190609_LM_vs_PT.tsv.gz", col = "logFC", name = "log2FC_GSE190609_LM_vs_PT"),
       list(file = "GSE50760_LM_vs_PT.tsv.gz",  col = "logFC", name = "log2FC_GSE50760_LM_vs_PT"),
       list(file = "GSE41258_LM_vs_PT.tsv.gz",  col = "logFC", name = "log2FC_GSE41258_LM_vs_PT"),
       list(file = "GSE41568_LM_vs_PT.tsv.gz",  col = "logFC", name = "log2FC_GSE41568_LM_vs_PT"),
       list(file = "GSE190609_PM_vs_LM.tsv.gz", col = "logFC", name = "log2FC_GSE190609_PM_vs_LM")),
  "LM_core_purity_retained.tsv",
  list(file = file.path(revd, "lmcore_gse50760_purity.tsv"), col = "sig_adj",
       name = "retained_after_purity_adj_GSE50760"),
  "TableS2_LM_core_57genes.csv")

## ---------- S3 ----------
het_pm <- fread(file.path(proc, "meta/revision_heterogeneity_per_gene_PM.tsv.gz"))[, program := "PM_vs_PT"]
het_lm <- fread(file.path(proc, "meta/revision_heterogeneity_per_gene_LM.tsv.gz"))[, program := "LM_vs_PT"]
s3 <- rbind(het_pm, het_lm)
setcolorder(s3, c("program", "gene"))
fwrite(s3, file.path(outdir, "TableS3_full_meta_output.csv"))
cat("S3 rows:", nrow(s3), "\n")

## ---------- S4 ----------
main4 <- fread(file.path(proc, "mskmet/mutation_by_organ.tsv"))
sens  <- fread("mskmet/revision_ovary_excluded_sensitivity.tsv")
sens_s <- sens[level == "sample_all"]
# 主表（18 基因，IA+Ovary vs Liver）：从敏感性文件为 APC/SMAD4/RNF43 回补计数与 CI
pool <- sens_s[pm_def == "IA+Ovary", .(gene, PM_mut, PM_n, Liver_mut, Liver_n, CI_lo, CI_hi)]
s4_main <- merge(main4, pool, by = "gene", all.x = TRUE)
s4_main[, `:=`(analysis = "sample_level_PMpool_IAplusOvary_vs_Liver", fdr_within = fdr)]
# 敏感性行（IA_only，3 基因），OR/p 直接取自敏感性文件
s4_sens <- sens_s[pm_def == "IA_only",
  .(gene, PM_mut, PM_n, Liver_mut, Liver_n, CI_lo, CI_hi,
    OR_PM_vs_Liver = OR, p_PM_vs_Liver = p, fdr_within = fdr_within_def)]
s4_sens[, `:=`(analysis = "sample_level_intraabdominal_only_vs_Liver_ovary_excluded",
               liver_pct = NA_real_, lung_pct = NA_real_, intraab_pct = NA_real_,
               ovary_pct = NA_real_)]
cols4 <- c("analysis", "gene", "PM_mut", "PM_n", "Liver_mut", "Liver_n",
           "liver_pct", "lung_pct", "intraab_pct", "ovary_pct",
           "OR_PM_vs_Liver", "CI_lo", "CI_hi", "p_PM_vs_Liver", "fdr_within")
s4 <- rbind(s4_main[, ..cols4], s4_sens[, ..cols4])
setorder(s4, analysis, p_PM_vs_Liver)
fwrite(s4, file.path(outdir, "TableS4_MSKMET_samplelevel_mutation_by_organ.csv"))
cat("S4 rows:", nrow(s4), "\n")

## ---------- S5 ----------
pl <- fread(file.path(proc, "mskmet/patient_mutation_peritoneal.tsv"))
ci <- t(vapply(seq_len(nrow(pl)), function(i) {
  ft <- fisher.test(matrix(c(pl$PM_mut[i], pl$PM_n[i] - pl$PM_mut[i],
                             pl$nonPM_mut[i], pl$nonPM_n[i] - pl$nonPM_mut[i]),
                           nrow = 2, byrow = TRUE))
  c(ft$conf.int[1], ft$conf.int[2])
}, numeric(2)))
pl[, `:=`(OR_ci_lo = ci[, 1], OR_ci_hi = ci[, 2])]
fwrite(pl, file.path(outdir, "TableS5_MSKMET_patientlevel_fisher.csv"))
cat("S5 rows:", nrow(pl), "\n")

## ---------- S6–S10：直接整合复制 ----------
copy_tsv <- function(src, out) {
  d <- fread(src)
  fwrite(d, file.path(outdir, out))
  cat(out, "rows:", nrow(d), "\n")
}
copy_tsv(file.path(proc, "mskmet/revision_cox_samesubset.tsv"),
         "TableS6_cox_models_full_coefficients.csv")
copy_tsv(file.path(revd, "cms_patient_level.tsv"),
         "TableS7_CMS4_GLMM_and_duplicateCorrelation.csv")
copy_tsv(file.path(proc, "ssgsea/revision_permutation_auc.tsv"),
         "TableS8_permutation_baseline_1000random_sets.csv")
copy_tsv("mskmet/revision_apc_logistic.tsv",
         "TableS9_MSKMET_logistic_models1_2.csv")
copy_tsv(file.path(revd, "signature_overlap.tsv"),
         "TableS10_external_signature_overlap.csv")
cat("done\n")
