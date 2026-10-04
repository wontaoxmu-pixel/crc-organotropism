# A6 / RA-M6: meta 分析异质性报告（REML vs FE, leave-one-cohort-out）
# 输入: data/processed/crc_organotropism/de/*.tsv.gz + meta/peritoneal_tropism_genes.tsv
# 输出: data/processed/crc_organotropism/meta/revision_heterogeneity.tsv (汇总)
#       data/processed/crc_organotropism/meta/revision_heterogeneity_per_gene_{PM,LM}.tsv.gz (逐基因)
#       data/processed/crc_organotropism/meta/revision_heterogeneity_notes.md (k=2 REML 局限性声明)
# 效应量与方差定义与 scripts/12_M2_layer2_meta.R 完全一致: yi = logFC, vi = (logFC/t)^2
set.seed(123)
suppressPackageStartupMessages({
  library(data.table); library(metafor)
})

proc <- "data/processed/crc_organotropism"
de_dir <- file.path(proc, "de")
meta_dir <- file.path(proc, "meta")

read_de <- function(f) fread(file.path(de_dir, f))

tt609_pm <- read_de("GSE190609_PM_vs_PT.tsv.gz")
tt609_lm <- read_de("GSE190609_LM_vs_PT.tsv.gz")
tt225    <- read_de("GSE225182_PM_vs_PT.tsv.gz")
tt760    <- read_de("GSE50760_LM_vs_PT.tsv.gz")
tt258_lm <- read_de("GSE41258_LM_vs_PT.tsv.gz")
tt568_lm <- read_de("GSE41568_LM_vs_PT.tsv.gz")

# 逐基因 meta：同一基因同时跑 REML 与 FE，收集异质性统计量
meta_het <- function(lst, cohort_names, tag) {
  genes <- Reduce(intersect, lapply(lst, function(x) x$symbol[!is.na(x$symbol)]))
  cat(tag, "共同基因:", length(genes), " k =", length(lst), "\n")
  rows <- lapply(genes, function(g) {
    ys <- vapply(lst, function(x) x[match(g, x$symbol), logFC], numeric(1))
    ts <- vapply(lst, function(x) x[match(g, x$symbol), t], numeric(1))
    vi <- (ys / ts)^2
    ok <- is.finite(vi) & vi > 0
    if (sum(ok) < 2) return(NULL)
    f_re <- tryCatch(rma(yi = ys[ok], vi = vi[ok], method = "REML"), error = function(e) NULL)
    f_fe <- tryCatch(rma(yi = ys[ok], vi = vi[ok], method = "FE"),   error = function(e) NULL)
    if (is.null(f_re) || is.null(f_fe)) return(NULL)
    data.table(gene = g, k = sum(ok),
               beta_re = as.numeric(f_re$beta), pval_re = f_re$pval,
               tau2 = f_re$tau2, I2 = f_re$I2, Q = f_re$QE, Q_p = f_re$QEp,
               beta_fe = as.numeric(f_fe$beta), pval_fe = f_fe$pval,
               sign_concord = mean(sign(ys[ok]) == sign(as.numeric(f_re$beta))))
  })
  res <- rbindlist(rows)
  res[, fdr_re := p.adjust(pval_re, "BH")]
  res[, fdr_fe := p.adjust(pval_fe, "BH")]
  setattr(res, "cohorts", cohort_names)
  res
}

summarize_meta <- function(res, tag) {
  data.table(
    analysis = tag,
    n_genes = nrow(res),
    k_cohorts = max(res$k),
    tau2_median = median(res$tau2), tau2_q1 = quantile(res$tau2, 0.25),
    tau2_q3 = quantile(res$tau2, 0.75), tau2_zero_pct = mean(res$tau2 == 0) * 100,
    I2_median = median(res$I2), I2_q1 = quantile(res$I2, 0.25),
    I2_q3 = quantile(res$I2, 0.75),
    I2_ge50_pct = mean(res$I2 >= 50) * 100, I2_ge75_pct = mean(res$I2 >= 75) * 100,
    Q_sig_pct = mean(res$Q_p < 0.05) * 100,
    n_re_fdr05_absbeta1 = sum(res$fdr_re < 0.05 & abs(res$beta_re) >= 1),
    n_fe_fdr05_absbeta1 = sum(res$fdr_fe < 0.05 & abs(res$beta_fe) >= 1),
    cor_beta_re_fe = cor(res$beta_re, res$beta_fe)
  )
}

pm_res <- meta_het(list(tt609_pm, tt225), c("GSE190609", "GSE225182"), "PM_vs_PT (k=2)")
lm_res <- meta_het(list(tt609_lm, tt760, tt258_lm, tt568_lm),
                   c("GSE190609", "GSE50760", "GSE41258", "GSE41568"), "LM_vs_PT (k=4)")

fwrite(pm_res, file.path(meta_dir, "revision_heterogeneity_per_gene_PM.tsv.gz"), sep = "\t")
fwrite(lm_res, file.path(meta_dir, "revision_heterogeneity_per_gene_LM.tsv.gz"), sep = "\t")

sum_pm <- summarize_meta(pm_res, "PM_vs_PT_meta")
sum_lm <- summarize_meta(lm_res, "LM_vs_PT_meta")

## ---------- PM leave-one-cohort-out（k=2 时即单队列分析）----------
pm_core <- fread(file.path(meta_dir, "peritoneal_tropism_genes.tsv"))
core_genes <- pm_core$gene
cat("PM_core 基因数:", length(core_genes), "\n")

loco_one <- function(tt, cohort) {
  sub <- tt[match(core_genes, tt$symbol), ]
  n_found <- sum(!is.na(sub$logFC))
  up <- sub$logFC > 0
  data.table(
    analysis = paste0("PM_core_LOCO_", cohort),
    n_genes = n_found,
    direction_up_n = sum(up, na.rm = TRUE),
    direction_up_pct = mean(up, na.rm = TRUE) * 100,
    rawP05_up_n = sum(sub$P.Value < 0.05 & up, na.rm = TRUE),
    rawP05_up_pct = mean(sub$P.Value < 0.05 & up, na.rm = TRUE) * 100,
    fdr05_up_n = sum(sub$adj.P.Val < 0.05 & up, na.rm = TRUE),
    fdr05_up_pct = mean(sub$adj.P.Val < 0.05 & up, na.rm = TRUE) * 100,
    abslogFC1_up_n = sum(sub$logFC >= 1, na.rm = TRUE),
    abslogFC1_up_pct = mean(sub$logFC >= 1, na.rm = TRUE) * 100,
    median_logFC = median(sub$logFC, na.rm = TRUE)
  )
}
loco_609 <- loco_one(tt609_pm, "GSE190609")
loco_225 <- loco_one(tt225, "GSE225182")

# 非循环校验：PM meta 全部显著上调基因（不设 sign_concord 过滤）在单队列中的方向一致率
pm_sig <- pm_res[fdr_re < 0.05 & beta_re >= 1]
noncirc <- rbindlist(lapply(list(GSE190609 = tt609_pm, GSE225182 = tt225), function(tt) {
  sub <- tt[match(pm_sig$gene, tt$symbol), ]
  data.table(n = sum(!is.na(sub$logFC)),
             dir_up_pct = mean(sub$logFC > 0, na.rm = TRUE) * 100,
             fdr05_up_pct = mean(sub$adj.P.Val < 0.05 & sub$logFC > 0, na.rm = TRUE) * 100)
}), idcol = "cohort")
cat("PM meta 显著上调基因(不设方向过滤) n =", nrow(pm_sig), "\n")
print(noncirc)

## ---------- 核心基因集自身的异质性 ----------
lm_core_genes <- fread(file.path(meta_dir, "liver_tropism_genes.tsv"))$gene
core_het <- function(res, core, tag) {
  a <- res[gene %in% core]
  data.table(
    analysis = tag, n_genes = nrow(a),
    tau2_median = median(a$tau2), tau2_max = max(a$tau2),
    I2_median = median(a$I2), I2_max = max(a$I2),
    Q_sig_n = sum(a$Q_p < 0.05),
    Q_sig_genes = paste(a[Q_p < 0.05]$gene, collapse = ",")
  )
}
core_pm <- core_het(pm_res, core_genes, "PM_core_heterogeneity")
core_lm <- core_het(lm_res, lm_core_genes, "LM_core_heterogeneity")
print(core_pm); print(core_lm)

## ---------- 汇总输出 ----------
# 长表：section / item / metric / value
long <- rbindlist(list(
  melt(sum_pm, id.vars = "analysis", variable.name = "metric", value.name = "value")[, section := "meta_heterogeneity"],
  melt(sum_lm, id.vars = "analysis", variable.name = "metric", value.name = "value")[, section := "meta_heterogeneity"],
  melt(loco_609, id.vars = "analysis", variable.name = "metric", value.name = "value")[, section := "PM_core_LOCO"],
  melt(loco_225, id.vars = "analysis", variable.name = "metric", value.name = "value")[, section := "PM_core_LOCO"],
  melt(core_pm, id.vars = "analysis", variable.name = "metric", value.name = "value")[, section := "core_gene_heterogeneity"],
  melt(core_lm, id.vars = "analysis", variable.name = "metric", value.name = "value")[, section := "core_gene_heterogeneity"],
  data.table(section = "PM_meta_sig_noncircular", analysis = paste0("PM_sig_up_", noncirc$cohort),
             metric = "n_genes", value = noncirc$n),
  data.table(section = "PM_meta_sig_noncircular", analysis = paste0("PM_sig_up_", noncirc$cohort),
             metric = "dir_up_pct", value = noncirc$dir_up_pct),
  data.table(section = "PM_meta_sig_noncircular", analysis = paste0("PM_sig_up_", noncirc$cohort),
             metric = "fdr05_up_pct", value = noncirc$fdr05_up_pct)
), fill = TRUE)
setcolorder(long, c("section", "analysis", "metric", "value"))
fwrite(long, file.path(meta_dir, "revision_heterogeneity.tsv"), sep = "\t")

print(sum_pm); print(sum_lm); print(loco_609); print(loco_225)

## ---------- k=2 REML 局限性声明（供 Methods/Discussion 引用）----------
notes <- sprintf(
'# A6 meta 异质性补充说明（%s）

## k=2 REML 局限性声明（建议措辞，供 Methods/Discussion 引用）

英文：
"For the peritoneal meta-analysis, only two RNA-seq cohorts (GSE190609 and GSE225182) were available. With k = 2 studies, between-study variance (tau2) estimated by restricted maximum likelihood (REML) is notoriously imprecise, and heterogeneity statistics (tau2, I2, Cochran\'s Q) should be interpreted as descriptive rather than inferential. We therefore (i) additionally report fixed-effect model estimates as a sensitivity analysis, (ii) performed leave-one-cohort-out analyses, which for k = 2 reduce to single-cohort re-analysis, and (iii) required sign concordance across cohorts for inclusion in the PM_core signature, so that no gene enters the signature on the strength of a single cohort."

中文对照：
"腹膜转移 meta 分析仅纳入 2 个 RNA-seq 队列（GSE190609、GSE225182）。当 k=2 时，REML 估计的研究间方差（tau2）精度很差，tau2、I2、Cochran Q 等异质性统计量只能作描述性解读。因此本研究：(i) 同时报告固定效应模型作为敏感性分析；(ii) 进行 leave-one-cohort-out 分析（k=2 时退化为单队列重分析）；(iii) PM_core 基因集要求跨队列方向一致，确保没有任何基因仅凭单一队列的证据进入签名。"

## 方法学备注
- 效应量 yi = limma logFC，抽样方差 vi = (logFC/t)^2，与 scripts/12_M2_layer2_meta.R 一致。
- PM_core 的 49 个基因在定义时已要求两个队列方向一致（sign_concord = 1），因此 LOCO 中的方向一致率存在结构性循环（恒为 100%%）；有意义的稳健性指标是单队列显著率（尤其 GSE225182 作为相对独立验证）与 |logFC| >= 1 的比例。非循环版本见 revision_heterogeneity.tsv 的 PM_meta_sig_noncircular 段。
- LM meta（k=4）混合了 RNA-seq（GSE190609、GSE50760）与芯片（GSE41258、GSE41568）队列，平台差异是异质性的预期来源之一。
', format(Sys.Date(), "%Y-%m-%d"))
writeLines(notes, file.path(meta_dir, "revision_heterogeneity_notes.md"))

cat("\nsessionInfo:\n"); print(sessionInfo())
