# 修订 A9（RA-m3）：LM_core 纯度校正在 GSE50760 的跨队列重复
# 动机：审稿人质疑 GSE190609 中 LM_core 定义与纯度校正评估同队列（57/57=100% 保留）存在循环。
#       本脚本在独立配对队列 GSE50760（18 对 PT-LM，FPKM）上重复 scripts/18_fix_A_purity.R 的
#       上皮分数校正流程（~ 0 + organ + patient + epi，epi = EPCAM/KRT8/KRT19 均值 log 表达），
#       报告 LM_core 57 基因中可检出子集的校正后保留率。
# 注：LM_core 由 k=4 队列 meta 定义（含 GSE50760），故本重复仍非完全独立验证，结论中如实说明。
# 输出: results/revision/lmcore_gse50760_purity.tsv（逐基因）
#       results/revision/lmcore_gse50760_purity_summary.tsv（汇总）
set.seed(123)
suppressMessages({
  library(data.table); library(limma)
})
proc <- "data/processed/crc_organotropism"
outdir <- "results/revision"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ---------- 读入 GSE50760（与 12 脚本一致：FPKM → log2(x+0.5)，PT+LM 36 样本 18 对）----------
m <- fread(file.path(proc, "GSE50760_meta.tsv"))
e <- fread(file.path(proc, "GSE50760_expr_FPKM.tsv.gz"))
emat <- log2(as.matrix(e[, -1], rownames = e[[1]]) + 0.5)
m <- m[match(colnames(emat), m$gsm), ]
stopifnot(!any(is.na(m$gsm)))
use <- m$tissue %in% c("primary colorectal cancer", "metastatic colorectal cancer to the liver")
sub <- m[use, ]
sub$pair <- sub(".*(AMC_[0-9]+)[-.].*", "\\1", sub$title)
sub$organ <- factor(ifelse(sub$tissue == "metastatic colorectal cancer to the liver", "LM", "PT"),
                    levels = c("PT", "LM"))
emat_u <- emat[, sub$gsm]
stopifnot(all(table(sub$pair) == 2L))  # 确认 18 对完整配对

## ---------- 过滤 + 去重（与 12 脚本对该队列的过滤阈值一致：>=9 样本表达）----------
keep <- rowSums(emat_u >= log2(1.5)) >= 9
emat_f <- emat_u[keep, ]
# 同一 symbol 多行（Excel 错位或重复）时按平均表达最高保留，与 12 脚本 collapse_symbol 一致
if (anyDuplicated(rownames(emat_f))) {
  dt <- data.table(sym = rownames(emat_f), mean_expr = rowMeans(emat_f), emat_f)
  dt <- dt[order(-mean_expr)][!duplicated(sym)]
  emat_f <- as.matrix(dt[, -c("sym", "mean_expr")], rownames = dt$sym)
}
cat(sprintf("GSE50760 过滤后基因数: %d；样本: %d（%d 对）\n",
            nrow(emat_f), ncol(emat_f), length(unique(sub$pair))))

## ---------- 上皮含量代理（EPCAM/KRT8/KRT19 均值 log 表达，与 18 脚本一致）----------
epi_genes <- intersect(c("EPCAM", "KRT8", "KRT19"), rownames(emat_f))
cat("上皮代理基因可用:", paste(epi_genes, collapse = ", "), "\n")
stopifnot(length(epi_genes) == 3)
epi <- colMeans(emat_f[epi_genes, ])

## ---------- 配对 limma：原始 vs 纯度校正（固定患者效应，与 18 脚本同构）----------
run_de <- function(epi_cov = FALSE) {
  if (epi_cov) {
    des <- model.matrix(~ 0 + organ + pair + epi, data = sub)
  } else {
    des <- model.matrix(~ 0 + organ + pair, data = sub)
  }
  fit <- lmFit(emat_f, des)
  cm <- makeContrasts(LM_vs_PT = "organLM - organPT", levels = des)
  fit2 <- eBayes(contrasts.fit(fit, cm), trend = TRUE)
  tt <- topTable(fit2, coef = "LM_vs_PT", n = Inf)
  tt$symbol <- rownames(tt)
  as.data.table(tt)
}
tt_raw <- run_de(FALSE)
tt_adj <- run_de(TRUE)

## ---------- LM_core 可检出子集的保留率 ----------
lm_core <- fread(file.path(proc, "meta", "liver_tropism_genes.tsv"))$gene
cat(sprintf("LM_core 总数: %d\n", length(lm_core)))
det <- intersect(lm_core, rownames(emat_f))
missing <- setdiff(lm_core, rownames(emat_f))
cat(sprintf("GSE50760 可检出: %d / %d；未检出: %s\n",
            length(det), length(lm_core), paste(missing, collapse = ", ")))

res <- data.table(
  gene = det,
  logFC_raw = tt_raw$logFC[match(det, tt_raw$symbol)],
  adjP_raw  = tt_raw$adj.P.Val[match(det, tt_raw$symbol)],
  logFC_adj = tt_adj$logFC[match(det, tt_adj$symbol)],
  adjP_adj  = tt_adj$adj.P.Val[match(det, tt_adj$symbol)]
)
res[, sig_raw := adjP_raw < 0.05 & logFC_raw > 0]
res[, sig_adj := adjP_adj < 0.05 & logFC_adj > 0]
fwrite(res, file.path(outdir, "lmcore_gse50760_purity.tsv"), sep = "\t")

n_det <- nrow(res)
n_raw <- sum(res$sig_raw)
n_adj <- sum(res$sig_adj)
cat(sprintf("LM_core 可检出子集（n=%d）：原始显著 LM>PT %d（%.1f%%）→ 纯度校正后 %d（%.1f%%）\n",
            n_det, n_raw, 100 * n_raw / n_det, n_adj, 100 * n_adj / n_det))
cat("校正后仍显著基因:", paste(res$gene[res$sig_adj], collapse = ", "), "\n")
cat("校正后丢失基因:", paste(res$gene[!res$sig_adj], collapse = ", "), "\n")

summary_dt <- data.table(
  metric = c("cohort", "pairs", "genes_after_filter", "lm_core_total", "lm_core_detectable",
             "lm_core_missing", "sig_raw_n", "sig_raw_pct", "sig_adj_n", "sig_adj_pct",
             "gse190609_reference_retention"),
  value = c("GSE50760", length(unique(sub$pair)), nrow(emat_f), length(lm_core), n_det,
            paste(missing, collapse = ","), n_raw, round(100 * n_raw / n_det, 1),
            n_adj, round(100 * n_adj / n_det, 1), "57/57 (100%)")
)
fwrite(summary_dt, file.path(outdir, "lmcore_gse50760_purity_summary.tsv"), sep = "\t")
