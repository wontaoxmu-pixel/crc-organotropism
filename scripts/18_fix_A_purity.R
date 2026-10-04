# 修复 A：纯度控制——GSE190609 配对 DE 加入上皮含量协变量
# 动机：PM_core/LM_core 与间质/肝细胞代理强相关（ρ=0.76/0.67），需量化控制纯度后器官差异剩余多少
# 设计：~0+organ+patient+epi（epi = EPCAM/KRT8/KRT19 均值 log 表达，上皮含量代理）
set.seed(123)
suppressMessages({
  library(data.table); library(limma); library(org.Hs.eg.db); library(AnnotationDbi)
})
proc <- "data/processed/crc_organotropism"

e <- fread(file.path(proc, "GSE190609_expr_RPKM.tsv.gz"))
ens <- sub("\\..*$", "", e[[1]])
emat <- log2(as.matrix(e[, -1], rownames = ens) + 0.5)
m <- fread(file.path(proc, "GSE190609_meta.tsv"))
m <- m[match(colnames(emat), m$title)]

# 上皮含量代理
sym <- mapIds(org.Hs.eg.db, keys = rownames(emat), column = "SYMBOL",
              keytype = "ENSEMBL", multiVals = "first")
epi_genes <- names(sym)[sym %in% c("EPCAM", "KRT8", "KRT19")]
epi <- colMeans(emat[epi_genes, ])

# 过滤 + 符号映射（与 11 脚本一致）
keep <- rowSums(emat >= log2(1.5)) >= 10
emat <- emat[keep, ]; sym <- sym[keep]
emat <- emat[!is.na(sym), ]; sym <- sym[!is.na(sym)]
emat <- emat[!duplicated(sym), ]; rownames(emat) <- sym[!duplicated(sym)]

run_de <- function(epi_cov = FALSE) {
  m2 <- copy(m)
  m2[, organ_f := factor(make.names(organ))]
  if (epi_cov) {
    des <- model.matrix(~ 0 + organ_f + patient + epi, data = m2)
  } else {
    des <- model.matrix(~ 0 + organ_f + patient, data = m2)
  }
  fit <- lmFit(emat, des)
  cn <- c("PM_vs_PT" = "organ_fperitoneal.metastasis - organ_fprimary.tumor",
          "LM_vs_PT" = "organ_fliver.metastasis - organ_fprimary.tumor",
          "PM_vs_LM" = "organ_fperitoneal.metastasis - organ_fliver.metastasis")
  cm <- makeContrasts(contrasts = cn, levels = des)
  fit2 <- eBayes(contrasts.fit(fit, cm), trend = TRUE)
  lapply(colnames(cm), function(cn2) {
    tt <- topTable(fit2, coef = cn2, n = Inf)
    tt$symbol <- rownames(tt)
    tt
  }) |> setNames(c("PM_vs_PT", "LM_vs_PT", "PM_vs_LM"))
}

de_raw <- run_de(FALSE)
de_adj <- run_de(TRUE)
for (cn in names(de_raw)) {
  n_raw <- sum(de_raw[[cn]]$adj.P.Val < 0.05 & abs(de_raw[[cn]]$logFC) >= 1)
  n_adj <- sum(de_adj[[cn]]$adj.P.Val < 0.05 & abs(de_adj[[cn]]$logFC) >= 1)
  cat(sprintf("%s: 原始显著 %d → 纯度校正后 %d（保留 %.0f%%）\n",
              cn, n_raw, n_adj, 100 * n_adj / max(n_raw, 1)))
}
for (cn in names(de_adj)) {
  fwrite(de_adj[[cn]], file.path(proc, "de", paste0("GSE190609_", cn, "_purityAdj.tsv.gz")), sep = "\t")
}

# 趋向基因集的纯度校正后保留率
pm_core <- fread(file.path(proc, "meta", "peritoneal_tropism_genes.tsv"))$gene
lm_core <- fread(file.path(proc, "meta", "liver_tropism_genes.tsv"))$gene
tt_pm <- as.data.table(de_adj[["PM_vs_LM"]])
sig_up <- tt_pm[adj.P.Val < 0.05 & logFC > 0]$symbol
sig_dn <- tt_pm[adj.P.Val < 0.05 & logFC < 0]$symbol
pm_keep <- intersect(pm_core, sig_up)
lm_keep <- intersect(lm_core, sig_dn)
cat(sprintf("PM_core 49 个中纯度校正后仍在 PM>LM 方向显著: %d（%.0f%%）\n",
            length(pm_keep), 100 * length(pm_keep) / 49))
cat(sprintf("LM_core 57 个中纯度校正后仍在 LM>PM 方向显著: %d（%.0f%%）\n",
            length(lm_keep), 100 * length(lm_keep) / 57))
fwrite(data.table(gene = pm_keep), file.path(proc, "meta", "PM_core_purity_retained.tsv"), sep = "\t")
fwrite(data.table(gene = lm_keep), file.path(proc, "meta", "LM_core_purity_retained.tsv"), sep = "\t")
cat("PM_core 保留基因:", paste(pm_keep, collapse = ", "), "\n")
cat("LM_core 保留基因:", paste(lm_keep, collapse = ", "), "\n")
