# 修复 B：层次 3 的留一队列（LOCO）无循环验证
# 循环性问题：PM_core/LM_core 的推导队列包含被打分队列 → AUC 高估
# LOCO 设计：
#   T1: LM_core 仅用芯片队列推导（GSE41258+GSE41568 meta）→ 在 RNA-seq PT（190609+225182 vs 50760）上打分
#       —— 50760 与 190609 均不在推导集，且打分只用 RNA-seq 避免芯片/RNA-seq 混合
#   T2: PM_core 仅用 GSE190609 推导 → 在 GSE225182 PT vs GSE50760 PT 上打分（两者均不在推导集）
# 敏感性：基因集限制为三个 RNA-seq 队列 PT 中均可检出（中位数>0）的基因，排除低端秩次偏移
set.seed(123)
suppressMessages({
  library(data.table); library(metafor); library(org.Hs.eg.db); library(AnnotationDbi)
})
raw <- "data/raw/crc_organotropism"
proc <- "data/processed/crc_organotropism"
dir.create(file.path(proc, "ssgsea"), showWarnings = FALSE)

ssgsea_manual <- function(mat, gsets) {
  N <- nrow(mat)
  es <- sapply(gsets, function(gs) {
    inset <- rownames(mat) %in% gs
    Ns <- sum(inset)
    if (Ns < 5) return(rep(NA_real_, ncol(mat)))
    apply(mat, 2, function(x) {
      o <- order(x, decreasing = TRUE)
      sum(cumsum(ifelse(inset[o], 1 / Ns, -1 / (N - Ns))))
    })
  })
  es <- apply(es, 2, function(v) (v - min(v)) / (max(v) - min(v)))
  t(es)
}
auc_report <- function(case, ctrl, tag) {
  n1 <- length(case); n2 <- length(ctrl)
  r <- rank(c(case, ctrl))[seq_len(n1)]
  auc <- (sum(r) - n1 * (n1 + 1) / 2) / (n1 * n2)
  p <- wilcox.test(case, ctrl)$p.value
  cat(sprintf("%s: AUC=%.3f (n=%d vs %d), Wilcoxon P=%.3g\n", tag, auc, n1, n2, p))
  list(auc = auc, p = p, n1 = n1, n2 = n2, tag = tag)
}

## ---------- 加载各队列 PT 表达（符号级） ----------
e <- fread(file.path(proc, "GSE190609_expr_RPKM.tsv.gz"))
ens <- sub("\\..*$", "", e[[1]])
emat <- log2(as.matrix(e[, -1], rownames = ens) + 0.5)
sym <- mapIds(org.Hs.eg.db, keys = rownames(emat), column = "SYMBOL",
              keytype = "ENSEMBL", multiVals = "first")
keep <- !is.na(sym); emat <- emat[keep, ]; sym <- sym[keep]
emat <- emat[!duplicated(sym), ]; rownames(emat) <- sym[!duplicated(sym)]
m190 <- fread(file.path(proc, "GSE190609_meta.tsv"))
m190_pt <- emat[, m190[organ == "primary tumor"]$title]
colnames(m190_pt) <- paste0("PT_PM190|", colnames(m190_pt))

e2 <- fread(file.path(proc, "GSE225182_expr_log2CPM.tsv.gz"))
g2 <- as.matrix(e2[, -1], rownames = e2[[1]])
m225_pt <- g2[, grep("_PT$", colnames(g2))]
colnames(m225_pt) <- paste0("PT_PM225|", colnames(m225_pt))

e5 <- fread(file.path(proc, "GSE50760_expr_FPKM.tsv.gz"))
g5 <- log2(as.matrix(e5[, -1], rownames = e5[[1]]) + 0.5)
m5 <- fread(file.path(proc, "GSE50760_meta.tsv"))
m760_pt <- g5[, m5[tissue == "primary colorectal cancer"]$gsm]
colnames(m760_pt) <- paste0("PT_LM760|", colnames(m760_pt))

## ---------- LOCO 基因集推导 ----------
# T1: LM_core 仅用芯片 meta（41258 LM + 41568 LM）
meta_de <- function(lst) {
  genes <- Reduce(intersect, lapply(lst, function(x) x$symbol[!is.na(x$symbol)]))
  rows <- lapply(genes, function(g) {
    ys <- vapply(lst, function(x) x[match(g, x$symbol)]$logFC[1], numeric(1))
    ts <- vapply(lst, function(x) x[match(g, x$symbol)]$t[1], numeric(1))
    vi <- (ys / ts)^2
    ok <- is.finite(vi) & vi > 0
    if (sum(ok) < 2) return(NULL)
    fit <- tryCatch(rma(yi = ys[ok], vi = vi[ok], method = "REML"), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    data.table(gene = g, beta = as.numeric(fit$beta), pval = fit$pval,
               sign_concord = mean(sign(ys[ok]) == sign(as.numeric(fit$beta))))
  })
  res <- rbindlist(rows)
  res[, fdr := p.adjust(pval, "BH")]
  res
}
tt258 <- fread(file.path(proc, "de", "GSE41258_LM_vs_PT.tsv.gz"))
tt568 <- fread(file.path(proc, "de", "GSE41568_LM_vs_PT.tsv.gz"))
meta_lm_chip <- meta_de(list(tt258, tt568))
lm_core_loco <- meta_lm_chip[fdr < 0.05 & beta >= 1 & sign_concord == 1]$gene
cat("LM_core_LOCO（纯芯片推导）:", length(lm_core_loco), "个\n")

# T2: PM_core 仅用 GSE190609（PM_vs_PT 上调 ∩ PM_vs_LM 上调）
tt609pm <- fread(file.path(proc, "de", "GSE190609_PM_vs_PT.tsv.gz"))
tt609pvlm <- fread(file.path(proc, "de", "GSE190609_PM_vs_LM.tsv.gz"))
pm_core_loco <- intersect(
  tt609pm[adj.P.Val < 0.05 & logFC >= 1]$symbol,
  tt609pvlm[adj.P.Val < 0.05 & logFC > 0.5]$symbol)
pm_core_loco <- pm_core_loco[!is.na(pm_core_loco)]
cat("PM_core_LOCO（仅 GSE190609 推导）:", length(pm_core_loco), "个\n")

## ---------- 打分：仅 RNA-seq 三队列 PT ----------
common <- Reduce(intersect, list(rownames(m190_pt), rownames(m225_pt), rownames(m760_pt)))
pool <- cbind(m190_pt[common, ], m225_pt[common, ], m760_pt[common, ])
cat("RNA-seq PT 合并:", nrow(pool), "基因 x", ncol(pool), "样本\n")

run_tests <- function(pool, lm_set, pm_set, suffix) {
  sc <- ssgsea_manual(pool, list(LM_core = intersect(lm_set, rownames(pool)),
                                 PM_core = intersect(pm_set, rownames(pool))))
  sdf <- data.table(sample = colnames(pool), grp = sub("\\|.*", "", colnames(pool)),
                    LM = sc["LM_core", ], PM = sc["PM_core", ])
  list(
    auc_report(sdf[grp == "PT_LM760"]$LM, sdf[grp %in% c("PT_PM190", "PT_PM225")]$LM,
               paste0("T1", suffix, ": LM_core(芯片推导) 肝患者PT vs 腹膜患者PT")),
    auc_report(sdf[grp == "PT_PM225"]$PM, sdf[grp == "PT_LM760"]$PM,
               paste0("T2", suffix, ": PM_core(仅609推导) 225182PT vs 50760PT"))
  )
}

res <- run_tests(pool, lm_core_loco, pm_core_loco, "")

# 敏感性：仅保留三个队列 PT 中均可检出（各队列中位数>0）的基因集成员
detected <- common[rowSums(pool[common, ] > 0) == ncol(pool)]
det2 <- common[apply(pool[common, ], 1, median) > 0]
cat("可检出基因（合并中位数>0）:", length(det2), "/", length(common), "\n")
res_sens <- run_tests(pool[det2, ], lm_core_loco, pm_core_loco, "-敏感性(可检出基因)")

res <- rbindlist(lapply(c(res, res_sens), as.data.frame))
print(res)
fwrite(res, file.path(proc, "ssgsea", "layer3_LOCO_auc.tsv"), sep = "\t")
fwrite(data.table(gene = lm_core_loco), file.path(proc, "meta", "liver_tropism_genes_LOCO_chip.tsv"), sep = "\t")
fwrite(data.table(gene = pm_core_loco), file.path(proc, "meta", "peritoneal_tropism_genes_LOCO_609.tsv"), sep = "\t")
