# 修订 A4：ssGSEA per-sample 基准（RA-M1b / S6 / RB-m7）
# ① M3（scripts/14）的 ssgsea_manual = 样本内秩次游走 + 基因集跨样本极差归一化（es 行第 33 行）
# ② 本脚本实现严格 per-sample 版本：同样的样本内秩次游走统计量，但完全不做跨样本归一化
#    （每个样本的分数只由该样本自身的秩次决定，与其它样本无关 —— Barbie 2009 原义）
# ③ 两版逐样本分数 Spearman 相关 + 层次 3 全部 AUC 配置前后对照
# 输出：data/processed/crc_organotropism/ssgsea/revision_persample_benchmark.tsv
set.seed(123)
suppressMessages({
  library(data.table)
  library(hgu133a.db); library(hgu133plus2.db); library(org.Hs.eg.db); library(AnnotationDbi)
})
raw <- "data/raw/crc_organotropism"
proc <- "data/processed/crc_organotropism"
dir.create(file.path(proc, "ssgsea"), showWarnings = FALSE)

pm_core <- fread(file.path(proc, "meta", "peritoneal_tropism_genes.tsv"))$gene
lm_core <- fread(file.path(proc, "meta", "liver_tropism_genes.tsv"))$gene
gsets <- list(PM_core = pm_core, LM_core = lm_core)

## M3 原版：样本内秩次游走 + 跨样本极差归一化（与 scripts/14 逐行一致）
ssgsea_manual_m3 <- function(mat, gsets) {
  N <- nrow(mat)
  es <- sapply(gsets, function(gs) {
    inset <- rownames(mat) %in% gs
    Ns <- sum(inset)
    apply(mat, 2, function(x) {
      o <- order(x, decreasing = TRUE)
      hit <- inset[o]
      sum(cumsum(ifelse(hit, 1 / Ns, -1 / (N - Ns))))
    })
  })
  es <- apply(es, 2, function(v) (v - min(v)) / (max(v) - min(v)))  # 跨样本归一化
  t(es)
}

## per-sample 版：同样统计量，去掉跨样本归一化；每样本分数仅由自身秩次决定
ssgsea_persample <- function(mat, gsets) {
  N <- nrow(mat)
  es <- sapply(gsets, function(gs) {
    inset <- rownames(mat) %in% gs
    Ns <- sum(inset)
    apply(mat, 2, function(x) {
      o <- order(x, decreasing = TRUE)
      hit <- inset[o]
      sum(cumsum(ifelse(hit, 1 / Ns, -1 / (N - Ns))))
    })
  })
  t(es)  # 无任何跨样本变换
}

auc_report <- function(score_case, score_ctrl, tag) {
  n1 <- length(score_case); n2 <- length(score_ctrl)
  r <- rank(c(score_case, score_ctrl))[seq_len(n1)]
  auc <- (sum(r) - n1 * (n1 + 1) / 2) / (n1 * n2)
  p <- wilcox.test(score_case, score_ctrl)$p.value
  list(auc = auc, p = p, n1 = n1, n2 = n2, tag = tag)
}

read_chip <- function(gse) {
  lines <- readLines(gzfile(file.path(raw, paste0(gse, "_series_matrix.txt.gz"))))
  b <- which(lines == "!series_matrix_table_begin")
  en <- which(lines == "!series_matrix_table_end")
  dt <- fread(text = lines[(b + 1):(en - 1)])
  ids <- dt[[1]]
  mat <- as.matrix(as.data.frame(dt)[, -1, drop = FALSE])
  rownames(mat) <- ids
  storage.mode(mat) <- "numeric"
  if (max(mat, na.rm = TRUE) > 100) mat <- log2(mat + 1)
  mat
}
collapse_symbol <- function(mat, symbols) {
  dt <- data.table(probe = rownames(mat), sym = symbols, mat)
  dt <- dt[!is.na(sym) & sym != ""]
  dt[, mean_expr := rowMeans(.SD), .SDcols = colnames(mat)]
  dt <- dt[order(-mean_expr)][!duplicated(sym)]
  as.matrix(dt[, .SD, .SDcols = colnames(mat)], rownames = dt$sym)
}
chip_gsm <- function(gse) {
  l <- readLines(pipe(sprintf("gunzip -c '%s/%s_series_matrix.txt.gz' | grep '^!Sample_geo_accession'", raw, gse)))
  gsub('^"|"$', "", strsplit(l, "\t", fixed = TRUE)[[1]][-1])
}
chip_char <- function(gse, key) {
  lines <- readLines(pipe(sprintf("gunzip -c '%s/%s_series_matrix.txt.gz' | grep '^!Sample_characteristics_ch1'", raw, gse)))
  l <- grep(paste0('"', key, ": "), lines, value = TRUE, fixed = TRUE)[1]
  sub(paste0("^", key, ": "), "", gsub('^"|"$', "", strsplit(l, "\t", fixed = TRUE)[[1]][-1]))
}

## ---------- 五队列 PT 表达矩阵（与 scripts/14 完全相同的数据加载） ----------
# a) GSE190609 PT（腹膜转移患者）
e <- fread(file.path(proc, "GSE190609_expr_RPKM.tsv.gz"))
ens <- sub("\\..*$", "", e[[1]])
emat <- log2(as.matrix(e[, -1], rownames = ens) + 0.5)
sym <- mapIds(org.Hs.eg.db, keys = rownames(emat), column = "SYMBOL",
              keytype = "ENSEMBL", multiVals = "first")
keep <- !is.na(sym); emat <- emat[keep, ]; sym <- sym[keep]
emat <- emat[!duplicated(sym), ]; rownames(emat) <- sym[!duplicated(sym)]
m190 <- fread(file.path(proc, "GSE190609_meta.tsv"))
pt190 <- m190[organ == "primary tumor"]$title
m190_pt <- emat[, pt190]
colnames(m190_pt) <- paste0("GSE190609|", colnames(m190_pt))

# b) GSE225182 PT（腹膜转移患者）
e2 <- fread(file.path(proc, "GSE225182_expr_log2CPM.tsv.gz"))
g2 <- as.matrix(e2[, -1], rownames = e2[[1]])
pt225 <- grep("_PT$", colnames(g2), value = TRUE)
m225_pt <- g2[, pt225]
colnames(m225_pt) <- paste0("GSE225182|", colnames(m225_pt))

# c) GSE50760 PT（肝转移患者）
e5 <- fread(file.path(proc, "GSE50760_expr_FPKM.tsv.gz"))
g5 <- log2(as.matrix(e5[, -1], rownames = e5[[1]]) + 0.5)
m5 <- fread(file.path(proc, "GSE50760_meta.tsv"))
pt760 <- m5[tissue == "primary colorectal cancer"]$gsm
m760_pt <- g5[, pt760]
colnames(m760_pt) <- paste0("GSE50760|", colnames(m760_pt))

# d) GSE41258 PT（未选择参考组）
mat <- read_chip("GSE41258")
symb <- mapIds(hgu133a.db, keys = rownames(mat), column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
g258 <- collapse_symbol(mat, symb)
gsm258 <- chip_gsm("GSE41258")
tis258 <- chip_char("GSE41258", "tissue")
inc258 <- chip_char("GSE41258", "included in analysis")
pt258 <- gsm258[tis258 == "Primary Tumor" & inc258 != "No"]
m258_pt <- g258[, pt258]
colnames(m258_pt) <- paste0("GSE41258|", colnames(m258_pt))

# e) GSE41568 PRI（未选择参考组）
mat5 <- read_chip("GSE41568")
symb5 <- mapIds(hgu133plus2.db, keys = rownames(mat5), column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
g568 <- collapse_symbol(mat5, symb5)
gsm568 <- chip_gsm("GSE41568")
lines5t <- readLines(pipe(sprintf("gunzip -c '%s/GSE41568_series_matrix.txt.gz' | grep '^!Sample_title'", raw)))
title568 <- gsub('^"|"$', "", strsplit(lines5t, "\t", fixed = TRUE)[[1]][-1])
pri568 <- gsm568[grepl("_PRI_", title568)]
m568_pt <- g568[, pri568]
colnames(m568_pt) <- paste0("GSE41568|", colnames(m568_pt))

common <- Reduce(intersect, list(rownames(m190_pt), rownames(m225_pt), rownames(m760_pt),
                                 rownames(m258_pt), rownames(m568_pt)))
pool <- cbind(m190_pt[common, ], m225_pt[common, ], m760_pt[common, ],
              m258_pt[common, ], m568_pt[common, ])
cat("合并 PT 矩阵:", nrow(pool), "基因 x", ncol(pool), "样本\n")

grp <- sub("\\|.*", "", colnames(pool))
lab <- fcase(grp %in% c("GSE190609", "GSE225182"), "PT_PM",
             grp == "GSE50760", "PT_LM",
             default = "PT_unselected")

## ---------- 两版打分 ----------
sc_m3 <- ssgsea_manual_m3(pool, gsets)   # M3 版（含跨样本归一化）
sc_ps <- ssgsea_persample(pool, gsets)   # per-sample 版

# 与已落盘的 M3 输出核对一致性
old <- fread(file.path(proc, "ssgsea", "pooled_PT_ssgsea.tsv.gz"))
old_mat <- as.matrix(old[, -1], rownames = old$geneset)
diff_chk <- max(abs(old_mat[rownames(sc_m3), colnames(sc_m3)] - sc_m3))
cat(sprintf("复算 M3 版与既有 pooled_PT_ssgsea.tsv.gz 最大绝对差: %.3g\n", diff_chk))

fwrite(as.data.table(sc_ps, keep.rownames = "geneset"),
       file.path(proc, "ssgsea", "revision_persample_scores.tsv.gz"), sep = "\t")

## ---------- 逐样本分数 Spearman 相关 ----------
samp <- colnames(pool)
cor_rows <- lapply(rownames(sc_m3), function(gs) {
  rho <- cor(sc_m3[gs, samp], sc_ps[gs, samp], method = "spearman")
  data.table(analysis = "score_correlation", tag = paste0("合并PT池: ", gs, " M3归一化版 vs per-sample版"),
             n1 = length(samp), n2 = NA_integer_, auc_m3 = NA_real_, p_m3 = NA_real_,
             auc_persample = NA_real_, p_persample = NA_real_, delta_auc = NA_real_,
             spearman_rho_scores = rho)
})

## ---------- 层次 3 全部 AUC 配置（两版对照） ----------
cfg <- list(
  list(case = "PT_PM", ctrl = "PT_LM", gs = "PM_core", tag = "主检验1: PM_core | 腹膜患者PT vs 肝转移患者PT"),
  list(case = "PT_LM", ctrl = "PT_PM", gs = "LM_core", tag = "主检验2: LM_core | 肝转移患者PT vs 腹膜患者PT"),
  list(case = "PT_PM", ctrl = "PT_unselected", gs = "PM_core", tag = "参考: PM_core | 腹膜患者PT vs 未选择PT"),
  list(case = "PT_LM", ctrl = "PT_unselected", gs = "LM_core", tag = "参考: LM_core | 肝转移患者PT vs 未选择PT"),
  list(case = "PT_LM", ctrl = "PT_unselected", gs = "PM_core", tag = "阴性对照: PM_core | 肝转移患者PT vs 未选择PT")
)
auc_rows <- lapply(cfg, function(cc) {
  a_m3 <- auc_report(sc_m3[cc$gs, samp[lab == cc$case]], sc_m3[cc$gs, samp[lab == cc$ctrl]], cc$tag)
  a_ps <- auc_report(sc_ps[cc$gs, samp[lab == cc$case]], sc_ps[cc$gs, samp[lab == cc$ctrl]], cc$tag)
  cat(sprintf("%s | AUC M3=%.3f (P=%.3g) vs per-sample=%.3f (P=%.3g), Δ=%.4f\n",
              cc$tag, a_m3$auc, a_m3$p, a_ps$auc, a_ps$p, a_ps$auc - a_m3$auc))
  data.table(analysis = "auc", tag = cc$tag, n1 = a_m3$n1, n2 = a_m3$n2,
             auc_m3 = a_m3$auc, p_m3 = a_m3$p,
             auc_persample = a_ps$auc, p_persample = a_ps$p,
             delta_auc = a_ps$auc - a_m3$auc, spearman_rho_scores = NA_real_)
})

## ---------- 附：GSE190609 内部探索配置（两版都在该队列内单独打分，对齐 scripts/14） ----------
metmap <- m190[organ != "primary tumor", .(hasLM = any(organ == "liver metastasis")), by = patient]
pt <- merge(m190[organ == "primary tumor"], metmap, by = "patient")
g190_m3 <- ssgsea_manual_m3(emat[, pt$title], gsets)
g190_ps <- ssgsea_persample(emat[, pt$title], gsets)
a_m3i <- auc_report(g190_m3["LM_core", pt$title][pt$hasLM], g190_m3["LM_core", pt$title][!pt$hasLM],
                    "附: GSE190609 PT 内部 LM_core PM+LM vs PM-only（探索性）")
a_psi <- auc_report(g190_ps["LM_core", pt$title][pt$hasLM], g190_ps["LM_core", pt$title][!pt$hasLM],
                    "附: GSE190609 PT 内部 LM_core PM+LM vs PM-only（探索性）")
cat(sprintf("%s | AUC M3=%.3f (P=%.3g) vs per-sample=%.3f (P=%.3g), Δ=%.4f\n",
            a_m3i$tag, a_m3i$auc, a_m3i$p, a_psi$auc, a_psi$p, a_psi$auc - a_m3i$auc))
int_auc <- data.table(analysis = "auc", tag = a_m3i$tag, n1 = a_m3i$n1, n2 = a_m3i$n2,
                      auc_m3 = a_m3i$auc, p_m3 = a_m3i$p,
                      auc_persample = a_psi$auc, p_persample = a_psi$p,
                      delta_auc = a_psi$auc - a_m3i$auc, spearman_rho_scores = NA_real_)
int_cor <- lapply(rownames(g190_m3), function(gs) {
  data.table(analysis = "score_correlation",
             tag = paste0("GSE190609内部PT: ", gs, " M3归一化版 vs per-sample版"),
             n1 = ncol(g190_m3), n2 = NA_integer_, auc_m3 = NA_real_, p_m3 = NA_real_,
             auc_persample = NA_real_, p_persample = NA_real_, delta_auc = NA_real_,
             spearman_rho_scores = cor(g190_m3[gs, ], g190_ps[gs, ], method = "spearman"))
})

out <- rbindlist(c(cor_rows, auc_rows, list(int_auc), int_cor))
gp_file <- file.path(proc, "ssgsea", "revision_persample_benchmark_gseapy.tsv")
if (file.exists(gp_file)) {
  gp <- fread(gp_file)
  gp_auc <- gp[analysis == "auc_gseapy", .(tag, auc_gseapy, p_gseapy)]
  out <- merge(out, gp_auc, by = "tag", all.x = TRUE)
  gp_cor <- gp[analysis == "score_correlation" & grepl("R per-sample", tag),
               .(geneset = sub(".*(PM_core|LM_core).*", "\\1", tag), rho_gseapy_vs_persample = spearman_rho_scores)]
  out[, geneset := sub(".*(PM_core|LM_core).*", "\\1", tag)]
  out <- merge(out, gp_cor, by = "geneset", all.x = TRUE)
  out[!grepl("^合并PT池", tag), rho_gseapy_vs_persample := NA_real_]
  out[, geneset := NULL]
}
fwrite(out, file.path(proc, "ssgsea", "revision_persample_benchmark.tsv"), sep = "\t")
print(out)

## ---------- 导出合并矩阵供 GSEApy 第二基准使用 ----------
pool_out <- data.table(gene = rownames(pool), pool)
fwrite(pool_out, file.path(proc, "ssgsea", "revision_pool_matrix_for_gseapy.tsv.gz"), sep = "\t")
fwrite(data.table(sample = samp, cohort = grp, label = lab),
       file.path(proc, "ssgsea", "revision_pool_labels.tsv"), sep = "\t")
cat("done\n")
