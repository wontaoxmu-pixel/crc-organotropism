# M3 层次 3：器官趋向基因集在原发灶的回推验证（ssGSEA 打分 + ROC）
# 队列实况（2026-10-04 核实）：
#   - GSE41258：仅 5 名患者同时有原发灶+转移灶 → 无法按患者标注，PT 只能作「未选择原发灶」参考组
#   - GSE41568：PRI 与 MET 用不同编号体系（CRC### vs CRC###X/DH###），配对关系不成立 → 同上作参考组
#   - 有明确转移部位标签的原发灶：GSE190609（全部腹膜转移患者）、GSE225182（腹膜转移患者）、GSE50760（肝转移患者）
# 设计：五队列 PT 样本合并（符号交集），单次 ssGSEA（秩次法，天然抗平台效应），
#   检验 ① PM_core：腹膜患者PT vs 肝转移患者PT ② LM_core：肝转移患者PT vs 腹膜患者PT ③ 各自 vs 未选择 PT 参考
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

# ssGSEA（Barbie et al. 2009）：样本内秩次 KS 随机游走 + 基因集跨样本极差归一化
ssgsea_manual <- function(mat, gsets) {
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
  es <- apply(es, 2, function(v) (v - min(v)) / (max(v) - min(v)))
  t(es)
}

auc_report <- function(score_case, score_ctrl, tag) {
  n1 <- length(score_case); n2 <- length(score_ctrl)
  r <- rank(c(score_case, score_ctrl))[seq_len(n1)]
  auc <- (sum(r) - n1 * (n1 + 1) / 2) / (n1 * n2)
  p <- wilcox.test(score_case, score_ctrl)$p.value
  cat(sprintf("%s: AUC=%.3f (n=%d vs %d), Wilcoxon P=%.3g\n", tag, auc, n1, n2, p))
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

## ---------- 构建五队列 PT 表达矩阵（符号级） ----------
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

# d) GSE41258 PT（未选择参考组，滤掉 included=No 复扫）
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

## ---------- 合并 + ssGSEA ----------
common <- Reduce(intersect, list(rownames(m190_pt), rownames(m225_pt), rownames(m760_pt),
                                 rownames(m258_pt), rownames(m568_pt)))
cat("五队列共同基因:", length(common), "\n")
pool <- cbind(m190_pt[common, ], m225_pt[common, ], m760_pt[common, ],
              m258_pt[common, ], m568_pt[common, ])
cat("合并 PT 矩阵:", nrow(pool), "基因 x", ncol(pool), "样本\n")
sc <- ssgsea_manual(pool, gsets)
fwrite(as.data.table(sc, keep.rownames = "geneset"),
       file.path(proc, "ssgsea", "pooled_PT_ssgsea.tsv.gz"), sep = "\t")

grp <- sub("\\|.*", "", colnames(pool))
sdf <- data.table(sample = colnames(pool), cohort = grp,
                  PM_core = sc["PM_core", ], LM_core = sc["LM_core", ])
sdf[, label := fcase(cohort %in% c("GSE190609", "GSE225182"), "PT_PM",
                     cohort == "GSE50760", "PT_LM",
                     default = "PT_unselected")]
print(sdf[, .(n = .N, PM_core = round(median(PM_core), 3), LM_core = round(median(LM_core), 3)),
          by = .(label, cohort)])
fwrite(sdf, file.path(proc, "ssgsea", "pooled_PT_scores.tsv"), sep = "\t")

## ---------- ROC ----------
res <- list(
  auc_report(sdf[label == "PT_PM"]$PM_core, sdf[label == "PT_LM"]$PM_core,
             "主检验1: PM_core | 腹膜患者PT vs 肝转移患者PT"),
  auc_report(sdf[label == "PT_LM"]$LM_core, sdf[label == "PT_PM"]$LM_core,
             "主检验2: LM_core | 肝转移患者PT vs 腹膜患者PT"),
  auc_report(sdf[label == "PT_PM"]$PM_core, sdf[label == "PT_unselected"]$PM_core,
             "参考: PM_core | 腹膜患者PT vs 未选择PT"),
  auc_report(sdf[label == "PT_LM"]$LM_core, sdf[label == "PT_unselected"]$LM_core,
             "参考: LM_core | 肝转移患者PT vs 未选择PT"),
  # 阴性对照：PM_core 不应区分 肝转移患者PT vs 未选择PT
  auc_report(sdf[label == "PT_LM"]$PM_core, sdf[label == "PT_unselected"]$PM_core,
             "阴性对照: PM_core | 肝转移患者PT vs 未选择PT")
)
res <- rbindlist(lapply(res, as.data.frame))
print(res)
fwrite(res, file.path(proc, "ssgsea", "layer3_auc_summary.tsv"), sep = "\t")

## ---------- 附：GSE190609 内部探索（全部腹膜患者，LM_core 区分是否同时有肝转移） ----------
metmap <- m190[organ != "primary tumor", .(hasLM = any(organ == "liver metastasis")), by = patient]
pt <- merge(m190[organ == "primary tumor"], metmap, by = "patient")
g190 <- ssgsea_manual(emat[, pt$title], gsets)
pt$LM_core <- g190["LM_core", pt$title]
auc_report(pt$LM_core[pt$hasLM], pt$LM_core[!pt$hasLM],
           "附: GSE190609 PT 内部 LM_core PM+LM vs PM-only（探索性）")
