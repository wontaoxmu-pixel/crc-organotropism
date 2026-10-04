# 修订任务 A3（RA-M1a / RA-m6）：表观与 LOCO AUC 的置换基线 + bootstrap 95% CI + 队列分层
# 核心问题：GSE190609 全为腹膜转移患者、GSE50760 全为肝转移患者，部位标签与队列/平台近乎完全混杂。
#   零假设：与真实基因集同大小、同表达检出率分位结构的随机基因集，在相同 ssGSEA 评分流程下
#   即可获得相当 AUC（即信号仅反映队列/批次结构）。
# 设计：每个关键 AUC 配置抽取 1,000 个随机基因集（set.seed(123)），报告零分布均值/区间、
#   经验 P 值（单侧 greater）与 z 分数；真实 AUC 附 2,000 次 bootstrap percentile 95% CI；
#   可分层的对比按队列分层（GSE190609 / GSE225182 分别 vs GSE50760）。
# 输入复用 scripts/14_M3_layer3_ssgsea.R 与 scripts/17_fix_B_LOCO.R 的数据加载代码（逐字复制）。
set.seed(123)
suppressMessages({
  library(data.table); library(org.Hs.eg.db); library(AnnotationDbi)
  library(hgu133a.db); library(hgu133plus2.db)
})
raw <- "data/raw/crc_organotropism"
proc <- "data/processed/crc_organotropism"
outdir <- file.path(proc, "ssgsea")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
B_PERM <- 1000L
B_BOOT <- 2000L

## ---------- 与 script 14/17 一致的 ssGSEA 手工实现（用于校验） ----------
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

# 快速等价实现：ES = sum_k w_k*(N-k+1)，其中 w_k=1/Ns（命中）或 -1/(N-Ns)。
# 因此 ES 只依赖「集合基因在该样本降序秩次上的位置集合」，可一次矩阵乘法给全部随机集合打分。
precompute_W <- function(mat) {
  N <- nrow(mat)
  O <- apply(mat, 2, order, decreasing = TRUE)
  R <- apply(O, 2, order)  # R[g, s] = 基因 g 在样本 s 降序序列中的位次
  list(W = N - R + 1, N = N, Ttot = N * (N + 1) / 2)
}
score_sets_fast <- function(pw, univ, sets) {
  H <- vapply(sets, function(s) as.numeric(univ %in% s), numeric(length(univ)))
  Ns <- colSums(H)
  sumH <- t(H) %*% pw$W                       # B x S
  ES <- sumH / Ns - (pw$Ttot - sumH) / (pw$N - Ns)
  rownames(ES) <- names(sets)
  ES
}

auc_val <- function(case, ctrl) {
  n1 <- length(case); n2 <- length(ctrl)
  r <- rank(c(case, ctrl))[seq_len(n1)]
  (sum(r) - n1 * (n1 + 1) / 2) / (n1 * n2)
}
boot_ci <- function(case, ctrl, B = B_BOOT) {
  b <- replicate(B, auc_val(sample(case, replace = TRUE), sample(ctrl, replace = TRUE)))
  as.numeric(quantile(b, c(0.025, 0.975), na.rm = TRUE))
}

# 检出率分位匹配抽样：按检出率十分位分层，每层抽取与真实集合相同的个数
draw_matched_sets <- function(det, real, B, prefix) {
  real <- intersect(real, names(det))
  u <- names(det)
  qs <- unique(quantile(det, probs = seq(0, 1, length.out = 11)))
  sets <- vector("list", B)
  if (length(qs) >= 3) {
    bin <- cut(det, breaks = qs, include.lowest = TRUE)
    cnt <- table(bin[match(real, u)])
    for (i in seq_len(B)) {
      sets[[i]] <- unlist(lapply(names(cnt), function(b) {
        poolb <- u[!is.na(bin) & bin == b]
        k <- min(cnt[[b]], length(poolb))
        poolb[sample.int(length(poolb), k)]
      }), use.names = FALSE)
    }
  } else {
    for (i in seq_len(B)) sets[[i]] <- u[sample.int(length(u), length(real))]
  }
  names(sets) <- paste0(prefix, "_null", seq_len(B))
  sets
}

## ---------- 数据加载（逐字复用 script 14/17） ----------
pm_core <- fread(file.path(proc, "meta", "peritoneal_tropism_genes.tsv"))$gene
lm_core <- fread(file.path(proc, "meta", "liver_tropism_genes.tsv"))$gene
lm_core_loco <- fread(file.path(proc, "meta", "liver_tropism_genes_LOCO_chip.tsv"))$gene
pm_core_loco <- fread(file.path(proc, "meta", "peritoneal_tropism_genes_LOCO_609.tsv"))$gene

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

# a) GSE190609 PT
e <- fread(file.path(proc, "GSE190609_expr_RPKM.tsv.gz"))
ens <- sub("\\..*$", "", e[[1]])
emat <- log2(as.matrix(e[, -1], rownames = ens) + 0.5)
sym <- mapIds(org.Hs.eg.db, keys = rownames(emat), column = "SYMBOL",
              keytype = "ENSEMBL", multiVals = "first")
keep <- !is.na(sym); emat <- emat[keep, ]; sym <- sym[keep]
emat <- emat[!duplicated(sym), ]; rownames(emat) <- sym[!duplicated(sym)]
m190 <- fread(file.path(proc, "GSE190609_meta.tsv"))
m190_pt <- emat[, m190[organ == "primary tumor"]$title]
colnames(m190_pt) <- paste0("GSE190609|", colnames(m190_pt))

# b) GSE225182 PT
e2 <- fread(file.path(proc, "GSE225182_expr_log2CPM.tsv.gz"))
g2 <- as.matrix(e2[, -1], rownames = e2[[1]])
m225_pt <- g2[, grep("_PT$", colnames(g2))]
colnames(m225_pt) <- paste0("GSE225182|", colnames(m225_pt))

# c) GSE50760 PT
e5 <- fread(file.path(proc, "GSE50760_expr_FPKM.tsv.gz"))
g5 <- log2(as.matrix(e5[, -1], rownames = e5[[1]]) + 0.5)
m5 <- fread(file.path(proc, "GSE50760_meta.tsv"))
m760_pt <- g5[, m5[tissue == "primary colorectal cancer"]$gsm]
colnames(m760_pt) <- paste0("GSE50760|", colnames(m760_pt))

# d) GSE41258 PT
mat <- read_chip("GSE41258")
symb <- mapIds(hgu133a.db, keys = rownames(mat), column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
g258 <- collapse_symbol(mat, symb)
gsm258 <- chip_gsm("GSE41258")
tis258 <- chip_char("GSE41258", "tissue")
inc258 <- chip_char("GSE41258", "included in analysis")
pt258 <- gsm258[tis258 == "Primary Tumor" & inc258 != "No"]
m258_pt <- g258[, pt258]
colnames(m258_pt) <- paste0("GSE41258|", colnames(m258_pt))

# e) GSE41568 PRI
mat5 <- read_chip("GSE41568")
symb5 <- mapIds(hgu133plus2.db, keys = rownames(mat5), column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
g568 <- collapse_symbol(mat5, symb5)
gsm568 <- chip_gsm("GSE41568")
lines5t <- readLines(pipe(sprintf("gunzip -c '%s/GSE41568_series_matrix.txt.gz' | grep '^!Sample_title'", raw)))
title568 <- gsub('^"|"$', "", strsplit(lines5t, "\t", fixed = TRUE)[[1]][-1])
pri568 <- gsm568[grepl("_PRI_", title568)]
m568_pt <- g568[, pri568]
colnames(m568_pt) <- paste0("GSE41568|", colnames(m568_pt))

# 五队列合并（script 14 表观 AUC 的打分池）
common5 <- Reduce(intersect, list(rownames(m190_pt), rownames(m225_pt), rownames(m760_pt),
                                  rownames(m258_pt), rownames(m568_pt)))
pool5 <- cbind(m190_pt[common5, ], m225_pt[common5, ], m760_pt[common5, ],
               m258_pt[common5, ], m568_pt[common5, ])
grp5 <- sub("\\|.*", "", colnames(pool5))
cat("五队列 PT 池:", nrow(pool5), "基因 x", ncol(pool5), "样本\n")

# 三队列 RNA-seq 合并（script 17 LOCO 的打分池；列名前缀沿用 script 17）
m190_pt3 <- m190_pt; colnames(m190_pt3) <- sub("^GSE190609\\|", "PT_PM190|", colnames(m190_pt3))
m225_pt3 <- m225_pt; colnames(m225_pt3) <- sub("^GSE225182\\|", "PT_PM225|", colnames(m225_pt3))
m760_pt3 <- m760_pt; colnames(m760_pt3) <- sub("^GSE50760\\|", "PT_LM760|", colnames(m760_pt3))
common3 <- Reduce(intersect, list(rownames(m190_pt3), rownames(m225_pt3), rownames(m760_pt3)))
pool3 <- cbind(m190_pt3[common3, ], m225_pt3[common3, ], m760_pt3[common3, ])
grp3 <- sub("\\|.*", "", colnames(pool3))
det2 <- common3[apply(pool3[common3, ], 1, median) > 0]
cat("RNA-seq 三队列 PT 池:", nrow(pool3), "基因 x", ncol(pool3), "样本；可检出基因", length(det2), "\n")

## ---------- 校验：快速打分与原实现一致 + 复现已落盘 AUC ----------
chk <- ssgsea_manual(pool5, list(LM_core = lm_core))["LM_core", ]
pw5_chk <- precompute_W(pool5)
fast_chk <- score_sets_fast(pw5_chk, common5, list(LM_core = lm_core))["LM_core", ]
# 手工实现做过 min-max 归一（单调），比较未归一的快速 ES 需先归一
fast_chk_n <- (fast_chk - min(fast_chk)) / (max(fast_chk) - min(fast_chk))
stopifnot(max(abs(fast_chk_n - chk[colnames(pool5)])) < 1e-8)
cat("快速 ssGSEA 与 script 14 手工实现一致（max abs diff < 1e-8）\n")

stored14 <- fread(file.path(outdir, "layer3_auc_summary.tsv"))
stored17 <- fread(file.path(outdir, "layer3_LOCO_auc.tsv"))

## ---------- 配置定义与运行 ----------
run_config <- function(cfg_id, geneset_name, real_set, mat, grp,
                       case_levels, ctrl_levels, strata = list()) {
  univ <- rownames(mat)
  rs <- intersect(real_set, univ)
  det <- setNames(rowMeans(mat > 0), univ)
  pw <- precompute_W(mat)
  null_sets <- draw_matched_sets(det, rs, B_PERM, cfg_id)
  ES_null <- score_sets_fast(pw, univ, null_sets)
  es_real <- score_sets_fast(pw, univ, setNames(list(rs), geneset_name))[1, ]

  contrasts <- c(list(list(name = "main", case = which(grp %in% case_levels),
                           ctrl = which(grp %in% ctrl_levels))), strata)
  rows <- lapply(contrasts, function(ct) {
    auc_real <- auc_val(es_real[ct$case], es_real[ct$ctrl])
    null_aucs <- apply(ES_null, 1, function(v) auc_val(v[ct$case], v[ct$ctrl]))
    ci <- boot_ci(es_real[ct$case], es_real[ct$ctrl])
    data.table(
      config_id = cfg_id, contrast = ct$name,
      case = paste(case_levels, collapse = "+"), ctrl = paste(ctrl_levels, collapse = "+"),
      geneset = geneset_name, set_size = length(rs), n_universe = length(univ),
      det_median_real = round(median(det[rs]), 4),
      det_median_universe = round(median(det[univ]), 4),
      n_case = length(ct$case), n_ctrl = length(ct$ctrl),
      auc_real = auc_real, boot_ci_lo = ci[1], boot_ci_hi = ci[2],
      null_mean = mean(null_aucs), null_sd = sd(null_aucs),
      null_q025 = as.numeric(quantile(null_aucs, 0.025)),
      null_q975 = as.numeric(quantile(null_aucs, 0.975)),
      null_min = min(null_aucs), null_max = max(null_aucs),
      emp_p_greater = (1 + sum(null_aucs >= auc_real)) / (B_PERM + 1),
      z_score = (auc_real - mean(null_aucs)) / sd(null_aucs)
    )
  })
  null_dt <- rbindlist(lapply(seq_along(contrasts), function(j) {
    ct <- contrasts[[j]]
    data.table(config_id = cfg_id, contrast = ct$name, iter = seq_len(B_PERM),
               auc = apply(ES_null, 1, function(v) auc_val(v[ct$case], v[ct$ctrl])))
  }))
  list(rows = rbindlist(rows), nulls = null_dt)
}

# 分层对比（三队列内 GSE190609、GSE225182 分别 vs GSE50760；方向随各配置主对比）
strata5_LM <- list(
  list(name = "strata_50760_vs_190609", case = which(grp5 == "GSE50760"), ctrl = which(grp5 == "GSE190609")),
  list(name = "strata_50760_vs_225182", case = which(grp5 == "GSE50760"), ctrl = which(grp5 == "GSE225182"))
)
strata5_PM <- list(
  list(name = "strata_190609_vs_50760", case = which(grp5 == "GSE190609"), ctrl = which(grp5 == "GSE50760")),
  list(name = "strata_225182_vs_50760", case = which(grp5 == "GSE225182"), ctrl = which(grp5 == "GSE50760"))
)
strata3 <- list(
  list(name = "strata_LM760_vs_PM190", case = which(grp3 == "PT_LM760"), ctrl = which(grp3 == "PT_PM190")),
  list(name = "strata_LM760_vs_PM225", case = which(grp3 == "PT_LM760"), ctrl = which(grp3 == "PT_PM225"))
)

res <- list(); nulls <- list()
cat("\n== A_apparent_LM_core: 表观 LM_core，50760PT(肝患者) vs 190609+225182PT(腹膜患者) ==\n")
r <- run_config("A_apparent_LMcore", "LM_core", lm_core, pool5, grp5,
                "GSE50760", c("GSE190609", "GSE225182"), strata5_LM)
res <- c(res, list(r$rows)); nulls <- c(nulls, list(r$nulls))

cat("== B_apparent_PM_core: 表观 PM_core，腹膜患者PT vs 50760PT ==\n")
r <- run_config("B_apparent_PMcore", "PM_core", pm_core, pool5, grp5,
                c("GSE190609", "GSE225182"), "GSE50760", strata5_PM)
res <- c(res, list(r$rows)); nulls <- c(nulls, list(r$nulls))

cat("== C_LOCO_main: LM_core(芯片推导)，50760 vs 190609+225182（RNA-seq 池） ==\n")
r <- run_config("C_LOCO_main_LMcore", "LM_core_LOCO", lm_core_loco, pool3, grp3,
                "PT_LM760", c("PT_PM190", "PT_PM225"), strata3)
res <- c(res, list(r$rows)); nulls <- c(nulls, list(r$nulls))

cat("== D_LOCO_sens: 同 C，但打分池限可检出基因 ==\n")
r <- run_config("D_LOCO_sens_LMcore", "LM_core_LOCO", lm_core_loco, pool3[det2, ], grp3,
                "PT_LM760", c("PT_PM190", "PT_PM225"), strata3)
res <- c(res, list(r$rows)); nulls <- c(nulls, list(r$nulls))

cat("== E_LOCO_T2: PM_core(仅609推导)，225182PT vs 50760PT ==\n")
r <- run_config("E_LOCO_T2_PMcore", "PM_core_LOCO", pm_core_loco, pool3, grp3,
                "PT_PM225", "PT_LM760")
res <- c(res, list(r$rows)); nulls <- c(nulls, list(r$nulls))

res <- rbindlist(res)
nulls <- rbindlist(nulls)

## ---------- 与已落盘 AUC 核对 ----------
expected <- c(A_apparent_LMcore = stored14$auc[2], B_apparent_PMcore = stored14$auc[1],
              C_LOCO_main_LMcore = stored17$auc[1], D_LOCO_sens_LMcore = stored17$auc[3],
              E_LOCO_T2_PMcore = stored17$auc[2])
got <- res[contrast == "main", setNames(auc_real, config_id)]
for (nm in names(expected)) {
  ok <- abs(got[[nm]] - expected[[nm]]) < 1e-9
  cat(sprintf("核对 %s: 复算 %.6f vs 已落盘 %.6f -> %s\n", nm, got[[nm]], expected[[nm]],
              ifelse(ok, "PASS", "FAIL")))
}

## ---------- 输出 ----------
fwrite(res, file.path(outdir, "revision_permutation_auc.tsv"), sep = "\t")
fwrite(nulls, file.path(outdir, "revision_permutation_null_distributions.tsv.gz"), sep = "\t")

cat("\n===== 主结果 =====\n")
print(res[, .(config_id, contrast, geneset, set_size, n_case, n_ctrl,
              auc_real = round(auc_real, 4), boot_ci_lo = round(boot_ci_lo, 4),
              boot_ci_hi = round(boot_ci_hi, 4), null_mean = round(null_mean, 4),
              null_q025 = round(null_q025, 4), null_q975 = round(null_q975, 4),
              emp_p_greater = signif(emp_p_greater, 4), z_score = round(z_score, 3))])
cat("\n输出：", file.path(outdir, "revision_permutation_auc.tsv"), "\n")
