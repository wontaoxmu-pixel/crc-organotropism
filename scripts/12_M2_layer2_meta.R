# M2 层次2：各队列差异分析 + 随机效应 meta 整合
# 输入: data/processed/crc_organotropism/ + data/raw/crc_organotropism/
# 输出: data/processed/crc_organotropism/de/ 各队列 DE 表 + meta/ 整合结果
set.seed(123)
suppressPackageStartupMessages({
  library(data.table); library(limma); library(metafor)
})

proc <- "data/processed/crc_organotropism"
raw  <- "data/raw/crc_organotropism"
dir.create(file.path(proc, "meta"), showWarnings = FALSE, recursive = TRUE)

run_limma <- function(emat, group, block = NULL) {
  # emat: genes x samples (log 尺度); group: 两水平因子; block: 配对因子或 NULL
  design <- model.matrix(~ group)
  colnames(design) <- c("intercept", "case")
  if (is.null(block)) {
    fit <- lmFit(emat, design)
  } else {
    corfit <- duplicateCorrelation(emat, design, block = block)
    fit <- lmFit(emat, design, block = block, correlation = corfit$consensus)
  }
  fit <- eBayes(fit, trend = TRUE)
  tt <- topTable(fit, coef = "case", number = Inf, sort.by = "P")
  tt$gene <- rownames(tt)
  tt
}

summ <- function(tt, tag) {
  sig <- tt$adj.P.Val < 0.05 & abs(tt$logFC) >= 1
  cat(sprintf("%s: n_genes=%d, FDR<0.05&|logFC|>=1: %d (up %d / down %d)\n",
              tag, nrow(tt), sum(sig), sum(sig & tt$logFC > 0), sum(sig & tt$logFC < 0)))
}

## ---------- GSE225182: PM vs PT（配对 7 对）----------
m <- fread(file.path(proc, "GSE225182_meta.tsv"))
e <- fread(file.path(proc, "GSE225182_expr_log2CPM.tsv.gz"))
emat <- as.matrix(e[, -1], rownames = e[[1]])
m$xlsx_col <- make.names(m$xlsx_col)          # 与 fread 列名一致
m <- m[match(colnames(emat), m$xlsx_col), ]
paired <- m[m$subject %in% m$subject[m$tissue == "metastasized tumor into peritoneal"] &
            m$tissue %in% c("metastasized tumor into peritoneal", "Primary colorectal tumor"), ]
paired <- paired[!duplicated(paste0(paired$subject, "_", paired$tissue)), ]
emat_p <- emat[, paired$xlsx_col]
grp <- factor(ifelse(paired$tissue == "metastasized tumor into peritoneal", "PM", "PT"), levels = c("PT", "PM"))
keep <- rowSums(emat_p >= 1) >= 4   # log2CPM >= 1 至少 4 样本
tt225 <- run_limma(emat_p[keep, ], grp, block = factor(paired$subject))
summ(tt225, "GSE225182 PM_vs_PT (paired)")
tt225$symbol <- tt225$gene
fwrite(tt225, file.path(proc, "de", "GSE225182_PM_vs_PT.tsv.gz"), sep = "\t")

## ---------- GSE50760: LM vs PT（配对 18 对）----------
m <- fread(file.path(proc, "GSE50760_meta.tsv"))
e <- fread(file.path(proc, "GSE50760_expr_FPKM.tsv.gz"))
emat <- log2(as.matrix(e[, -1], rownames = e[[1]]) + 0.5)
m <- m[match(colnames(emat), m$gsm), ]
use <- m$tissue %in% c("primary colorectal cancer", "metastatic colorectal cancer to the liver")
sub <- m[use, ]
# 配对: 文件名 AMC_x.1=primary? 用 title 内 AMC 编号配对
sub$pair <- sub(".*(AMC_[0-9]+)[-.].*", "\\1", sub$title)
emat_u <- emat[, use]
grp <- factor(ifelse(sub$tissue == "metastatic colorectal cancer to the liver", "LM", "PT"), levels = c("PT", "LM"))
keep <- rowSums(emat_u >= log2(1.5)) >= 9
tt760 <- run_limma(emat_u[keep, ], grp, block = factor(sub$pair))
summ(tt760, "GSE50760 LM_vs_PT (paired)")
tt760$symbol <- tt760$gene
fwrite(tt760, file.path(proc, "de", "GSE50760_LM_vs_PT.tsv.gz"), sep = "\t")

## ---------- 芯片队列通用读取 ----------
read_chip <- function(gse) {
  path <- file.path(raw, paste0(gse, "_series_matrix.txt.gz"))
  lines <- readLines(gzfile(path))
  b <- which(lines == "!series_matrix_table_begin")
  e <- which(lines == "!series_matrix_table_end")
  dt <- fread(text = lines[(b + 1):(e - 1)])
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
  mat2 <- as.matrix(dt[, .SD, .SDcols = colnames(mat)], rownames = dt$sym)
  mat2
}

chip_de <- function(gse, meta_tissue_col = "tissue", case_label, case_name) {
  mat <- read_chip(gse)
  if (gse == "GSE41258") {
    stopifnot(requireNamespace("hgu133a.db", quietly = TRUE))
    sym <- AnnotationDbi::mapIds(hgu133a.db::hgu133a.db, keys = rownames(mat),
                                 column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
  } else {
    stopifnot(requireNamespace("hgu133plus2.db", quietly = TRUE))
    sym <- AnnotationDbi::mapIds(hgu133plus2.db::hgu133plus2.db, keys = rownames(mat),
                                 column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
  }
  gmat <- collapse_symbol(mat, sym)
  lines <- readLines(pipe(sprintf("gunzip -c '%s/%s_series_matrix.txt.gz' | grep '^!Sample_'", raw, gse)))
  getv <- function(pat) {
    l <- grep(pat, lines, value = TRUE)[1]
    v <- strsplit(l, "\t", fixed = TRUE)[[1]][-1]
    gsub('^"|"$', "", v)
  }
  char_lines <- grep("^!Sample_characteristics_ch1", lines, value = TRUE)
  getchar <- function(key) {
    l <- grep(paste0('"', key, ": "), char_lines, value = TRUE, fixed = TRUE)[1]
    v <- strsplit(l, "\t", fixed = TRUE)[[1]][-1]
    sub(paste0("^", key, ": "), "", gsub('^"|"$', "", v))
  }
  gsm <- getv("geo_accession")
  if (gse == "GSE41258") {
    tissue <- getchar("tissue")
    included <- getchar("included in analysis")
    ok_inc <- is.na(included) | included != "No"
    is_case <- tissue == case_label & ok_inc
    is_ctrl <- tissue == "Primary Tumor" & ok_inc
  } else {
    title <- getv("title")
    site  <- getchar("metastatic tumor site")
    is_case <- grepl("_MET_", title) & site == case_label
    is_ctrl <- grepl("_PRI_", title)
  }
  sel <- is_case | is_ctrl
  gsm_sel <- gsm[sel]
  gmat <- gmat[, gsm_sel]
  grp <- factor(ifelse(is_case[sel], "MET", "PT"), levels = c("PT", "MET"))
  tt <- run_limma(gmat, grp)
  summ(tt, sprintf("%s %s_vs_PT (unpaired, n=%d vs %d)", gse, case_name, sum(is_case[sel]), sum(is_ctrl[sel])))
  tt$symbol <- tt$gene
  fwrite(tt, file.path(proc, "de", sprintf("%s_%s_vs_PT.tsv.gz", gse, case_name)), sep = "\t")
  tt
}

tt258_lm   <- chip_de("GSE41258", case_label = "Liver Metastasis", case_name = "LM")
tt258_lung <- chip_de("GSE41258", case_label = "Lung Metastasis",  case_name = "LungM")
tt568_lm   <- chip_de("GSE41568", case_label = "Liver",            case_name = "LM")
tt568_lung <- chip_de("GSE41568", case_label = "Lung",             case_name = "LungM")

## ---------- 随机效应 meta ----------
tt609_pm <- fread(file.path(proc, "de", "GSE190609_PM_vs_PT.tsv.gz"))
tt609_lm <- fread(file.path(proc, "de", "GSE190609_LM_vs_PT.tsv.gz"))

meta_de <- function(lst, tag) {
  # lst: list of topTable(含 symbol, logFC, t)
  genes <- Reduce(intersect, lapply(lst, function(x) x$symbol[!is.na(x$symbol)]))
  cat(tag, "共同基因:", length(genes), "\n")
  rows <- lapply(genes, function(g) {
    ys <- vapply(lst, function(x) { r <- x[match(g, x$symbol), ]; r$logFC }, numeric(1))
    ts <- vapply(lst, function(x) { r <- x[match(g, x$symbol), ]; r$t }, numeric(1))
    vi <- (ys / ts)^2
    vi[vi <= 0 | !is.finite(vi)] <- NA
    ok <- !is.na(vi)
    if (sum(ok) < 2) return(NULL)
    fit <- tryCatch(rma(yi = ys[ok], vi = vi[ok], method = "REML"), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    b <- as.numeric(fit$beta)
    data.table(gene = g, beta = b, pval = fit$pval,
               k = sum(ok), sign_concord = mean(sign(ys[ok]) == sign(b)))
  })
  res <- rbindlist(rows)
  res[, fdr := p.adjust(pval, "BH")]
  fwrite(res, file.path(proc, "meta", paste0("meta_", tag, ".tsv.gz")), sep = "\t")
  cat(tag, "FDR<0.05 & |beta|>=1:", sum(res$fdr < 0.05 & abs(res$beta) >= 1),
      " 且方向一致:", sum(res$fdr < 0.05 & abs(res$beta) >= 1 & res$sign_concord == 1), "\n")
  res
}

meta_pm <- meta_de(list(tt609_pm, tt225), "PM_vs_PT_RNAseq")
meta_lm <- meta_de(list(tt609_lm, tt760, tt258_lm, tt568_lm), "LM_vs_PT_mixed")

## ---------- 器官趋向基因集 ----------
tt609_pvlm <- fread(file.path(proc, "de", "GSE190609_PM_vs_LM.tsv.gz"))
pm_core <- meta_pm[fdr < 0.05 & beta >= 1 & sign_concord == 1]
pm_core <- pm_core[gene %in% tt609_pvlm[adj.P.Val < 0.05 & logFC > 0.5]$symbol]
lm_core <- meta_lm[fdr < 0.05 & beta >= 1 & sign_concord == 1]
lm_core <- lm_core[gene %in% tt609_pvlm[adj.P.Val < 0.05 & logFC < -0.5]$symbol]
cat("\n== 器官趋向基因集 ==\n")
cat("腹膜趋向 (PM_core):", nrow(pm_core), " 肝趋向 (LM_core):", nrow(lm_core), "\n")
fwrite(pm_core, file.path(proc, "meta", "peritoneal_tropism_genes.tsv"), sep = "\t")
fwrite(lm_core, file.path(proc, "meta", "liver_tropism_genes.tsv"), sep = "\t")
cat("PM_core top20:", paste(head(pm_core[order(fdr)]$gene, 20), collapse = ", "), "\n")
cat("LM_core top20:", paste(head(lm_core[order(fdr)]$gene, 20), collapse = ", "), "\n")

cat("\nsessionInfo:\n"); print(sessionInfo())
