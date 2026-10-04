# M1: CRC organotropism — GEO 队列元数据解析、表达矩阵对齐与 QC
# 输入: data/raw/crc_organotropism/  输出: data/processed/crc_organotropism/
set.seed(123)
suppressPackageStartupMessages({ library(data.table); library(readxl) })

raw <- "data/raw/crc_organotropism"
out <- "data/processed/crc_organotropism"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

parse_meta <- function(path) {
  lines <- readLines(pipe(sprintf("gunzip -c '%s' | grep '^!Sample_'", path)))
  keys <- sub("\t.*", "", lines)
  vals <- lapply(lines, function(l) {
    v <- strsplit(l, "\t", fixed = TRUE)[[1]][-1]
    gsub('^"|"$', "", v)
  })
  names(vals) <- make.unique(keys)
  n <- max(vapply(vals, length, 1L))
  vals <- lapply(vals, function(v) { length(v) <- n; v })
  as.data.frame(vals, stringsAsFactors = FALSE, check.names = FALSE)
}

pick <- function(meta, pattern) {
  k <- grep(pattern, colnames(meta), value = TRUE)
  if (length(k) == 0) return(rep(NA_character_, nrow(meta)))
  meta[[k[1]]]
}

pickv <- function(meta, prefix) {
  ch <- grep("characteristics_ch1", colnames(meta), value = TRUE)
  for (k in ch) {
    v <- meta[[k]]
    vv <- v[!is.na(v) & v != ""]
    if (length(vv) && mean(startsWith(vv, prefix)) > 0.8) {
      return(sub(paste0("^", prefix, "\\s*"), "", ifelse(is.na(v), NA_character_, v)))
    }
  }
  rep(NA_character_, nrow(meta))
}

qc <- list()

## ---- GSE190609 (RNA-seq, Ensembl, 主队列) ----
m1 <- parse_meta(file.path(raw, "GSE190609_series_matrix.txt.gz"))
meta1 <- data.frame(
  gsm     = pick(m1, "geo_accession"),
  title   = pick(m1, "title"),
  organ   = pick(m1, "source_name_ch1"),
  patient = pickv(m1, "patient:"),
  tnm     = pickv(m1, "tnmstage:"),
  driver  = pickv(m1, "driver_mutation:"),
  stringsAsFactors = FALSE
)
e1 <- fread(cmd = sprintf("gunzip -c '%s'", file.path(raw, "GSE190609_Laoukili_readCounts_rpkm_rv.txt.gz")))
e1_ids <- colnames(e1)[-1]
meta1 <- meta1[match(e1_ids, meta1$title), ]
rownames(meta1) <- NULL
qc$GSE190609 <- list(genes = nrow(e1), samples_expr = length(e1_ids),
                     aligned = sum(!is.na(meta1$gsm)),
                     organ_counts = table(meta1$organ))
fwrite(meta1, file.path(out, "GSE190609_meta.tsv"), sep = "\t")
fwrite(e1, file.path(out, "GSE190609_expr_RPKM.tsv.gz"), sep = "\t")

## ---- GSE50760 (RNA-seq, 基因符号 FPKM, 每样本单文件) ----
m2 <- parse_meta(file.path(raw, "GSE50760_series_matrix.txt.gz"))
meta2 <- data.frame(
  gsm     = pick(m2, "geo_accession"),
  title   = pick(m2, "title"),
  tissue  = pickv(m2, "tissue:"),
  stringsAsFactors = FALSE
)
d2 <- file.path(raw, "GSE50760_RAW")
fs <- list.files(d2, pattern = "_FPKM.txt.gz$", full.names = TRUE)
mats <- lapply(fs, function(f) {
  x <- fread(cmd = sprintf("gunzip -c '%s'", f))
  setnames(x, c("gene", "val"))
  x <- x[, .(val = mean(val, na.rm = TRUE)), by = gene]
  setnames(x, "val", sub("_FPKM.txt.gz$", "", basename(f)))
  x
})
expr2 <- Reduce(function(a, b) merge(a, b, by = "gene", all = FALSE), mats)
colnames(expr2)[-1] <- sub("^([A-Z0-9]+)_.*$", "\\1", colnames(expr2)[-1])
meta2 <- meta2[match(colnames(expr2)[-1], meta2$gsm), ]
qc$GSE50760 <- list(genes = nrow(expr2), samples_expr = ncol(expr2) - 1L,
                    aligned = sum(!is.na(meta2$title)),
                    tissue_counts = table(meta2$tissue))
fwrite(meta2, file.path(out, "GSE50760_meta.tsv"), sep = "\t")
fwrite(expr2, file.path(out, "GSE50760_expr_FPKM.tsv.gz"), sep = "\t")

## ---- GSE225182 (RNA-seq, log2CPM xlsx) ----
m3 <- parse_meta(file.path(raw, "GSE225182_series_matrix.txt.gz"))
meta3 <- data.frame(
  gsm     = pick(m3, "geo_accession"),
  title   = pick(m3, "title"),
  tissue  = pickv(m3, "tissue:"),
  subject = pickv(m3, "subject:"),
  stringsAsFactors = FALSE
)
e3 <- as.data.table(read_excel(file.path(raw, "GSE225182_PMCRClog2CPM.xlsx")))
suffix <- c("metastasized tumor into peritoneal" = "MT",
            "Normal colorectal mucosae" = "N",
            "Primary colorectal tumor" = "PT")
meta3$xlsx_col <- paste0(meta3$subject, "_", unname(suffix[meta3$tissue]))
xlsx_cols <- colnames(e3)[-1]
qc$GSE225182 <- list(genes = nrow(e3), samples_expr = ncol(e3) - 1L,
                     meta_n = nrow(meta3),
                     tissue_counts = table(meta3$tissue),
                     xlsx_aligned = sum(meta3$xlsx_col %in% xlsx_cols),
                     unmatched_meta = meta3$xlsx_col[!meta3$xlsx_col %in% xlsx_cols],
                     unmatched_xlsx = xlsx_cols[!xlsx_cols %in% meta3$xlsx_col])
fwrite(meta3, file.path(out, "GSE225182_meta.tsv"), sep = "\t")
fwrite(e3, file.path(out, "GSE225182_expr_log2CPM.tsv.gz"), sep = "\t")

## ---- QC 汇总 ----
sink(file.path(out, "M1_QC_report.txt"))
for (nm in names(qc)) { cat("=====", nm, "=====\n"); print(qc[[nm]]); cat("\n") }
sink()
print(qc)
cat("\nsessionInfo:\n"); print(sessionInfo())
