# M2 层次1：GSE190609 同患者多器官配对差异分析（limma-trend）
# 主对比: peritoneal vs liver / peritoneal vs primary / liver vs primary
set.seed(123)
suppressPackageStartupMessages({
  library(data.table); library(limma); library(org.Hs.eg.db); library(AnnotationDbi)
})

out <- "data/processed/crc_organotropism"
dir.create(file.path(out, "de"), showWarnings = FALSE, recursive = TRUE)

meta <- fread(file.path(out, "GSE190609_meta.tsv"))
expr <- fread(file.path(out, "GSE190609_expr_RPKM.tsv.gz"))
genes <- expr[[1]]
emat <- log2(as.matrix(expr[, -1], rownames = genes) + 0.5)

stopifnot(all(colnames(emat) == meta$title))
meta$organ <- factor(meta$organ,
  levels = c("primary tumor", "liver metastasis", "lymph node metastasis",
             "peritoneal metastasis", "ovarium metastasis"))

# 低表达过滤: RPKM >= 1 (log2(1.5)=0.585) 在 >= 10 个样本
keep <- rowSums(emat >= log2(1.5)) >= 10
emat_f <- emat[keep, ]
cat("genes after filter:", nrow(emat_f), "/", nrow(emat), "\n")

# Ensembl -> symbol
ens <- rownames(emat_f)
sym <- mapIds(org.Hs.eg.db, keys = ens, column = "SYMBOL",
              keytype = "ENSEMBL", multiVals = "first")
cat("mapped symbols:", sum(!is.na(sym)), "/", length(sym), "\n")

design <- model.matrix(~ 0 + organ + patient, data = meta)
colnames(design) <- make.names(colnames(design))
colnames(design) <- sub("^organ", "", colnames(design))

fit <- lmFit(emat_f, design)
cm <- makeContrasts(
  PM_vs_PT  = peritoneal.metastasis - primary.tumor,
  LM_vs_PT  = liver.metastasis - primary.tumor,
  PM_vs_LM  = peritoneal.metastasis - liver.metastasis,
  levels = design)
fit2 <- eBayes(contrasts.fit(fit, cm), trend = TRUE)

res <- list()
for (cn in colnames(cm)) {
  tt <- topTable(fit2, coef = cn, number = Inf, sort.by = "P")
  tt$ensembl <- rownames(tt)
  tt$symbol <- sym[rownames(tt)]
  res[[cn]] <- tt
  sig <- tt$adj.P.Val < 0.05 & abs(tt$logFC) >= 1
  cat(sprintf("%s: 总基因 %d, FDR<0.05 且 |logFC|>=1: %d (上调 %d / 下调 %d)\n",
              cn, nrow(tt), sum(sig), sum(sig & tt$logFC > 0), sum(sig & tt$logFC < 0)))
  fwrite(tt, file.path(out, "de", paste0("GSE190609_", cn, ".tsv.gz")), sep = "\t")
}

# sanity: 已知腹膜/间质-EMT 标记在 PM_vs_PT 的方向
markers <- c("VIM","ZEB1","ZEB2","SNAI1","SNAI2","TWIST1","COL1A1","COL3A1",
             "FAP","ACTA2","SPARC","TGFBI","TAGLN","FN1","CDH1","CDH2",
             "EPCAM","MUC1","CLDN3","CLDN4","CEACAM5","KRT8","KRT18")
tt <- res$PM_vs_PT
ttm <- tt[!is.na(tt$symbol) & tt$symbol %in% markers, ]
ttm <- ttm[match(markers[markers %in% ttm$symbol], ttm$symbol), ]
cat("\n== sanity: PM_vs_PT 间质/EMT/上皮标记 ==\n")
print(ttm[, c("symbol","logFC","adj.P.Val")])
fwrite(ttm, file.path(out, "de", "GSE190609_PM_vs_PT_markers.tsv"), sep = "\t")

cat("\nsessionInfo:\n"); print(sessionInfo())
