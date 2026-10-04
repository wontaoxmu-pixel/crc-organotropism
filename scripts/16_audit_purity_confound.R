# 核对用补充分析：趋向基因集的肿瘤纯度/TME 混淆量化
# 1) GSE190609 转移灶中 PM_core/LM_core 得分与上皮含量代理（EPCAM+KRT8/18 均值）的相关
# 2) 基因集细胞来源标注：间质/脂肪/肝细胞 vs 肿瘤上皮内在
set.seed(123)
suppressMessages({library(data.table); library(org.Hs.eg.db); library(AnnotationDbi)})
proc <- "data/processed/crc_organotropism"

pm_core <- fread(file.path(proc, "meta", "peritoneal_tropism_genes.tsv"))$gene
lm_core <- fread(file.path(proc, "meta", "liver_tropism_genes.tsv"))$gene

e <- fread(file.path(proc, "GSE190609_expr_RPKM.tsv.gz"))
ens <- sub("\\..*$", "", e[[1]])
emat <- log2(as.matrix(e[, -1], rownames = ens) + 0.5)
sym <- mapIds(org.Hs.eg.db, keys = rownames(emat), column = "SYMBOL",
              keytype = "ENSEMBL", multiVals = "first")
keep <- !is.na(sym); emat <- emat[keep, ]; sym <- sym[keep]
emat <- emat[!duplicated(sym), ]; rownames(emat) <- sym[!duplicated(sym)]
m <- fread(file.path(proc, "GSE190609_meta.tsv"))

ssgsea_manual <- function(mat, gsets) {
  N <- nrow(mat)
  es <- sapply(gsets, function(gs) {
    inset <- rownames(mat) %in% gs
    Ns <- sum(inset)
    apply(mat, 2, function(x) {
      o <- order(x, decreasing = TRUE)
      sum(cumsum(ifelse(inset[o], 1 / Ns, -1 / (N - Ns))))
    })
  })
  apply(es, 2, function(v) (v - min(v)) / (max(v) - min(v)))
}

met_titles <- m[organ != "primary tumor"]$title
mm <- emat[, met_titles]
sc <- ssgsea_manual(mm, list(PM_core = pm_core, LM_core = lm_core))
# 上皮含量代理：EPCAM/KRT8/KRT19 均值；间质代理：DCN/COL1A1/FAP 均值；肝细胞代理：ALB/APOB
epi <- colMeans(mm[c("EPCAM", "KRT8", "KRT19"), ])
stro <- colMeans(mm[c("DCN", "COL1A1", "FAP"), ])
hepa <- colMeans(mm[c("ALB", "APOB"), ])
org <- m$organ[match(met_titles, m$title)]
df <- data.table(organ = org, PM = sc[, "PM_core"], LM = sc[, "LM_core"],
                 epi = epi, stromal = stro, hepatocyte = hepa)
cat("== GSE190609 转移灶（n=", ncol(mm), "）相关（Spearman）==\n")
cat("PM_core vs 上皮代理:", round(cor(df$PM, df$epi, method = "spearman"), 3),
    " | vs 间质代理:", round(cor(df$PM, df$stromal, method = "spearman"), 3), "\n")
cat("LM_core vs 上皮代理:", round(cor(df$LM, df$epi, method = "spearman"), 3),
    " | vs 肝细胞代理:", round(cor(df$LM, df$hepatocyte, method = "spearman"), 3), "\n")
cat("\n按器官均值:\n")
print(df[, .(PM = round(mean(PM), 2), LM = round(mean(LM), 2), epi = round(mean(epi), 2),
             stromal = round(mean(stromal), 2), hepatocyte = round(mean(hepatocyte), 2)), by = organ])

# 肝细胞代理在肝转移 vs 其他转移中的差异（污染证据）
cat("\nALB/APOB 均值：肝转移", round(mean(df$hepatocyte[df$organ == "liver metastasis"]), 2),
    " vs 腹膜", round(mean(df$hepatocyte[df$organ == "peritoneal metastasis"]), 2),
    " vs 淋巴结", round(mean(df$hepatocyte[df$organ == "lymph node metastasis"]), 2), "\n")
w <- wilcox.test(df$hepatocyte[df$organ == "liver metastasis"],
                 df$hepatocyte[df$organ == "peritoneal metastasis"])
cat("肝细胞代理 肝 vs 腹膜 Wilcoxon P:", format(w$p.value, digits = 2), "\n")
