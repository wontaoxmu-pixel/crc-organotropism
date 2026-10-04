# M2: CMS 分型 + CMS4 × 腹膜富集阳性对照（GSE190609, n=113）
# 预期（Saris 2025 / Guinney 2015）：腹膜转移显著富集 CMS4
set.seed(123)
suppressMessages({
  library(data.table); library(CMSclassifier)
  library(org.Hs.eg.db); library(AnnotationDbi)
})
proc <- "data/processed/crc_organotropism"

e <- fread(file.path(proc, "GSE190609_expr_RPKM.tsv.gz"))
ens <- sub("\\..*$", "", e[[1]])
emat <- log2(as.matrix(e[, -1], rownames = ens) + 0.5)

# Ensembl -> Entrez（CMSclassifier RF 模型要求 Entrez 行名）
ent <- mapIds(org.Hs.eg.db, keys = rownames(emat), column = "ENTREZID",
              keytype = "ENSEMBL", multiVals = "first")
keep <- !is.na(ent)
emat <- emat[keep, ]; ent <- ent[keep]
emat <- emat[!duplicated(ent), ]
rownames(emat) <- ent[!duplicated(ent)]
cat("CMS 输入矩阵:", nrow(emat), "Entrez 基因 x", ncol(emat), "样本\n")

# 手工复现 classifyCMS.RF（包内 naImpute 对 matrix 输入有 names bug）
g <- listModelGenes("RF")
fm <- get("finalModel", envir = as.environment("package:CMSclassifier"))
mm <- mean(rowMeans(emat, na.rm = TRUE))
missing <- setdiff(g, rownames(emat))
if (length(missing) > 0) {
  add <- matrix(mm, nrow = length(missing), ncol = ncol(emat),
                dimnames = list(missing, colnames(emat)))
  emat <- rbind(emat, add)
}
cat("模型基因覆盖:", length(g) - length(missing), "/", length(g), "\n")
Ec <- emat[g, ]
Ec <- Ec - rowMeans(Ec)
prob <- predict(fm, t(Ec), type = "prob")
cms <- apply(prob, 1, function(z) paste(c("CMS1","CMS2","CMS3","CMS4")[z == max(z)], collapse = ","))
cms[apply(prob, 1, max) < 0.5] <- NA
cms_nearest <- apply(prob, 1, function(z) paste(c("CMS1","CMS2","CMS3","CMS4")[z == max(z)], collapse = ","))

# 逐样本后验概率矩阵落盘（补充图 S1 复现 minPosterior >= 0.5 筛选用；纯新增输出）
post <- as.data.table(prob, keep.rownames = "title")
setnames(post, c("title", "CMS1", "CMS2", "CMS3", "CMS4"))
post[, `:=`(nearestCMS = cms_nearest[title],
            max_posterior = pmax(CMS1, CMS2, CMS3, CMS4))]
fwrite(post, file.path(proc, "GSE190609_CMS_posteriors.tsv"), sep = "\t")
cat("后验概率矩阵已落盘: GSE190609_CMS_posteriors.tsv (", nrow(post), "样本 )\n", sep = "")

m <- fread(file.path(proc, "GSE190609_meta.tsv"))
m <- m[match(names(cms), m$title), ]
if (any(is.na(m$organ))) {  # 尝试 gsm 匹配兜底
  m <- fread(file.path(proc, "GSE190609_meta.tsv"))[match(names(cms), gsm), ]
}
stopifnot(all(!is.na(m$organ)))
m$cms <- cms

tab <- table(m$organ, m$cms)
print(tab)
fwrite(m[, .(gsm, title, organ, patient, cms)],
       file.path(proc, "GSE190609_CMS.tsv"), sep = "\t")

# 全样本 nearestCMS 版本（无后验阈值）
m2 <- m; m2$cms <- cms_nearest
cat("== nearestCMS（全 113 样本）==\n")
print(table(m2$organ, m2$cms))
fwrite(m2[, .(gsm, title, organ, patient, cms)],
       file.path(proc, "GSE190609_CMS_nearest.tsv"), sep = "\t")
run_fisher <- function(mm, tag) {
  isPM <- mm$organ == "peritoneal metastasis"
  is4  <- mm$cms == "CMS4"
  isLM <- mm$organ == "liver metastasis"
  t1 <- matrix(c(sum(isPM & is4), sum(isPM & !is4), sum(!isPM & is4), sum(!isPM & !is4)),
               nrow = 2, byrow = TRUE, dimnames = list(c("PM","non-PM"), c("CMS4","non-CMS4")))
  f1 <- fisher.test(t1)
  t2 <- matrix(c(sum(isPM & is4), sum(isPM & !is4), sum(isLM & is4), sum(isLM & !is4)),
               nrow = 2, byrow = TRUE, dimnames = list(c("PM","LM"), c("CMS4","non-CMS4")))
  f2 <- fisher.test(t2)
  cat(sprintf("%s | CMS4×PM(vs其余): OR=%.2f P=%.3g | CMS4×PM(vs肝): OR=%.2f P=%.3g | 腹膜CMS4=%.0f%% 肝CMS4=%.0f%%\n",
              tag, f1$estimate, f1$p.value, f2$estimate, f2$p.value,
              100*mean(is4[isPM]), 100*mean(is4[isLM])))
  print(t1)
}
run_fisher(m2, "nearestCMS")

# CMS4 × 腹膜 vs 其他器官 Fisher（仅用可分型样本）
ok <- !is.na(m$cms)
cat("可分型样本:", sum(ok), "/", nrow(m), "\n")
isPM <- m$organ[ok] == "peritoneal metastasis"
is4  <- m$cms[ok] == "CMS4"
pm4 <- matrix(c(sum(isPM & is4), sum(isPM & !is4),
                sum(!isPM & is4), sum(!isPM & !is4)),
              nrow = 2, byrow = TRUE,
              dimnames = list(c("PM", "non-PM"), c("CMS4", "non-CMS4")))
print(pm4)
ft <- fisher.test(pm4)
cat(sprintf("CMS4 × 腹膜富集: OR=%.2f, P=%.2e\n", ft$estimate, ft$p.value))
# 腹膜内 CMS4 比例 vs 肝内 CMS4 比例
cat("腹膜 CMS4 占比:", round(mean(m$cms[ok][isPM] == "CMS4"), 3),
    " 肝 CMS4 占比:", round(mean(m$cms[ok][m$organ[ok] == "liver metastasis"] == "CMS4"), 3), "\n")
