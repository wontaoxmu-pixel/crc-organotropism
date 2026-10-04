# A8 修订分析：已发表签名与 PM_core / LM_core 的超几何重叠检验
# 审稿意见 RC-M1/S2：本研究签名是否与已发表签名重叠，增量何在
#
# 输入：
#   data/processed/crc_organotropism/meta/peritoneal_tropism_genes.tsv  (PM_core, 49)
#   data/processed/crc_organotropism/meta/liver_tropism_genes.tsv       (LM_core, 57)
#   data/processed/crc_organotropism/meta/meta_PM_vs_PT_RNAseq.tsv.gz   (PM meta 宇宙)
#   data/processed/crc_organotropism/meta/meta_LM_vs_PT_mixed.tsv.gz    (LM meta 宇宙)
#   results/revision/external_signatures.tsv (已发表签名基因清单，来源见下)
#
# 签名来源（2026-10-04 实际下载并解析）：
#   Gelli2023_PM / Gelli2023_LM : Gelli et al., Cancers 2023;15:4418, Table S1
#       (61 DEG, 原发灶 CRC-Peritoneum vs CRC-Liver; logFC>0 = 腹膜侧 40 个, <0 = 肝侧 21 个)
#       PMC10648258 补充 PDF cancers-2534275-supplementary.pdf (Europe PMC supplementaryFiles)
#   Laoukili2022_PM366 / PT138 : Laoukili et al., Br J Cancer 2022;126:1824-33, Table S2
#       (GSE190609 源文献; PM 高于配对 PT 366 基因 / PT 高于 PM 138 基因; PMC9174226 MOESM2.docx)
#   Lenos2022_PMhigh_up : Lenos et al., Nat Commun 2022;13:4443, Supplementary Data 5
#       (in vivo PM (CCLE) PM_high vs PM_low 上调 855 基因; PMC9352687 MOESM7.xlsx)
#   Yu2025_PM_malignant / Yu2025_LM_malignant : Yu et al., Cancer Lett 2025;629:217880
#       (scRNA 恶性细胞亚群 marker 并集: PM 6 群 805 基因 / LM 6 群 1060 基因;
#        els-cdn 1-s2.0-S0304383525004483-mmc5.xlsx; FindAllMarkers 输出)
#
# 背景宇宙（声明）：主分析 N = 本 meta 实际可检出基因
#   （PM 检验用 PM meta 宇宙 = GSE190609 ∩ GSE225182 共同 symbol；LM 检验用 LM meta 宇宙
#    = GSE190609 ∩ GSE50760 ∩ GSE41258 ∩ GSE41568 共同 symbol）；
#   敏感性分析 N = 20,000 蛋白编码基因近似值。
# 符号匹配：大小写不敏感（toupper）；Laoukili 旧式克隆号逗号小数已在提取时转为点号。
set.seed(123)
suppressPackageStartupMessages(library(data.table))

proc <- "data/processed/crc_organotropism"
outdir <- "results/revision"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

pm_core <- fread(file.path(proc, "meta", "peritoneal_tropism_genes.tsv"))$gene
lm_core <- fread(file.path(proc, "meta", "liver_tropism_genes.tsv"))$gene
uni_pm  <- fread(file.path(proc, "meta", "meta_PM_vs_PT_RNAseq.tsv.gz"))$gene
uni_lm  <- fread(file.path(proc, "meta", "meta_LM_vs_PT_mixed.tsv.gz"))$gene
sig_dt  <- fread(file.path(outdir, "external_signatures.tsv"))

cat("PM_core:", length(pm_core), " LM_core:", length(lm_core), "\n")
cat("Universe PM:", length(uni_pm), " LM:", length(uni_lm), "\n")

up <- toupper
uni_pm_u <- unique(up(uni_pm)); uni_lm_u <- unique(up(uni_lm))
pm_core_u <- unique(up(pm_core)); lm_core_u <- unique(up(lm_core))
stopifnot(all(pm_core_u %in% uni_pm_u), all(lm_core_u %in% uni_lm_u))

sig_list <- split(sig_dt$gene, sig_dt$signature)
sig_list <- lapply(sig_list, function(g) unique(up(g)))

one_test <- function(sig_genes, core_u, uni_u, N_label, N_alt = 20000) {
  K <- sum(sig_genes %in% uni_u)          # 签名中落在宇宙内的基因数
  x <- sum(sig_genes %in% core_u)         # 重叠数
  n <- length(core_u); N <- length(uni_u)
  expected <- n * K / N
  p <- if (x > 0) phyper(x - 1, K, N - K, n, lower.tail = FALSE) else 1
  # 敏感性：N=20000（K、x 不变，近似全蛋白编码背景）
  expected20 <- n * K / N_alt
  p20 <- if (x > 0) phyper(x - 1, K, N_alt - K, n, lower.tail = FALSE) else 1
  data.table(universe = N_label, N = N, sig_total = length(sig_genes),
             sig_in_universe = K, core_size = n, overlap = x,
             expected = round(expected, 3),
             fold_enrichment = ifelse(expected > 0, round(x / expected, 2), NA),
             hypergeom_p = signif(p, 4),
             expected_N20000 = round(expected20, 3),
             hypergeom_p_N20000 = signif(p20, 4))
}

res <- rbindlist(lapply(names(sig_list), function(s) {
  g <- sig_list[[s]]
  r1 <- one_test(g, pm_core_u, uni_pm_u, "meta_PM_detectable")
  r1[, `:=`(signature = s, core = "PM_core",
            overlap_genes = paste(pm_core[toupper(pm_core) %in% g], collapse = ","))]
  r2 <- one_test(g, lm_core_u, uni_lm_u, "meta_LM_detectable")
  r2[, `:=`(signature = s, core = "LM_core",
            overlap_genes = paste(lm_core[toupper(lm_core) %in% g], collapse = ","))]
  rbind(r1, r2)
}))
res[, fdr_BH := signif(p.adjust(hypergeom_p, method = "BH"), 4)]
setcolorder(res, c("signature", "core"))
fwrite(res, file.path(outdir, "signature_overlap.tsv"), sep = "\t")
print(res[, .(signature, core, sig_total, sig_in_universe, core_size, overlap,
              expected, fold_enrichment, hypergeom_p, fdr_BH)],
      right = FALSE)

cat("\nsessionInfo:\n"); print(sessionInfo())
