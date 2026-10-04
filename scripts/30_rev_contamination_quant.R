# 任务 A7（RA-M3）：器官驻留细胞定量——肝/脂肪/间皮标志 score 验证 + LM_core/PM_core 成分分析
# 不依赖新包：z-score 标志评分 + 线性归因模型；xCell 未安装（跳过，报告中注明）
set.seed(123)
suppressMessages({
  library(data.table); library(org.Hs.eg.db); library(AnnotationDbi)
})
proc <- "data/processed/crc_organotropism"
dir.create("results/revision", showWarnings = FALSE, recursive = TRUE)

## ---- 1. 数据加载（与 scripts/18_fix_A_purity.R 一致） ----
e <- fread(file.path(proc, "GSE190609_expr_RPKM.tsv.gz"))
ens <- sub("\\..*$", "", e[[1]])
emat <- log2(as.matrix(e[, -1], rownames = ens) + 0.5)
m <- fread(file.path(proc, "GSE190609_meta.tsv"))
m <- m[match(colnames(emat), m$title)]
sym <- mapIds(org.Hs.eg.db, keys = rownames(emat), column = "SYMBOL",
              keytype = "ENSEMBL", multiVals = "first")
keep <- rowSums(emat >= log2(1.5)) >= 10
emat <- emat[keep, ]; sym <- sym[keep]
emat <- emat[!is.na(sym), ]; sym <- sym[!is.na(sym)]
emat <- emat[!duplicated(sym), ]; rownames(emat) <- sym[!duplicated(sym)]
cat(sprintf("表达矩阵: %d 基因 x %d 样本\n", nrow(emat), ncol(emat)))

## ---- 2. 器官驻留细胞标志 score（基因跨样本 z-score 后按样本取均值） ----
mk <- list(
  hepatocyte = c("ALB", "APOB", "APOA1", "APOA2", "SERPINA1", "HP", "FGG",
                 "TTR", "GC", "ORM1", "CYP2E1", "CYP3A4"),
  adipocyte  = c("ADIPOQ", "PLIN1", "FABP4", "CD36", "LPL", "LIPE",
                 "CIDEC", "ADIRF"),
  mesothelial = c("MSLN", "WT1", "UPK3B", "CALB2", "LRRN4"),
  epithelial = c("EPCAM", "KRT8", "KRT19")
)
z <- t(scale(t(emat)))  # 每基因跨样本 z-score
scores <- sapply(mk, function(gs) {
  g <- intersect(gs, rownames(z))
  cat(sprintf("marker 可用: %s\n", paste(names(mk)[sapply(mk, identical, gs)], paste(g, collapse = ","))))
  colMeans(z[g, , drop = FALSE])
})
scores <- as.data.frame(scores)
m <- cbind(m, scores[colnames(emat), , drop = FALSE])

out <- list()
rec <- function(section, item, value, detail = "")
  out[[length(out) + 1]] <<- data.table(section = section, item = item,
                                        value = value, detail = detail)

for (s in names(mk)) {
  present <- intersect(mk[[s]], rownames(z))
  rec("marker_coverage", s, paste0(length(present), "/", length(mk[[s]])),
      paste(present, collapse = ","))
}

## ---- 3. score 器官特异性验证（GSE190609） ----
org_lv <- c("primary tumor", "peritoneal metastasis", "liver metastasis",
            "lymph node metastasis", "ovarium metastasis")
for (s in c("hepatocyte", "adipocyte", "mesothelial", "epithelial")) {
  v <- m[[s]]; g <- factor(m$organ, levels = org_lv)
  kw <- kruskal.test(v, g)
  med <- tapply(v, g, median)
  rec("score_specificity", paste0(s, "_KW_p"), signif(kw$p.value, 4),
      paste(paste(names(med), round(med, 2), sep = "="), collapse = "; "))
  # 关键对比
  cmp <- if (s == "hepatocyte") c("liver metastasis", "primary tumor") else c("peritoneal metastasis", "primary tumor")
  w <- wilcox.test(v[m$organ == cmp[1]], v[m$organ == cmp[2]])
  d_med <- median(v[m$organ == cmp[1]]) - median(v[m$organ == cmp[2]])
  rec("score_specificity", paste0(s, "_", sub(" .*", "", cmp[1]), "_vs_PT"),
      sprintf("dmedian=%.2f, wilcox_p=%.3g", d_med, w$p.value), "")
}

## ---- 4. LM_core / PM_core 成分分类（知识注释 + 数据驱动肝/脂肪特异性） ----
lm_core <- fread(file.path(proc, "meta", "liver_tropism_genes.tsv"))
pm_core <- fread(file.path(proc, "meta", "peritoneal_tropism_genes.tsv"))

# 知识分类：肝分泌蛋白/急性期、肝代谢/异生物质代谢、其他
hep_secreted <- c("GC", "FGG", "FGA", "FGB", "F9", "F10", "CPB2", "APOA1", "APOA2",
                  "APOC3", "APOB", "APOM", "APOH", "HRG", "MBL2", "ORM1", "ORM2",
                  "C9", "C5", "C6", "C4BPA", "CP", "HP", "HPR", "ITIH3", "FGL1",
                  "IGFBP1", "COLEC11", "SERPINA1", "SERPINA3", "SERPINA5", "MST1",
                  "TTR", "PON3", "ASGR2")
hep_metab <- c("CYP2E1", "CYP2C8", "CYP4F3", "AKR1D1", "ARG1", "ALDH8A1", "GYS2",
               "FMO3", "PAH", "RGN", "CBS", "GATM", "CPS1", "AOX1", "ABCB4",
               "PIPOX", "ABAT", "MTTP", "HSD11B1", "CES1", "ABCG5", "ALDH8A1")
lm_core[, class := fifelse(gene %in% hep_secreted, "hepatic_secreted_acute_phase",
                    fifelse(gene %in% hep_metab, "hepatic_metabolic_xenobiotic", "other"))]
tab <- lm_core[, .N, by = class]
for (i in seq_len(nrow(tab)))
  rec("LM_core_composition", tab$class[i], sprintf("%d/57 (%.0f%%)", tab$N[i], 100 * tab$N[i] / 57),
      paste(lm_core$gene[lm_core$class == tab$class[i]], collapse = ","))

# 数据驱动特异性：基因表达与各驻留 score 的跨样本 Spearman 相关
gcor <- function(genes, score) sapply(intersect(genes, rownames(emat)),
                                      function(g) cor(emat[g, ], m[[score]], method = "spearman"))
lm_core[, hep_cor := gcor(gene, "hepatocyte")[gene]]
pm_core[, adipo_cor := gcor(gene, "adipocyte")[gene]]
pm_core[, meso_cor := gcor(gene, "mesothelial")[gene]]
pm_core[, hep_cor := gcor(gene, "hepatocyte")[gene]]
lm_core[, adipo_cor := gcor(gene, "adipocyte")[gene]]

# ③ 效应量与肝特异性相关
sp <- cor.test(lm_core$beta, lm_core$hep_cor, method = "spearman", exact = FALSE)
rec("LM_core_effect_vs_liverSpec", "spearman(beta, hep_cor)",
    sprintf("rho=%.3f, p=%.3g, n=%d", sp$estimate, sp$p.value, sum(!is.na(lm_core$hep_cor))), "")
rec("LM_core_hepCor_summary", "median (IQR)",
    sprintf("%.2f (%.2f-%.2f)", median(lm_core$hep_cor, na.rm = TRUE),
            quantile(lm_core$hep_cor, .25, na.rm = TRUE), quantile(lm_core$hep_cor, .75, na.rm = TRUE)),
    sprintf("%d/57 genes rho>0.5", sum(lm_core$hep_cor > 0.5, na.rm = TRUE)))

# ④ PM_core 同框架：脂肪/间皮来源
adipo_genes <- c("FABP4", "CD36", "LPL", "LIPE", "RBP7", "NNAT", "EBF2", "MEDAG", "PLAC9", "VIT")
fibro_ecm  <- c("PI16", "OGN", "MFAP5", "SFRP2", "OLFML2B", "SRPX", "ADAMTS15", "ADAMTS16",
                "COLEC12", "CPXM1", "CHRDL1", "GALNT13", "FLRT2", "SLIT3", "PDGFRL", "C16orf89")
meso_canonical <- c("MSLN", "WT1", "UPK3B", "CALB2")
pm_core[, class := fifelse(gene %in% adipo_genes, "adipocyte_lipid",
                    fifelse(gene %in% meso_canonical, "mesothelial_canonical",
                     fifelse(gene %in% fibro_ecm, "fibroblast_ECM_stromal", "other_neural_misc")))]
tab2 <- pm_core[, .(N = .N, med_beta = round(median(beta), 2)), by = class]
for (i in seq_len(nrow(tab2)))
  rec("PM_core_composition", tab2$class[i],
      sprintf("%d/49 (%.0f%%), median_beta=%.2f", tab2$N[i], 100 * tab2$N[i] / 49, tab2$med_beta[i]),
      paste(pm_core$gene[pm_core$class == tab2$class[i]], collapse = ","))
sp2 <- cor.test(pm_core$beta, pm_core$adipo_cor, method = "spearman", exact = FALSE)
rec("PM_core_effect_vs_adipoSpec", "spearman(beta, adipo_cor)",
    sprintf("rho=%.3f, p=%.3g, n=%d", sp2$estimate, sp2$p.value, sum(!is.na(pm_core$adipo_cor))), "")
rec("PM_core_adipoCor_summary", "median (IQR)",
    sprintf("%.2f (%.2f-%.2f)", median(pm_core$adipo_cor, na.rm = TRUE),
            quantile(pm_core$adipo_cor, .25, na.rm = TRUE), quantile(pm_core$adipo_cor, .75, na.rm = TRUE)),
    sprintf("%d/49 genes rho>0.5", sum(pm_core$adipo_cor > 0.5, na.rm = TRUE)))
rec("PM_core_mesoCor_summary", "median (IQR)",
    sprintf("%.2f (%.2f-%.2f)", median(pm_core$meso_cor, na.rm = TRUE),
            quantile(pm_core$meso_cor, .25, na.rm = TRUE), quantile(pm_core$meso_cor, .75, na.rm = TRUE)),
    sprintf("%d/49 genes rho>0.5", sum(pm_core$meso_cor > 0.5, na.rm = TRUE)))

# 与纯度校正保留清单对照
lm_ret <- fread(file.path(proc, "meta", "LM_core_purity_retained.tsv"))$gene
pm_ret <- fread(file.path(proc, "meta", "PM_core_purity_retained.tsv"))$gene
rec("purity_retention", "LM_core", sprintf("%d/57 (%.0f%%) retained after EPCAM/KRT8/KRT19 adjustment", length(lm_ret), 100 * length(lm_ret) / 57),
    paste(setdiff(lm_core$gene, lm_ret), collapse = ","))
rec("purity_retention", "PM_core", sprintf("%d/49 (%.0f%%) retained after EPCAM/KRT8/KRT19 adjustment", length(pm_ret), 100 * length(pm_ret) / 49),
    paste(setdiff(pm_core$gene, pm_ret), collapse = ","))

## ---- 5. 定量归因：LM_core 信号中可归因于肝细胞混入的比例 ----
# 模型：每基因 expr ~ hepatocyte score（全 113 样本）；predicted ΔLM-PT = slope * ΔhepScore；fraction = predicted / observed
is_lm <- m$organ == "liver metastasis"; is_pt <- m$organ == "primary tumor"
d_hep <- median(m$hepatocyte[is_lm]) - median(m$hepatocyte[is_pt])
g57 <- intersect(lm_core$gene, rownames(emat))
attr_tab <- rbindlist(lapply(g57, function(g) {
  y <- emat[g, ]
  sl <- coef(lm(y ~ m$hepatocyte))[2]
  obs <- mean(y[is_lm]) - mean(y[is_pt])
  pred <- sl * d_hep
  data.table(gene = g, slope = sl, observed_d = obs, predicted_d = pred,
             frac = ifelse(obs > 0, pred / obs, NA_real_))
}))
frac_pos <- attr_tab[observed_d > 0]
rec("LM_attribution", "per_gene_fraction_median",
    sprintf("%.2f (IQR %.2f-%.2f), n=%d", median(frac_pos$frac),
            quantile(frac_pos$frac, .25), quantile(frac_pos$frac, .75), nrow(frac_pos)),
    sprintf("%d/%d genes frac>=1; %d/%d genes frac>=0.5",
            sum(frac_pos$frac >= 1), nrow(frac_pos), sum(frac_pos$frac >= 0.5), nrow(frac_pos)))
# score 层面：LM_core score Δobs vs 由肝 score 预测的 Δ
core_sc <- colMeans(z[g57, , drop = FALSE])
sc_obs <- mean(core_sc[is_lm]) - mean(core_sc[is_pt])
sl_sc <- coef(lm(core_sc ~ m$hepatocyte))[2]
sc_pred <- sl_sc * d_hep
rec("LM_attribution", "score_level",
    sprintf("observed_d=%.2f, predicted_from_hepScore=%.2f, frac=%.2f", sc_obs, sc_pred, sc_pred / sc_obs),
    sprintf("hepScore d(LM-PT)=%.2f (z units)", d_hep))
# 对照：PM_core 用脂肪 score 同法归因（PM vs PT）
is_pm <- m$organ == "peritoneal metastasis"
d_adi <- median(m$adipocyte[is_pm]) - median(m$adipocyte[is_pt])
g49 <- intersect(pm_core$gene, rownames(emat))
attr_pm <- rbindlist(lapply(g49, function(g) {
  y <- emat[g, ]
  sl <- coef(lm(y ~ m$adipocyte))[2]
  obs <- mean(y[is_pm]) - mean(y[is_pt])
  pred <- sl * d_adi
  data.table(gene = g, slope = sl, observed_d = obs, predicted_d = pred,
             frac = ifelse(obs > 0, pred / obs, NA_real_))
}))
frac_pm <- attr_pm[observed_d > 0]
rec("PM_attribution_adipocyte", "per_gene_fraction_median",
    sprintf("%.2f (IQR %.2f-%.2f), n=%d", median(frac_pm$frac),
            quantile(frac_pm$frac, .25), quantile(frac_pm$frac, .75), nrow(frac_pm)),
    sprintf("%d/%d genes frac>=1; %d/%d genes frac>=0.5; adipoScore d(PM-PT)=%.2f",
            sum(frac_pm$frac >= 1), nrow(frac_pm), sum(frac_pm$frac >= 0.5), nrow(frac_pm), d_adi))
core_pm <- colMeans(z[g49, , drop = FALSE])
pm_obs <- mean(core_pm[is_pm]) - mean(core_pm[is_pt])
pm_pred <- coef(lm(core_pm ~ m$adipocyte))[2] * d_adi
rec("PM_attribution_adipocyte", "score_level",
    sprintf("observed_d=%.2f, predicted_from_adipoScore=%.2f, frac=%.2f", pm_obs, pm_pred, pm_pred / pm_obs), "")

# 反向对照：LM 的脂肪 score、PM 的肝 score 不应有同等解释力（特异性检验）
sl_cross <- coef(lm(core_sc ~ m$adipocyte))[2]
d_adi_lm <- median(m$adipocyte[is_lm]) - median(m$adipocyte[is_pt])
rec("LM_attribution_cross", "adipoScore_control",
    sprintf("predicted_from_adipoScore=%.2f vs observed=%.2f (frac=%.2f); adipoScore d(LM-PT)=%.2f",
            sl_cross * d_adi_lm, sc_obs, (sl_cross * d_adi_lm) / sc_obs, d_adi_lm), "")

fwrite(rbindlist(out), "results/revision/contamination_quant.tsv", sep = "\t")
fwrite(attr_tab, "results/revision/contamination_quant_LM_per_gene.tsv", sep = "\t")
fwrite(attr_pm, "results/revision/contamination_quant_PM_per_gene.tsv", sep = "\t")
fwrite(lm_core[, .(gene, beta, fdr, class, hep_cor)], "results/revision/LM_core_annotation.tsv", sep = "\t")
fwrite(pm_core[, .(gene, beta, fdr, class, adipo_cor, meso_cor)], "results/revision/PM_core_annotation.tsv", sep = "\t")

cat("\n==== 摘要 ====\n")
print(rbindlist(out))
