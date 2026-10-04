# Revision A1: 患者级 CMS4 敏感性分析（GSE190609）
# 回应审稿意见 RA-M2 / RC-M4 / RA-m7：样本级 59 PM 来自 12 患者，存在伪重复风险
# 内容：
#   (1) 队列结构核对 + 三对比实际进入模型的患者数/样本数 + scripts/11 design 复核
#   (2) 患者级 CMS4：每患者 PM CMS4 比例、>=1 CMS4 PM、LM 来源患者
#   (3) 患者级检验：配对 Wilcoxon、McNemar 精确检验、GLMM (lme4)
#   (4) limma 患者阻断敏感性：duplicateCorrelation vs 固定效应（PM_vs_PT）
# set.seed(123)
set.seed(123)
suppressPackageStartupMessages({
  library(data.table); library(limma); library(lme4)
})

proc <- "data/processed/crc_organotropism"
outdir <- "results/revision"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

meta  <- fread(file.path(proc, "GSE190609_meta.tsv"))
cmsN  <- fread(file.path(proc, "GSE190609_CMS_nearest.tsv"))   # 无阈值，113 全分型（含平局）
cmsT  <- fread(file.path(proc, "GSE190609_CMS.tsv"))           # 后验>=0.5，低置信为空
cmsT[cms == "", cms := NA_character_]

m <- merge(meta, cmsN[, .(title, cmsN = cms)], by = "title")
m <- merge(m, cmsT[, .(title, cmsT = cms)], by = "title")
m$patient <- factor(m$patient)
stopifnot(nrow(m) == 113)

# 平局样本（nearestCMS 仍给出 "CMSx,CMSy"）
ties <- m[grepl(",", cmsN)]
cat("nearestCMS 平局样本数:", nrow(ties), "\n")
print(ties[, .(title, organ, patient, cmsN)])

m[, isPM := organ == "peritoneal metastasis"]
m[, isPT := organ == "primary tumor"]
m[, isLM := organ == "liver metastasis"]
# 与 scripts/13 一致：严格等值判定 CMS4（平局不算 CMS4）
m[, cms4N := cmsN == "CMS4"]
m[, cms4T := !is.na(cmsT) & cmsT == "CMS4"]

results <- list()  # 汇总行
add <- function(section, analysis, term, n_pat = NA, n_samp = NA,
                est = NA, lo = NA, hi = NA, p = NA, detail = "") {
  results[[length(results) + 1L]] <<- data.table(
    section = section, analysis = analysis, term = term,
    n_patients = n_pat, n_samples = n_samp, estimate = est,
    ci_lo = lo, ci_hi = hi, p_value = p, detail = detail)
}

## ---------- (1) 队列结构 + 三对比进入模型的患者/样本数 ----------
org_tab <- m[, .N, by = organ][order(-N)]
print(org_tab)
org_pat <- m[, .(n_samples = .N), by = .(organ, patient)][order(organ, patient)]
cat("\norgan x patient:\n"); print(org_pat)

# scripts/11 design 复核：~0 + organ + patient，患者固定效应阻断（此处仅验秩，与行序无关）
meta11 <- copy(m)
meta11$organ <- factor(meta11$organ,
  levels = c("primary tumor", "liver metastasis", "lymph node metastasis",
             "peritoneal metastasis", "ovarium metastasis"))
design11 <- model.matrix(~ 0 + organ + patient, data = meta11)
rk <- qr(design11)$rank
cat(sprintf("\nscripts/11 design: %d 样本, %d 参数, 秩=%d (%s)\n",
            nrow(design11), ncol(design11), rk,
            ifelse(rk == ncol(design11), "满秩", "秩亏!")))
add("design_check", "scripts11_design", "~0+organ+patient",
    n_pat = nlevels(m$patient), n_samp = nrow(design11),
    est = rk, detail = paste0("ncol=", ncol(design11),
      "; 患者作为加性固定效应阻断；同患者多区域(含35个PT区域)共享患者截距；",
      ifelse(rk == ncol(design11), "满秩可估计", "秩亏")))

# 三对比：直接贡献两器官水平的样本/患者；模型拟合用全部113样本
contrasts <- list(PM_vs_PT = c("peritoneal metastasis", "primary tumor"),
                  LM_vs_PT = c("liver metastasis", "primary tumor"),
                  PM_vs_LM = c("peritoneal metastasis", "liver metastasis"))
for (cn in names(contrasts)) {
  oo <- contrasts[[cn]]
  sub <- m[organ %in% oo]
  n_s <- nrow(sub); n_p <- uniqueN(sub$patient)
  n_s1 <- nrow(m[organ == oo[1]]); n_p1 <- uniqueN(m[organ == oo[1]]$patient)
  n_s2 <- nrow(m[organ == oo[2]]); n_p2 <- uniqueN(m[organ == oo[2]]$patient)
  cat(sprintf("%s: 直接相关 %d 样本/%d 患者 (%s %d样本/%d患者 vs %s %d样本/%d患者); 全模型113样本/12患者\n",
              cn, n_s, n_p, oo[1], n_s1, n_p1, oo[2], n_s2, n_p2))
  add("contrast_composition", cn, "samples/patients contributing to contrast",
      n_pat = n_p, n_samp = n_s,
      detail = sprintf("%s:%d样本/%d患者; %s:%d样本/%d患者; limma模型同时拟合全部113样本/12患者",
                       oo[1], n_s1, n_p1, oo[2], n_s2, n_p2))
}

## ---------- (2) 患者级 CMS4 描述 ----------
per_pat_pm <- m[isPM == TRUE, .(
  n_PM = .N,
  n_PM_CMS4_N = sum(cms4N),
  prop_CMS4_N = mean(cms4N),
  n_PM_typable_T = sum(!is.na(cmsT)),
  n_PM_CMS4_T = sum(cms4T),
  prop_CMS4_T = sum(cms4T) / sum(!is.na(cmsT))
), by = patient][order(patient)]
per_pat_pm[, any_CMS4_N := n_PM_CMS4_N >= 1]
per_pat_pm[, majority_CMS4_N := prop_CMS4_N >= 0.5]
cat("\n== 每患者 PM CMS4（nearestCMS）==\n"); print(per_pat_pm)

per_pat_pt <- m[isPT == TRUE, .(
  n_PT = .N, n_PT_CMS4_N = sum(cms4N), prop_CMS4_PT_N = mean(cms4N)
), by = patient][order(patient)]
per_pat_pt[, majority_CMS4_PT_N := prop_CMS4_PT_N >= 0.5]

lm_tab <- m[isLM == TRUE, .(title, patient, cmsN, cmsT)]
cat("\n== LM 样本（6个）==\n"); print(lm_tab)
n_lm_pat <- uniqueN(lm_tab$patient)

# 写 per-patient 明细
for (i in seq_len(nrow(per_pat_pm))) {
  r <- per_pat_pm[i]
  add("per_patient_PM", paste0("patient_", r$patient), "PM CMS4 (nearest)",
      n_pat = 1, n_samp = r$n_PM, est = round(r$prop_CMS4_N, 4),
      detail = sprintf("n_PM=%d; n_CMS4=%d; any_CMS4=%s; majority_CMS4=%s; threshold版 %d/%d typable CMS4",
                       r$n_PM, r$n_PM_CMS4_N, r$any_CMS4_N, r$majority_CMS4_N,
                       r$n_PM_CMS4_T, r$n_PM_typable_T))
}
add("summary", "PM_CMS4_per_patient", "patients with >=1 CMS4 PM",
    n_pat = nrow(per_pat_pm), est = sum(per_pat_pm$any_CMS4_N),
    detail = sprintf("%d/%d 患者至少1个CMS4 PM; %d/%d 患者PM中CMS4占多数(>=50%%)",
                     sum(per_pat_pm$any_CMS4_N), nrow(per_pat_pm),
                     sum(per_pat_pm$majority_CMS4_N), nrow(per_pat_pm)))
add("summary", "PM_CMS4_proportions", "patient-level mean vs sample-level",
    est = round(mean(per_pat_pm$prop_CMS4_N), 4),
    detail = sprintf("患者级均值=%.3f; 样本级=%d/%d=%.3f; LM样本级=%d/%d=%.3f; LM来自%d个患者",
                     mean(per_pat_pm$prop_CMS4_N),
                     sum(m$isPM & m$cms4N), sum(m$isPM),
                     mean(m$cms4N[m$isPM]),
                     sum(m$isLM & m$cms4N), sum(m$isLM),
                     mean(m$cms4N[m$isLM]), n_lm_pat))

## ---------- (3) 患者级检验 ----------
# 合并每患者 PM/PT 比例（nearestCMS）
pp <- merge(per_pat_pm[, .(patient, prop_PM = prop_CMS4_N, maj_PM = majority_CMS4_N)],
            per_pat_pt[, .(patient, prop_PT = prop_CMS4_PT_N, maj_PT = majority_CMS4_PT_N)],
            by = "patient")

# 3a. PM vs PT 配对 Wilcoxon（12 患者）
wt1 <- suppressWarnings(wilcox.test(pp$prop_PM, pp$prop_PT, paired = TRUE, exact = FALSE))
cat(sprintf("\nPM vs PT 患者级配对Wilcoxon: V=%.1f, P=%.4f (n=%d 患者)\n",
            wt1$statistic, wt1$p.value, nrow(pp)))
add("patient_level_test", "PM_vs_PT", "paired Wilcoxon on per-patient CMS4 proportion",
    n_pat = nrow(pp), est = round(median(pp$prop_PM - pp$prop_PT), 4),
    p = signif(wt1$p.value, 4),
    detail = sprintf("中位差(PM-PT)=%.3f; %d正/%d负/%d零",
                     median(pp$prop_PM - pp$prop_PT),
                     sum(pp$prop_PM > pp$prop_PT), sum(pp$prop_PM < pp$prop_PT),
                     sum(pp$prop_PM == pp$prop_PT)))

# 3b. PM vs PT McNemar 精确（多数派二分）
disc <- pp[maj_PM != maj_PT]
b1 <- sum(pp$maj_PM & !pp$maj_PT); c1 <- sum(!pp$maj_PM & pp$maj_PT)
mc1 <- binom.test(b1, b1 + c1)
cat(sprintf("PM vs PT McNemar精确: b=%d c=%d, P=%.4f\n", b1, c1, mc1$p.value))
add("patient_level_test", "PM_vs_PT", "McNemar exact on majority-CMS4",
    n_pat = nrow(pp), p = signif(mc1$p.value, 4),
    detail = sprintf("PM多数CMS4且PT非=%d; 反向=%d", b1, c1))

# 3c. PM vs LM（6 个有 LM 的患者）配对 Wilcoxon
lm6 <- lm_tab[, .(lm_cms4 = any(cmsN == "CMS4")), by = patient]
pp6 <- merge(per_pat_pm[, .(patient, prop_PM = prop_CMS4_N, maj_PM = majority_CMS4_N)],
             lm6, by = "patient")
wt2 <- suppressWarnings(wilcox.test(pp6$prop_PM, as.numeric(pp6$lm_cms4),
                                    paired = TRUE, exact = FALSE))
cat(sprintf("PM vs LM 患者级配对Wilcoxon: V=%.1f, P=%.4f (n=%d 患者, LM各1样本)\n",
            wt2$statistic, wt2$p.value, nrow(pp6)))
add("patient_level_test", "PM_vs_LM", "paired Wilcoxon (PM prop vs LM 0/1)",
    n_pat = nrow(pp6), est = round(median(pp6$prop_PM - as.numeric(pp6$lm_cms4)), 4),
    p = signif(wt2$p.value, 4),
    detail = sprintf("LM样本=%d 来自%d患者, CMS4 LM=%d; PM比例>LM: %d/%d",
                     nrow(lm_tab), n_lm_pat, sum(lm_tab$cmsN == "CMS4"),
                     sum(pp6$prop_PM > as.numeric(pp6$lm_cms4)), nrow(pp6)))
b2 <- sum(pp6$maj_PM & !pp6$lm_cms4); c2 <- sum(!pp6$maj_PM & pp6$lm_cms4)
if (b2 + c2 > 0) {
  mc2 <- binom.test(b2, b2 + c2)
  add("patient_level_test", "PM_vs_LM", "McNemar exact (PM majority vs LM CMS4)",
      n_pat = nrow(pp6), p = signif(mc2$p.value, 4),
      detail = sprintf("PM多数CMS4且LM非=%d; 反向=%d", b2, c2))
}

# 3d. GLMM: cms4 ~ organ + (1|patient)
fit_glmm <- function(sub, tag, ref) {
  sub <- copy(sub)
  sub$organ <- relevel(factor(sub$organ), ref = ref)
  f <- tryCatch(
    glmer(cms4N ~ organ + (1 | patient), data = sub, family = binomial,
          control = glmerControl(optimizer = "bobyqa")),
    error = function(e) e)
  if (inherits(f, "error")) {
    cat("GLMM", tag, "失败:", conditionMessage(f), "\n")
    add("glmm", tag, "cms4 ~ organ + (1|patient)",
        n_pat = uniqueN(sub$patient), n_samp = nrow(sub), detail = paste("拟合失败:", conditionMessage(f)))
    return(invisible(NULL))
  }
  sing <- isSingular(f)
  sm <- summary(f)$coefficients
  ci <- suppressMessages(confint(f, method = "Wald"))
  for (term in rownames(sm)[-1]) {
    or <- exp(sm[term, "Estimate"])
    add("glmm", tag, term, n_pat = uniqueN(sub$patient), n_samp = nrow(sub),
        est = round(or, 3), lo = round(exp(ci[term, 1]), 3),
        hi = round(exp(ci[term, 2]), 3), p = signif(sm[term, "Pr(>|z|)"], 4),
        detail = sprintf("OR(CMS4|%s vs %s)=%.2f; singular=%s; 随机截距SD=%.2f",
                         gsub("^organ", "", term), ref, or, sing,
                         sqrt(as.numeric(VarCorr(f)$patient))))
    cat(sprintf("GLMM %s | %s: OR=%.2f [%.2f, %.2f], P=%.4g (singular=%s)\n",
                tag, term, or, exp(ci[term, 1]), exp(ci[term, 2]),
                sm[term, "Pr(>|z|)"], sing))
  }
}
fit_glmm(m[isPM | isPT], "PM_vs_PT", "primary tumor")
fit_glmm(m[isPM | isLM], "PM_vs_LM", "liver metastasis")

# 3e. 阈值版（posterior>=0.5）敏感性：每患者比例 + 配对 Wilcoxon
mt <- m[!is.na(cmsT)]
ppt_pm <- mt[isPM == TRUE, .(prop = mean(cms4T), n = .N), by = patient]
ppt_pt <- mt[isPT == TRUE, .(prop = mean(cms4T), n = .N), by = patient]
ppt <- merge(ppt_pm[, .(patient, prop_PM_T = prop, n_PM_T = n)],
             ppt_pt[, .(patient, prop_PT_T = prop, n_PT_T = n)], by = "patient")
wt3 <- suppressWarnings(wilcox.test(ppt$prop_PM_T, ppt$prop_PT_T, paired = TRUE, exact = FALSE))
cat(sprintf("\n阈值版敏感性: 可分型样本 %d/113; PM vs PT 患者级Wilcoxon P=%.4f (n=%d 患者)\n",
            nrow(mt), wt3$p.value, nrow(ppt)))
add("sensitivity_threshold", "PM_vs_PT", "paired Wilcoxon, posterior>=0.5 typable only",
    n_pat = nrow(ppt), n_samp = nrow(mt),
    est = round(median(ppt$prop_PM_T - ppt$prop_PT_T), 4), p = signif(wt3$p.value, 4),
    detail = sprintf("可分型样本%d/113 (%.0f%%); 患者级PM CMS4比例均值=%.3f vs PT=%.3f",
                     nrow(mt), 100 * nrow(mt) / 113,
                     mean(ppt$prop_PM_T), mean(ppt$prop_PT_T)))

## ---------- (4) scripts/11 患者阻断敏感性: duplicateCorrelation ----------
expr <- fread(file.path(proc, "GSE190609_expr_RPKM.tsv.gz"))
genes <- expr[[1]]
emat <- log2(as.matrix(expr[, -1], rownames = genes) + 0.5)
# merge 会按 title 排序，需把 meta 重排回表达矩阵列序
meta11 <- m[match(colnames(emat), m$title)]
stopifnot(all(!is.na(meta11$patient)), all(meta11$title == colnames(emat)))
meta11$organ <- factor(meta11$organ,
  levels = c("primary tumor", "liver metastasis", "lymph node metastasis",
             "peritoneal metastasis", "ovarium metastasis"))
design11 <- model.matrix(~ 0 + organ + patient, data = meta11)
colnames(design11) <- make.names(colnames(design11))
colnames(design11) <- sub("^organ", "", colnames(design11))
keep <- rowSums(emat >= log2(1.5)) >= 10
emat_f <- emat[keep, ]
cat("\nduplicateCorrelation 输入:", nrow(emat_f), "基因 x", ncol(emat_f), "样本\n")

# (a) 复现 scripts/11: 固定效应阻断
fit_fe <- lmFit(emat_f, design11)
cm <- makeContrasts(PM_vs_PT = peritoneal.metastasis - primary.tumor, levels = design11)
fit_fe <- eBayes(contrasts.fit(fit_fe, cm), trend = TRUE)
tt_fe <- topTable(fit_fe, number = Inf, sort.by = "none")
sig_fe <- tt_fe$adj.P.Val < 0.05 & abs(tt_fe$logFC) >= 1

# (b) duplicateCorrelation 随机效应阻断（limma 官方推荐的重复测量处理）
design_org <- model.matrix(~ 0 + organ, data = meta11)
colnames(design_org) <- make.names(colnames(design_org))
colnames(design_org) <- sub("^organ", "", colnames(design_org))
dc <- duplicateCorrelation(emat_f, design_org, block = meta11$patient)
cat("consensus correlation (patient):", round(dc$consensus, 4), "\n")
fit_dc <- lmFit(emat_f, design_org, block = meta11$patient, correlation = dc$consensus)
cm_org <- makeContrasts(PM_vs_PT = peritoneal.metastasis - primary.tumor, levels = design_org)
fit_dc <- eBayes(contrasts.fit(fit_dc, cm_org), trend = TRUE)
tt_dc <- topTable(fit_dc, number = Inf, sort.by = "none")
sig_dc <- tt_dc$adj.P.Val < 0.05 & abs(tt_dc$logFC) >= 1
sp <- cor(tt_fe$logFC, tt_dc$logFC, method = "spearman")
ov <- sum(sig_fe & sig_dc)
cat(sprintf("PM_vs_PT: 固定效应显著基因 %d; dupCor显著基因 %d; 交集 %d; logFC Spearman=%.4f\n",
            sum(sig_fe), sum(sig_dc), ov, sp))
add("blocking_sensitivity", "PM_vs_PT", "duplicateCorrelation vs fixed-effect patient",
    n_pat = 12, n_samp = 113, est = round(dc$consensus, 4),
    detail = sprintf("患者内consensus correlation=%.3f; 显著基因(FDR<0.05,|logFC|>=1): 固定效应%d vs dupCor%d (交集%d, Jaccard=%.3f); logFC Spearman=%.4f",
                     dc$consensus, sum(sig_fe), sum(sig_dc), ov,
                     ov / sum(sig_fe | sig_dc), sp))

# 已知间质/EMT 标记方向一致性（两种阻断下）
mk <- c("VIM","ZEB1","ZEB2","COL1A1","COL3A1","FAP","ACTA2","SPARC","TGFBI","TAGLN","FN1")
suppressPackageStartupMessages({library(org.Hs.eg.db); library(AnnotationDbi)})
sym <- mapIds(org.Hs.eg.db, keys = rownames(emat_f), column = "SYMBOL",
              keytype = "ENSEMBL", multiVals = "first")
idx <- match(mk, sym[rownames(emat_f)])
idx <- rownames(emat_f)[idx[!is.na(idx)]]
mk_tab <- data.table(symbol = sym[idx],
                     logFC_fixed = tt_fe[idx, "logFC"], logFC_dupcor = tt_dc[idx, "logFC"])
cat("\nEMT/stromal markers 两法 logFC:\n"); print(mk_tab)

## ---------- 落盘 ----------
res <- rbindlist(results, fill = TRUE)
fwrite(res, file.path(outdir, "cms_patient_level.tsv"), sep = "\t")
cat("\n写出:", file.path(outdir, "cms_patient_level.tsv"), "行数:", nrow(res), "\n")
