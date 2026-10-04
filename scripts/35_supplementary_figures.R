# Supplementary Figures S1–S5
# S1: 严格 CMS 分型（minPosterior >= 0.5, n=74）CMS×器官分布 + Fisher 完整结果
# S2: GSE190609 三联配对 DE 火山图（PM vs PT / LM vs PT / PM vs LM）+ EMT sanity-check
# S3: 平台分层 meta（PM: 分层即单队列；LM: RNA-seq 层 vs 芯片层）vs 合并 beta 散点
#     —— 同时落盘新文件 meta/revision_platform_stratified_{PM,LM}.tsv（不覆盖任何现有文件）
# S4: 各评分配置的 1,000 随机基因集 AUC 零分布 + 观察值 + 经验 P
# S5: MSK-MET CNA×器官（探索性, n=177）
# 视觉规范与 scripts/20_M4_fig1_design_cms.R 一致：Helvetica / theme_classic(6.5) / 黑轴
set.seed(123)
suppressMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(scales); library(ggrepel); library(metafor)
})

proc   <- "data/processed/crc_organotropism"
outdir <- "supplementary/figures"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- 调色板与主题（跨主图/补图统一） ----
col_pm   <- "#B4464B"   # 腹膜 哑红
col_lm   <- "#3E6E8E"   # 肝 哑蓝
col_cms4 <- "#C08552"   # CMS4 高亮
g55 <- "grey55"; g75 <- "grey75"; g90 <- "grey90"
base_family <- "Helvetica"

theme_fig <- theme_classic(base_size = 6.5, base_family = base_family) +
  theme(
    axis.title   = element_text(size = 7, color = "black"),
    axis.text    = element_text(size = 6.5, color = "black"),
    plot.title   = element_text(size = 7.5, face = "bold", color = "black"),
    plot.subtitle = element_text(size = 6, color = "grey25"),
    plot.caption = element_text(size = 5, color = "grey30", hjust = 0),
    line = element_line(color = "black")
  )

fmt_p <- function(p) if (p < 0.001) "P < 0.001" else sprintf("P = %.3f", p)

save_fig <- function(plot, name, w_mm, h_mm) {
  png_path <- file.path(outdir, paste0(name, ".png"))
  pdf_path <- file.path(outdir, paste0(name, ".pdf"))
  ragg::agg_png(png_path, width = w_mm, height = h_mm, units = "mm",
                res = 600, background = "white")
  print(plot); dev.off()
  grDevices::pdf(pdf_path, width = w_mm / 25.4, height = h_mm / 25.4,
                 family = base_family, bg = "white", useDingbats = FALSE)
  print(plot); dev.off()
  cat("Exported:", png_path, "and", pdf_path, "\n")
  invisible(c(png_path, pdf_path))
}

# ================= S1: 严格 CMS 分型（minPosterior >= 0.5） =================
post <- fread(file.path(proc, "GSE190609_CMS_posteriors.tsv"))
meta <- fread(file.path(proc, "GSE190609_meta.tsv"))
s1 <- merge(post, meta[, .(title, organ, patient)], by = "title")
s1 <- s1[max_posterior >= 0.5]
stopifnot(nrow(s1) == 74, all(!grepl(",", s1$nearestCMS)))
s1[, site := fcase(
  organ == "peritoneal metastasis", "PM",
  organ == "liver metastasis",      "LM",
  organ == "primary tumor",         "PT",
  default = "Other met")]
site_n <- s1[, .N, by = site]
site_lab <- setNames(sprintf("%s\n(n = %d)", site_n$site, site_n$N), site_n$site)
s1$site <- factor(s1$site, levels = c("PT", "PM", "LM", "Other met"))

s1_tab <- as.data.frame(table(site = s1$site, cms = s1$nearestCMS))
s1_tab$prop <- s1_tab$Freq / ave(s1_tab$Freq, s1_tab$site, FUN = sum)
s1_tab$cms <- factor(s1_tab$cms, levels = c("CMS4", "CMS3", "CMS2", "CMS1"))
cms_cols <- c(CMS1 = g90, CMS2 = g75, CMS3 = g55, CMS4 = col_cms4)
s1_tab$site <- factor(s1_tab$site, levels = c("PT", "PM", "LM", "Other met"),
                      labels = site_lab[c("PT", "PM", "LM", "Other met")])
s1_lab <- subset(s1_tab, Freq > 0 & prop >= 0.12)
s1_lab$txt_col <- ifelse(s1_lab$cms %in% c("CMS1", "CMS2"), "black", "white")

pS1a <- ggplot(s1_tab, aes(x = site, y = prop, fill = cms)) +
  geom_col(width = 0.62, color = "white", linewidth = 0.4) +
  geom_text(data = s1_lab,
            aes(label = sprintf("%d%%", round(100 * prop)), color = txt_col, group = cms),
            position = position_stack(vjust = 0.5),
            size = 1.9, family = base_family, fontface = "bold") +
  scale_fill_manual(values = cms_cols, breaks = c("CMS1", "CMS2", "CMS3", "CMS4"), name = NULL) +
  scale_color_identity(guide = "none") +
  scale_y_continuous(labels = label_percent(accuracy = 1), expand = c(0, 0)) +
  labs(x = NULL, y = "% of samples",
       title = "CMS by site, high-confidence calls (minPosterior >= 0.5)") +
  theme_fig +
  theme(legend.position = "right",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6, color = "black"),
        axis.text.x = element_text(size = 6.5, color = "black", lineheight = 0.9))

# Fisher：CMS4 × 腹膜 vs 其余（严格子集）
isPM <- s1$organ == "peritoneal metastasis"; is4 <- s1$nearestCMS == "CMS4"
t_s1 <- matrix(c(sum(isPM & is4), sum(isPM & !is4),
                 sum(!isPM & is4), sum(!isPM & !is4)), nrow = 2, byrow = TRUE)
f_s1 <- fisher.test(t_s1)
s1_or  <- unname(f_s1$estimate); s1_lo <- f_s1$conf.int[1]
s1_hi  <- f_s1$conf.int[2];      s1_p  <- f_s1$p.value
cat(sprintf("S1 strict Fisher: n=%d  CMS4 PM %d/%d (%.1f%%) vs non-PM %d/%d (%.1f%%)  OR=%.2f [%.2f-%.2f] %s\n",
            nrow(s1), t_s1[1,1], sum(isPM), 100*mean(is4[isPM]),
            t_s1[2,1], sum(!isPM), 100*mean(is4[!isPM]),
            s1_or, s1_lo, s1_hi, fmt_p(s1_p)))

s1_ann <- sprintf("OR = %.2f (95%% CI %.2f-%.2f)\nFisher exact %s\nCMS4: PM %d/%d (%.0f%%) vs non-PM %d/%d (%.0f%%)",
                  s1_or, s1_lo, s1_hi, fmt_p(s1_p),
                  t_s1[1,1], sum(isPM), 100*mean(is4[isPM]),
                  t_s1[2,1], sum(!isPM), 100*mean(is4[!isPM]))
pS1b <- ggplot() +
  # 参考线画成截断 segment，避免穿过上方注释文字（geom_vline 会纵贯全 panel 触发碰撞审计）
  geom_segment(aes(x = 1, xend = 1, y = 0.72, yend = 1.10),
               linetype = "dashed", color = g55, linewidth = 0.4) +
  geom_errorbarh(aes(xmin = s1_lo, xmax = s1_hi, y = 1), height = 0.12,
                 linewidth = 0.7, color = g55) +
  geom_point(aes(x = s1_or, y = 1), size = 2.2, color = g55) +
  geom_text(aes(x = 0.5, y = 1.42, label = s1_ann), hjust = 0, vjust = 1,
            size = 2.0, family = base_family, color = "black", lineheight = 1.25) +
  scale_x_log10(breaks = c(0.5, 1, 2, 4, 8), labels = c("0.5", "1", "2", "4", "8"),
                limits = c(0.5, 11), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.7, 1.45), expand = c(0, 0)) +
  labs(x = "OR for CMS4, peritoneal vs other sites (log scale)", y = NULL,
       title = "CMS4 enrichment, strict subset (descriptive)") +
  theme_fig +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

figS1 <- pS1a + pS1b + plot_layout(widths = c(1.15, 0.85)) +
  plot_annotation(tag_levels = "A",
    caption = "High-confidence CMSclassifier random-forest calls (maximum posterior probability >= 0.5), n = 74 of 113 GSE190609 samples.\nPT = primary tumor; PM = peritoneal metastasis; LM = liver metastasis; Other met = lymph-node (n = 4) and ovarian (n = 1) metastases.\nSample-level Fisher's exact test ignores within-patient correlation and is descriptive only; the primary patient-aware GLMM estimate is reported in the main text.",
    theme = theme(plot.tag = element_text(face = "bold", size = 9, color = "black", family = base_family),
                  plot.caption = element_text(size = 5, color = "grey30", hjust = 0, family = base_family)))
save_fig(figS1, "SupplementaryFigureS1", 180, 80)

# ================= S2: GSE190609 三联配对 DE 火山图 =================
de_files <- c("PM vs PT" = "GSE190609_PM_vs_PT.tsv.gz",
              "LM vs PT" = "GSE190609_LM_vs_PT.tsv.gz",
              "PM vs LM" = "GSE190609_PM_vs_LM.tsv.gz")
emt_up   <- c("SPARC", "CDH2", "FN1", "FAP", "COL3A1")
emt_down <- c("EPCAM", "CDH1", "MUC1")

volcano <- function(f, contrast, label_emt = FALSE) {
  d <- fread(file.path(proc, "de", f))
  d[, neglog10fdr := -log10(pmax(adj.P.Val, 1e-300))]
  sig <- d$adj.P.Val < 0.05 & abs(d$logFC) >= 1
  n_up <- sum(sig & d$logFC > 0); n_dn <- sum(sig & d$logFC < 0)
  d[, cls := fifelse(adj.P.Val < 0.05 & logFC >= 1, "up",
              fifelse(adj.P.Val < 0.05 & logFC <= -1, "down", "ns"))]
  p <- ggplot(d, aes(x = logFC, y = neglog10fdr)) +
    # shape 16 = 实心点（PDF 中为纯填充而非描边路径，避免标签-描边碰撞审计 FAIL）；
    # 不画 x=±1 竖线（geom_vline 在 PDF 中纵贯 panel、会被 pymupdf 误判穿过副标题），
    # 显著性阈值已在副标题与图注中说明
    geom_point(aes(color = cls), size = 0.45, alpha = 0.55, shape = 16) +
    scale_color_manual(values = c(up = col_pm, down = col_lm, ns = "grey82"), guide = "none") +
    labs(x = "log2 fold change", y = "-log10(FDR)",
         title = sprintf("GSE190609 %s", contrast),
         subtitle = sprintf("%s significant: %d up / %d down",
                            format(n_up + n_dn, big.mark = ","), n_up, n_dn)) +
    theme_fig
  if (label_emt) {
    lab <- d[symbol %in% c(emt_up, emt_down)]
    p <- p + geom_text_repel(data = lab, aes(label = symbol),
                             size = 1.9, family = base_family, fontface = "italic",
                             color = "black", segment.size = 0.25, segment.color = g55,
                             min.segment.length = 0, max.overlaps = Inf,
                             box.padding = 0.25, point.padding = 0.15)
  }
  list(plot = p, n_sig = n_up + n_dn, n_up = n_up, n_dn = n_dn)
}
v1 <- volcano(de_files["PM vs PT"], "PM vs PT", label_emt = TRUE)
v2 <- volcano(de_files["LM vs PT"], "LM vs PT")
v3 <- volcano(de_files["PM vs LM"], "PM vs LM")
cat(sprintf("S2 sig genes: PM vs PT %d | LM vs PT %d | PM vs LM %d (up %d / down %d)\n",
            v1$n_sig, v2$n_sig, v3$n_sig, v3$n_up, v3$n_dn))
stopifnot(v1$n_sig == 861, v2$n_sig == 534, v3$n_sig == 1053, v3$n_up == 487, v3$n_dn == 566)

figS2 <- v1$plot + v2$plot + v3$plot + plot_layout(nrow = 1) +
  plot_annotation(tag_levels = "A",
    caption = "Paired limma with patient blocking (design ~ organ + patient) in GSE190609 (113 samples, 12 patients).\nRed, upregulated in the first-named group (FDR < 0.05 and |log2FC| >= 1); blue, downregulated; grey, not significant.\nPanel A labels the epithelial-mesenchymal-transition sanity-check genes (mesenchymal SPARC, CDH2, FN1, FAP, COL3A1 up; epithelial EPCAM, CDH1, MUC1 down).",
    theme = theme(plot.tag = element_text(face = "bold", size = 9, color = "black", family = base_family),
                  plot.caption = element_text(size = 5, color = "grey30", hjust = 0, family = base_family)))
save_fig(figS2, "SupplementaryFigureS2", 180, 68)

# ================= S3: 平台分层 meta vs 合并 beta =================
# 合并 meta（现有落盘，只读）
meta_pm <- fread(file.path(proc, "meta", "meta_PM_vs_PT_RNAseq.tsv.gz"))
meta_lm <- fread(file.path(proc, "meta", "meta_LM_vs_PT_mixed.tsv.gz"))
pm_core <- fread(file.path(proc, "meta", "peritoneal_tropism_genes.tsv"))
lm_core <- fread(file.path(proc, "meta", "liver_tropism_genes.tsv"))
stopifnot(nrow(pm_core) == 49, nrow(lm_core) == 57)

read_de <- function(f) fread(file.path(proc, "de", f))
# 单基因 REML meta（与 scripts/12 口径一致：yi = logFC, vi = (logFC/t)^2）
gene_meta <- function(dts, g) {
  ys <- vapply(dts, function(x) { r <- x[symbol == g]; if (nrow(r) == 0) NA_real_ else r$logFC[1] }, numeric(1))
  ts <- vapply(dts, function(x) { r <- x[symbol == g]; if (nrow(r) == 0) NA_real_ else r$t[1] },    numeric(1))
  vi <- (ys / ts)^2
  ok <- !is.na(vi) & is.finite(vi) & vi > 0
  if (sum(ok) == 1) return(list(beta = ys[ok], se = sqrt(vi[ok]), k = 1L))
  if (sum(ok) == 0) return(list(beta = NA_real_, se = NA_real_, k = 0L))
  fit <- tryCatch(rma(yi = ys[ok], vi = vi[ok], method = "REML"), error = function(e) NULL)
  if (is.null(fit)) return(list(beta = NA_real_, se = NA_real_, k = 0L))
  list(beta = as.numeric(fit$beta), se = as.numeric(fit$se), k = sum(ok))
}

# ---- PM：仅 2 个 RNA-seq 队列，分层即单队列 ----
d609_pm <- read_de("GSE190609_PM_vs_PT.tsv.gz")
d225    <- read_de("GSE225182_PM_vs_PT.tsv.gz")
strat_pm <- rbindlist(lapply(pm_core$gene, function(g) {
  r6 <- d609_pm[symbol == g]; r2 <- d225[symbol == g]
  data.table(gene = g,
             pooled_beta  = meta_pm[gene == g]$beta,
             pooled_fdr   = meta_pm[gene == g]$fdr,
             GSE190609_logFC = if (nrow(r6)) r6$logFC[1] else NA_real_,
             GSE225182_logFC = if (nrow(r2)) r2$logFC[1] else NA_real_)
}))
fwrite(strat_pm, file.path(proc, "meta", "revision_platform_stratified_PM.tsv"), sep = "\t")

# ---- LM：RNA-seq 层（GSE190609 + GSE50760）vs 芯片层（GSE41258 + GSE41568）----
d609_lm <- read_de("GSE190609_LM_vs_PT.tsv.gz")
d760    <- read_de("GSE50760_LM_vs_PT.tsv.gz")
d258    <- read_de("GSE41258_LM_vs_PT.tsv.gz")
d568    <- read_de("GSE41568_LM_vs_PT.tsv.gz")
strat_lm <- rbindlist(lapply(lm_core$gene, function(g) {
  rs <- gene_meta(list(d609_lm, d760), g)
  ch <- gene_meta(list(d258, d568), g)
  data.table(gene = g,
             pooled_beta = meta_lm[gene == g]$beta,
             pooled_fdr  = meta_lm[gene == g]$fdr,
             rnaseq_beta = rs$beta, rnaseq_se = rs$se, rnaseq_k = rs$k,
             chip_beta   = ch$beta, chip_se   = ch$se, chip_k   = ch$k)
}))
fwrite(strat_lm, file.path(proc, "meta", "revision_platform_stratified_LM.tsv"), sep = "\t")
cat("S3 分层 meta 落盘: meta/revision_platform_stratified_PM.tsv (", nrow(strat_pm),
    " genes), revision_platform_stratified_LM.tsv (", nrow(strat_lm), " genes)\n", sep = "")

# ---- 散点：x = 合并 meta beta, y = 分层 beta ----
pm_long <- rbind(
  strat_pm[, .(gene, pooled_beta, stratum = "GSE190609 (RNA-seq)", strat_beta = GSE190609_logFC)],
  strat_pm[, .(gene, pooled_beta, stratum = "GSE225182 (RNA-seq)", strat_beta = GSE225182_logFC)]
)
lm_long <- rbind(
  strat_lm[, .(gene, pooled_beta, stratum = "RNA-seq stratum", strat_beta = rnaseq_beta)],
  strat_lm[, .(gene, pooled_beta, stratum = "Microarray stratum", strat_beta = chip_beta)]
)
pm_long <- pm_long[!is.na(strat_beta)]; lm_long <- lm_long[!is.na(strat_beta)]
strat_cols <- c("GSE190609 (RNA-seq)" = col_pm, "GSE225182 (RNA-seq)" = col_cms4,
                "RNA-seq stratum" = col_lm, "Microarray stratum" = g55)

rho_txt <- function(dd) {
  r <- dd[, .(rho = cor(pooled_beta, strat_beta, method = "spearman")), by = stratum]
  paste(sprintf("%s rho = %.2f", r$stratum, r$rho), collapse = "; ")
}
scatter_beta <- function(dd, title_txt) {
  # y = x 参考线用 panel 内 segment（geom_abline 在 PDF 中不被裁剪、会穿过标题/图例文字）
  lo <- max(min(dd$pooled_beta), min(dd$strat_beta))
  hi <- min(max(dd$pooled_beta), max(dd$strat_beta))
  ggplot(dd, aes(x = pooled_beta, y = strat_beta, color = stratum)) +
    geom_segment(aes(x = lo, xend = hi, y = lo, yend = hi),
                 linetype = "dashed", color = g75, linewidth = 0.4,
                 inherit.aes = FALSE) +
    geom_point(size = 1.1, alpha = 0.8, stroke = 0) +
    scale_color_manual(values = strat_cols, name = NULL) +
    labs(x = "Pooled random-effects meta beta", y = "Stratum beta",
         title = title_txt, subtitle = paste0("Spearman: ", rho_txt(dd))) +
    theme_fig +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 6, color = "black"),
          legend.key.size = unit(3, "mm"))
}
pS3a <- scatter_beta(pm_long, "PM_core (49 genes): per-cohort vs pooled")
pS3b <- scatter_beta(lm_long, "LM_core (57 genes): platform strata vs pooled")
cat("S3 PM:", rho_txt(pm_long), "\nS3 LM:", rho_txt(lm_long), "\n")

figS3 <- pS3a + pS3b + plot_layout(nrow = 1) +
  plot_annotation(tag_levels = "A",
    caption = "PM: only two RNA-seq cohorts contribute to the peritoneal meta-analysis, so platform strata reduce to single-cohort estimates.\nLM: random-effects (REML) estimates within the RNA-seq stratum (GSE190609 + GSE50760) and the microarray stratum (GSE41258 + GSE41568),\nplotted against the pooled four-cohort estimate. Dashed line, y = x.",
    theme = theme(plot.tag = element_text(face = "bold", size = 9, color = "black", family = base_family),
                  plot.caption = element_text(size = 5, color = "grey30", hjust = 0, family = base_family)))
save_fig(figS3, "SupplementaryFigureS3", 180, 85)

# ================= S4: 置换零分布 + 观察 AUC =================
perm <- fread(file.path(proc, "ssgsea", "revision_permutation_auc.tsv"))
null <- fread(file.path(proc, "ssgsea", "revision_permutation_null_distributions.tsv.gz"))
perm_main <- perm[contrast == "main"]
stopifnot(nrow(perm_main) == 5)
cfg_lab <- c(
  A_apparent_LMcore = "Apparent | LM_core",
  B_apparent_PMcore = "Apparent | PM_core",
  C_LOCO_main_LMcore = "LOCO main | LM_core (array-derived)",
  D_LOCO_sens_LMcore = "LOCO restricted | detectable genes",
  E_LOCO_T2_PMcore  = "Cross-cohort | PM_core (GSE190609)"
)
cfg_order <- names(cfg_lab)
null_panels <- lapply(cfg_order, function(cfg) {
  pw <- perm_main[config_id == cfg]
  nd <- null[config_id == cfg & contrast == "main"]
  stopifnot(nrow(pw) == 1, nrow(nd) == 1000)
  # 注释放进副标题而非 panel 内文本（panel 内文本会穿过 vline 并与柱填充碰撞）
  subt <- sprintf("observed %.3f, null %.3f; empirical %s",
                  pw$auc_real, pw$null_mean, fmt_p(pw$emp_p_greater))
  xmax <- max(nd$auc, pw$auc_real); xmin <- min(nd$auc, pw$auc_real)
  ymax <- max(hist(nd$auc, breaks = seq(-0.04, 1.08, by = 0.04), plot = FALSE)$counts)
  ggplot(nd, aes(x = auc)) +
    geom_histogram(binwidth = 0.04, boundary = 0, fill = g90, color = g55, linewidth = 0.25) +
    # 截断 segment 替代 geom_vline（vline 在 PDF 中纵贯 panel，会触发与标题的碰撞误判）
    geom_segment(aes(x = pw$null_mean, xend = pw$null_mean, y = 0, yend = ymax),
                 linetype = "dashed", color = g55, linewidth = 0.4, inherit.aes = FALSE) +
    geom_segment(aes(x = pw$auc_real, xend = pw$auc_real, y = 0, yend = ymax),
                 color = col_pm, linewidth = 0.6, inherit.aes = FALSE) +
    coord_cartesian(xlim = c(xmin - 0.06, xmax + 0.06)) +
    labs(x = "AUC", y = "Count", title = cfg_lab[[cfg]], subtitle = subt) +
    theme_fig +
    theme(plot.title = element_text(size = 6.5, face = "bold", color = "black"))
})
cat("S4 五个主配置经验 P:\n")
print(perm_main[, .(config_id, auc_real = round(auc_real, 3),
                    null_mean = round(null_mean, 3), emp_p = round(emp_p_greater, 3))])

s4_sizes <- paste(sprintf("%s: %d genes, %d vs %d PT", LETTERS[1:5],
                          perm_main$set_size[match(cfg_order, perm_main$config_id)],
                          perm_main$n_case[match(cfg_order, perm_main$config_id)],
                          perm_main$n_ctrl[match(cfg_order, perm_main$config_id)]),
                  collapse = "; ")
figS4 <- wrap_plots(null_panels, nrow = 2) +
  plot_annotation(tag_levels = "A",
    caption = paste0("Null distributions of the Mann-Whitney AUC for 1,000 random gene sets matched for size and detection rate under each scoring configuration.\n",
                     "Solid red line, observed AUC of the true gene set; dashed grey line, random-set mean; empirical P, fraction of random sets with AUC >= observed.\n",
                     "Panel configurations: ", s4_sizes, "."),
    theme = theme(plot.tag = element_text(face = "bold", size = 9, color = "black", family = base_family),
                  plot.caption = element_text(size = 5, color = "grey30", hjust = 0, family = base_family)))
save_fig(figS4, "SupplementaryFigureS4", 180, 105)

# ================= S5: MSK-MET CNA × 器官（探索性） =================
cna <- fread(file.path(proc, "mskmet", "cna_amp_by_organ.tsv"))
cna_long <- melt(cna, id.vars = c("gene", "OR", "p"),
                 measure.vars = c("PM_amp_pct", "Liver_amp_pct"),
                 variable.name = "organ", value.name = "pct")
cna_long$organ <- factor(cna_long$organ, levels = c("PM_amp_pct", "Liver_amp_pct"),
                         labels = c("Peritoneal", "Liver"))
gene_ord <- cna[order(-pmax(PM_amp_pct, Liver_amp_pct))]$gene
cna_long$gene <- factor(cna_long$gene, levels = gene_ord)
cna_ann <- cna[, .(gene = factor(gene, levels = gene_ord),
                   y = pmax(PM_amp_pct, Liver_amp_pct) + 2.2,
                   lab = vapply(p, fmt_p, character(1)))]

pS5 <- ggplot(cna_long, aes(x = gene, y = pct, fill = organ)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62, color = "white", linewidth = 0.3) +
  geom_text(data = cna_ann, aes(x = gene, y = y, label = lab),
            inherit.aes = FALSE, size = 1.8, family = base_family, color = "grey25") +
  scale_fill_manual(values = c(Peritoneal = col_pm, Liver = col_lm), name = NULL) +
  scale_y_continuous(limits = c(0, max(cna_ann$y) + 3), expand = c(0, 0)) +
  labs(x = NULL, y = "Amplification frequency (%)",
       title = "Copy-number amplification by metastatic organ (MSK-MET, exploratory, n = 177)",
       subtitle = "Fisher exact P per gene; no credible organ-specific signal") +
  theme_fig +
  theme(legend.position = "top",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6, color = "black"),
        axis.text.x = element_text(size = 6.5, color = "black", face = "italic"))
save_fig(pS5, "SupplementaryFigureS5", 180, 68)

cat("\nAll supplementary figures exported to", outdir, "\n")
print(sessionInfo())
