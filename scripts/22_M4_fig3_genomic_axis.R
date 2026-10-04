# Figure 3 — MSK-MET 基因组轴与生存（腹膜 vs 肝转移灶突变、患者级富集、Cox OS）
# 输出: main_figures/Figure3.png (600 dpi) + Figure3.pdf (cairo 矢量)
# 数据: data/processed/crc_organotropism/mskmet/{mutation_by_organ, patient_mutation_peritoneal,
#       cox_os_with_burden, revision_cox_samesubset}.tsv；Panel B 的 APC 患者级 Fisher 从原始 JSON 重算补齐。
# Panel A 频率为 pooled peritoneal（Intra-Abdominal + Ovary）转移灶样本级重算，与 pooled Fisher 检验口径一致；
# Panel C 为同子集（n=405）含/不含负荷项两模型对照。
set.seed(123)
suppressMessages({
  library(data.table); library(jsonlite)
  library(ggplot2); library(patchwork); library(scales)
})

raw  <- "data/raw/crc_organotropism"
proc <- "data/processed/crc_organotropism/mskmet"
dir.create("main_figures", showWarnings = FALSE)

## ---------- 全局规范 ----------
COL_PM    <- "#B4464B"   # 腹膜 哑红
COL_LIVER <- "#3E6E8E"   # 肝 哑蓝
COL_NS    <- "grey55"
COL_LINE  <- "grey75"
COL_FAINT <- "grey90"
FF <- "Helvetica"

fmt_p <- function(p) ifelse(p < 0.001, "P < 0.001", sprintf("P = %.3f", p))

theme_fig <- theme_classic(base_size = 6.5, base_family = FF) +
  theme(
    axis.text        = element_text(color = "black", size = 6),
    axis.title       = element_text(color = "black", size = 6.5),
    plot.title       = element_text(color = "black", size = 7, face = "bold",
                                    margin = margin(b = 2)),
    plot.subtitle    = element_text(color = "grey30", size = 5.5,
                                    margin = margin(b = 3)),
    plot.caption     = element_text(color = "black", size = 5.5, hjust = 0,
                                    margin = margin(t = 3)),
    legend.text      = element_text(color = "black", size = 6),
    legend.title     = element_blank(),
    legend.key.size  = unit(3, "mm"),
    legend.margin    = margin(0, 0, 0, 0),
    axis.line        = element_line(color = "black", linewidth = 0.3),
    axis.ticks       = element_line(color = "black", linewidth = 0.3),
    plot.margin      = margin(4, 5, 3, 4)
  )

## =====================================================================
## 数据读取与核验
## =====================================================================
mut_organ <- fread(file.path(proc, "mutation_by_organ.tsv"))
stopifnot(all(c("gene","liver_pct","intraab_pct","OR_PM_vs_Liver","p_PM_vs_Liver","fdr") %in%
              colnames(mut_organ)))

pat_mut_file <- fread(file.path(proc, "patient_mutation_peritoneal.tsv"))
stopifnot(all(c("gene","PM_mut","PM_n","nonPM_mut","nonPM_n","OR","p") %in%
              colnames(pat_mut_file)))

cox1 <- fread(file.path(proc, "cox_os_with_burden.tsv"))
stopifnot(all(c("term","exp(coef)","lower .95","upper .95","Pr(>|z|)") %in% colnames(cox1)))

cox_ss <- fread(file.path(proc, "revision_cox_samesubset.tsv"))
stopifnot(all(c("model","term","HR","p","HR_lo95","HR_hi95","n","events") %in% colnames(cox_ss)))

# 同源核验：M1（子集负荷校正）须与 cox_os_with_burden.tsv 逐行一致（HR/p 容差 1e-3）
cat("== 核验 Panel C：revision_cox_samesubset M1 vs cox_os_with_burden ==\n")
chk_m1 <- merge(cox_ss[model == "M1_subset_with_met_count", .(term, HR_new = HR, p_new = p)],
                cox1[, .(term, HR_file = `exp(coef)`, p_file = `Pr(>|z|)`)], by = "term")
print(chk_m1)
stopifnot(nrow(chk_m1) == nrow(cox1),
          all(abs(chk_m1$HR_new - chk_m1$HR_file) < 1e-3),
          all(abs(chk_m1$p_new - chk_m1$p_file) < 1e-3))

## ---------- Panel B: 患者级 Fisher 从原始数据重算（补齐 APC 并复核已有 8 基因）----------
crc_ids <- fromJSON(file.path(raw, "mskmet_crc_sample_ids.json"), simplifyVector = TRUE)
read_attr <- function(f) {
  x <- fromJSON(file.path(raw, f), simplifyVector = FALSE)
  rbindlist(lapply(x, function(r) data.table(sampleId = r$sampleId, patientId = r$patientId,
                                             value = r$value)))
}
site <- read_attr("mskmet_S_METASTATIC_SITE.json")
stype <- read_attr("mskmet_S_SAMPLE_TYPE.json")
setnames(site, "value", "met_site"); setnames(stype, "value", "sample_type")
clin <- merge(site, stype[, .(sampleId, sample_type)], by = "sampleId", all.x = TRUE)
clin <- clin[sampleId %in% crc_ids]

mut <- fromJSON(file.path(raw, "mskmet_crc_mutations.json"), simplifyVector = TRUE)
mut <- data.table(sampleId = mut$sampleId, patientId = mut$patientId,
                  gene = mut$gene$hugoGeneSymbol, mutType = mut$mutationType)
mut <- mut[mutType != "Silent"]
mut_p <- unique(mut[, .(patientId, gene)])

pw <- clin[sample_type == "Metastasis",
           .(hasPM = any(met_site %in% c("Intra-Abdominal", "Ovary"))), by = patientId]
cat("患者级队列: n =", nrow(pw), ", 腹膜受累 =", sum(pw$hasPM), "\n")

genesB <- c("RNF43","SMAD4","PTEN","KRAS","BRAF","NRAS","TP53","PIK3CA","APC")
resB <- rbindlist(lapply(genesB, function(g) {
  ids <- unique(mut_p[gene == g]$patientId)
  a <- sum(pw$hasPM & pw$patientId %in% ids);  b <- sum(pw$hasPM & !pw$patientId %in% ids)
  c1 <- sum(!pw$hasPM & pw$patientId %in% ids); d <- sum(!pw$hasPM & !pw$patientId %in% ids)
  ft <- fisher.test(matrix(c(a, b, c1, d), nrow = 2, byrow = TRUE))
  data.table(gene = g, PM_mut = a, PM_n = a + b, nonPM_mut = c1, nonPM_n = c1 + d,
             pct_PM = 100 * a / (a + b), pct_nonPM = 100 * c1 / (c1 + d),
             OR = as.numeric(ft$estimate), lo = ft$conf.int[1], hi = ft$conf.int[2],
             p = ft$p.value)
}))
resB[, fdr := p.adjust(p, "BH")]   # 9 基因集合统一 BH
resB[, `:=`(pct_PM = round(pct_PM, 1), pct_nonPM = round(pct_nonPM, 1))]

cat("\n== 核验 Panel B：重算 vs 文件（OR / p）==\n")
cmp <- merge(resB[, .(gene, OR_new = round(OR, 2), p_new = signif(p, 3))],
             pat_mut_file[, .(gene, OR_file = OR, p_file = p)], by = "gene", all = TRUE)
print(cmp)
cat("Panel B 重算全表（含 CI、9 基因 FDR）:\n")
print(resB[, .(gene, pct_PM, pct_nonPM, OR = round(OR, 2), lo = round(lo, 2),
               hi = round(hi, 2), p = signif(p, 3), fdr = signif(fdr, 3))])

cat("\n== 核验 Panel C：同子集两模型关键数值 ==\n")
print(cox_ss[model %in% c("M0_subset_no_burden","M1_subset_with_met_count"),
             .(model, term, HR = round(HR, 2), lo = round(HR_lo95, 2),
               hi = round(HR_hi95, 2), P = signif(p, 3), n, events)])

## =====================================================================
## Panel A — 转移灶突变频率哑铃图（腹膜 vs 肝）
## =====================================================================
genesA <- c("APC","TP53","KRAS","SMAD4","PIK3CA","BRAF","RNF43")

## 频率口径与检验一致：pooled peritoneal = Intra-Abdominal + Ovary 转移灶样本，从原始 JSON 重算
met_grp <- clin[sample_type == "Metastasis" & met_site %in% c("Intra-Abdominal","Ovary","Liver"),
                .(sampleId, grp = fifelse(met_site == "Liver", "liver", "peritoneal"))]
n_pm <- uniqueN(met_grp[grp == "peritoneal"]$sampleId)
n_lv <- uniqueN(met_grp[grp == "liver"]$sampleId)
freqA <- rbindlist(lapply(genesA, function(g) {
  mutg <- unique(mut[gene == g]$sampleId)
  data.table(gene = g,
             pm_n    = uniqueN(met_grp[grp == "peritoneal" & sampleId %in% mutg]$sampleId),
             liver_n = uniqueN(met_grp[grp == "liver"     & sampleId %in% mutg]$sampleId))
}))
freqA[, `:=`(pm = 100 * pm_n / n_pm, liver = 100 * liver_n / n_lv)]

cat("== 核验 Panel A：pooled peritoneal vs liver 转移灶样本级频率 ==\n")
cat("pooled peritoneal n =", n_pm, ", liver n =", n_lv, "\n")
print(freqA[, .(gene, pm_n, pm = round(pm, 1), liver_n, liver = round(liver, 1))])
print(mut_organ[gene %in% genesA, .(gene, OR_PM_vs_Liver, p_PM_vs_Liver, fdr)])  # OR/P/FDR 取自 pooled 检验
# 预期：n_pm = 124（81 腹腔 + 43 卵巢），n_lv = 624；APC 64.5 vs 82.2（80/124、513/624），fdr=5.148e-4
stopifnot(n_pm == 124, n_lv == 624,
          freqA[gene == "APC"]$pm_n == 80, freqA[gene == "APC"]$liver_n == 513,
          abs(freqA[gene == "APC"]$pm - 64.5) < 0.1,
          abs(freqA[gene == "APC"]$liver - 82.2) < 0.1)

dA <- merge(mut_organ[gene %in% genesA, .(gene, fdr, p = p_PM_vs_Liver)],
            freqA[, .(gene, pm, liver)], by = "gene")
dA[, ord := pmax(pm, liver)]
setorder(dA, -ord, -pm)
dA[, gene := factor(gene, levels = rev(dA$gene))]

dA_long <- melt(dA, id.vars = c("gene","fdr"), measure.vars = c("pm","liver"),
                variable.name = "site", value.name = "pct")
dA_long[, site := factor(ifelse(site == "pm", "Peritoneal", "Liver"),
                         levels = c("Peritoneal","Liver"))]
# 数值标签放在远离另一侧圆点的一端，避免碰撞；RNF43 肝值太小，标签改放圆点上方
dA[, `:=`(pm_lab_x    = ifelse(pm <= liver, pm - 2.2, pm + 2.2),
          liver_lab_x = ifelse(pm <= liver, liver + 2.2, liver - 2.2),
          pm_hj       = ifelse(pm <= liver, 1, 0),
          liver_hj    = ifelse(pm <= liver, 0, 1),
          pm_dy       = 0, liver_dy = 0)]
dA[, ynum := as.numeric(gene)]
dA[gene == "RNF43", `:=`(liver_lab_x = liver, liver_hj = 0.5, liver_dy = 0.45)]

star_lab <- dA[fdr < 0.05, .(gene, ynum, x = pmax(pm, liver) + 9.5,
                             lab = ifelse(fdr < 0.001, "***", ifelse(fdr < 0.01, "**", "*")))]

pA <- ggplot(dA) +
  geom_segment(aes(x = pm, xend = liver, y = gene, yend = gene),
               color = COL_LINE, linewidth = 0.7) +
  geom_point(data = dA_long,
             aes(x = pct, y = gene, color = site), size = 2.1) +
  geom_text(aes(x = pm_lab_x, y = ynum + pm_dy, label = sprintf("%.1f", pm), hjust = pm_hj),
            color = COL_PM, size = 5 / .pt, family = FF) +
  geom_text(aes(x = liver_lab_x, y = ynum + liver_dy, label = sprintf("%.1f", liver), hjust = liver_hj),
            color = COL_LIVER, size = 5 / .pt, family = FF) +
  geom_text(data = star_lab, aes(x = x, y = ynum, label = lab),
            color = "black", size = 6 / .pt, family = FF, fontface = "bold", hjust = 0) +
  scale_color_manual(values = c(Peritoneal = COL_PM, Liver = COL_LIVER)) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(title = "Mutation frequency in metastatic sites",
       subtitle = "Non-silent driver mutations, MSK-MET CRC metastasis samples",
       x = "Mutated samples (%)", y = NULL, color = NULL,
       caption = "*** FDR < 0.001 (Fisher, peritoneal vs liver, BH).\nPeritoneal = intra-abdominal + ovary pooled (n = 124);\nFisher test: same pooled group vs liver (n = 624).") +
  theme_fig +
  theme(axis.text.y = element_text(face = "italic", size = 6.5, color = "black"),
        legend.position = c(0.78, 0.16),
        legend.background = element_rect(fill = "white", color = NA),
        plot.caption = element_text(color = "grey30", size = 5, hjust = 0,
                                    margin = margin(t = 3)))

## =====================================================================
## Panel B — 患者级 Fisher 富集森林图
## =====================================================================
dB <- copy(resB)
setorder(dB, OR)
dB[, gene := factor(gene, levels = dB$gene)]          # OR 升序 → 大 OR 在顶部
dB[, sig_col := fifelse(p < 0.05 & OR > 1, COL_PM,
                        fifelse(p < 0.05 & OR < 1, COL_LIVER, COL_NS))]
dB[, lab_or := sprintf("%.2f (%.2f\u2013%.2f)", OR, lo, hi)]
dB[, lab_p  := fmt_p(p)]
dB[fdr < 0.01, lab_p := paste0(lab_p, " **")]
x_txt1 <- 5.2; x_txt2 <- 26

pB <- ggplot(dB, aes(y = gene)) +
  annotate("segment", x = 1, xend = 1, y = 0.5, yend = length(genesB) + 0.5,
           color = COL_NS, linetype = "dashed", linewidth = 0.35) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.22, linewidth = 0.45,
                 color = dB$sig_col) +
  geom_point(aes(x = OR), color = dB$sig_col, size = 1.9) +
  geom_text(aes(x = x_txt1, label = lab_or), hjust = 0,
            size = 5 / .pt, family = FF, color = "black") +
  geom_text(aes(x = x_txt2, label = lab_p), hjust = 0,
            size = 5 / .pt, family = FF, color = "black") +
  annotate("text", x = x_txt1, y = length(genesB) + 0.75, label = "OR (95% CI)",
           hjust = 0, size = 5.5 / .pt, family = FF, fontface = "bold", color = "black") +
  annotate("text", x = x_txt2, y = length(genesB) + 0.75, label = "P value",
           hjust = 0, size = 5.5 / .pt, family = FF, fontface = "bold", color = "black") +
  scale_x_log10(limits = c(0.25, 95), breaks = c(0.25, 0.5, 1, 2, 4),
                labels = c("0.25","0.5","1","2","4"),
                expand = expansion(mult = c(0.02, 0))) +
  scale_y_discrete(expand = expansion(add = c(0.7, 1.3))) +
  coord_cartesian(clip = "off") +
  labs(title = "Patient-level enrichment, peritoneal involvement",
       subtitle = "Fisher test: peritoneal-involved (n = 124) vs non-involved (n = 1,023)",
       x = "Odds ratio (log scale)", y = NULL,
       caption = "** FDR < 0.01 (BH across 9 genes). Red: significant enrichment with peritoneal\ninvolvement; blue: significant depletion; grey: not significant.") +
  theme_fig +
  theme(axis.text.y = element_text(face = "italic", size = 6.5, color = "black"),
        plot.caption = element_text(color = "grey30", size = 5, hjust = 0,
                                    margin = margin(t = 3)))

## =====================================================================
## Panel C — Cox OS 森林图（同子集：含 vs 不含负荷项）
## =====================================================================
term_lab <- c(hasPMTRUE    = "Peritoneal~involvement",
              hasLiverTRUE = "Liver~involvement",
              hasLungTRUE  = "Lung~involvement",
              BRAF_mutTRUE = "italic(BRAF)~mutation",
              age          = "Age~(per~year)",
              met_count    = "No.~of~metastatic~sites")
term_order <- c("hasPMTRUE","hasLiverTRUE","hasLungTRUE","BRAF_mutTRUE","age","met_count")

mk_cox <- function(dt, model) {
  dt[, .(term, HR, lo = HR_lo95, hi = HR_hi95, p)][, model := model]
}
dC <- rbind(mk_cox(cox_ss[model == "M1_subset_with_met_count"], "adj"),
            mk_cox(cox_ss[model == "M0_subset_no_burden"],     "noburden"))
dC[, term := factor(term, levels = rev(term_order))]
dC[, y := as.numeric(term) + ifelse(model == "adj", 0.16, -0.16)]
dC[, lab := sprintf("%.2f (%.2f\u2013%.2f); %s", HR, lo, hi, fmt_p(p))]
x_c1 <- 5.0; x_c2 <- 38
n_terms <- length(term_order)

pC <- ggplot(dC) +
  annotate("segment", x = 1, xend = 1, y = 0.5, yend = n_terms + 0.5,
           color = COL_NS, linetype = "dashed", linewidth = 0.35) +
  geom_errorbarh(data = dC[model == "adj"],
                 aes(xmin = lo, xmax = hi, y = y), height = 0.18,
                 linewidth = 0.45, color = "black") +
  geom_point(data = dC[model == "adj"], aes(x = HR, y = y),
             shape = 16, size = 1.9, color = "black") +
  geom_errorbarh(data = dC[model == "noburden"],
                 aes(xmin = lo, xmax = hi, y = y), height = 0.18,
                 linewidth = 0.45, color = COL_NS) +
  geom_point(data = dC[model == "noburden"], aes(x = HR, y = y),
             shape = 1, size = 1.9, color = COL_NS, stroke = 0.7) +
  geom_text(data = dC[model == "adj"],  aes(x = x_c1, y = as.numeric(term), label = lab),
            hjust = 0, size = 5 / .pt, family = FF, color = "black") +
  geom_text(data = dC[model == "noburden"], aes(x = x_c2, y = as.numeric(term), label = lab),
            hjust = 0, size = 5 / .pt, family = FF, color = COL_NS) +
  annotate("text", x = x_c1, y = n_terms + 0.95,
           label = "Burden-adjusted\nHR (95% CI), P",
           hjust = 0, size = 5.5 / .pt, family = FF, fontface = "bold", color = "black") +
  annotate("text", x = x_c2, y = n_terms + 0.95,
           label = "Without burden term\nHR (95% CI), P",
           hjust = 0, size = 5.5 / .pt, family = FF, fontface = "bold", color = "grey35") +
  scale_y_continuous(breaks = seq_len(n_terms),
                     labels = parse(text = unname(term_lab[rev(term_order)])),
                     limits = c(0.5, n_terms + 1.5),
                     expand = expansion(mult = c(0, 0))) +
  scale_x_log10(limits = c(0.35, 150), breaks = c(0.5, 1, 2, 3),
                labels = c("0.5","1","2","3"),
                expand = expansion(mult = c(0.02, 0))) +
  labs(title = "Overall survival, Cox proportional-hazards models",
       subtitle = "Filled black: burden-adjusted; open grey: without burden term — same subset n = 405 (226 events)",
       x = "Hazard ratio (log scale)", y = NULL,
       caption = "In the same subset, peritoneal involvement is not significant even without the burden term\n(HR = 1.38, 95% CI 0.86–2.22, P = 0.181); with burden, HR = 1.49 (0.92–2.39, P = 0.102), so the loss of\nsignificance mainly reflects subset power rather than burden confounding. Burden: HR = 1.08 per site, P < 0.001.\nFull-cohort model without burden term (n = 1,136): peritoneal involvement HR = 1.42 (1.05–1.92), P = 0.022.") +
  theme_fig

## =====================================================================
## 组装 + 导出
## =====================================================================
W_MM <- 180; H_MM <- 152
# 扁平 design（避免嵌套 patchwork），保证 panel_alignment.R 能识别三个 panel
fig <- pA + pB + pC +
  plot_layout(design = "12\n33", widths = c(0.92, 1.08), heights = c(1.22, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 8, color = "black"),
        plot.tag.position = c(0.005, 0.995))

png_path <- "main_figures/Figure3.png"
pdf_path <- "main_figures/Figure3.pdf"
ggsave(png_path, fig, device = ragg::agg_png, width = W_MM, height = H_MM,
       units = "mm", dpi = 600, bg = "white")
# 本机 cairo_pdf 依赖 /opt/X11 缺失不可用；改用 grDevices::pdf（矢量、文字可编辑）
ggsave(pdf_path, fig, device = grDevices::pdf, width = W_MM / 25.4, height = H_MM / 25.4,
       family = FF, useDingbats = FALSE)
cat("written:", png_path, "and", pdf_path, "\n")

## ---------- Panel 对齐 QA manifest ----------
## panel_alignment.R 的 convertWidth 在本机 grid 版本上把 null 单位转成 0，
## 这里按同一 schema 手动解析 gtable 尺寸（绝对单位转 pt + null 单位按比例分配剩余空间）
resolve_gtable_units <- function(u, total_pt, what = c("width", "height")) {
  what <- match.arg(what)
  conv <- if (what == "width") grid::convertWidth else grid::convertHeight
  n <- length(u)
  pt <- numeric(n); is_null <- logical(n)
  for (i in seq_len(n)) {
    one <- u[i]
    if (grid::unitType(one) == "null") {
      is_null[i] <- TRUE; pt[i] <- as.numeric(one)
    } else {
      pt[i] <- conv(one, "pt", valueOnly = TRUE)
    }
  }
  rest <- total_pt - sum(pt[!is_null])
  pt[is_null] <- rest * pt[is_null] / sum(pt[is_null])
  pt
}

write_fig3_manifest <- function(plot, manifest_path, width_in, height_in, panel_ids) {
  grob <- patchwork::patchworkGrob(plot)
  probe <- tempfile(fileext = ".pdf")
  grDevices::pdf(probe, width = width_in, height = height_in, useDingbats = FALSE)
  grid::grid.newpage(); grid::grid.draw(grob); grid::grid.force()
  widths_pt  <- resolve_gtable_units(grob$widths,  width_in * 72, "width")
  heights_pt <- resolve_gtable_units(grob$heights, height_in * 72, "height")
  grDevices::dev.off(); unlink(probe)
  stopifnot(all(is.finite(widths_pt)), all(is.finite(heights_pt)),
            abs(sum(widths_pt) - width_in * 72) < 0.01,
            abs(sum(heights_pt) - height_in * 72) < 0.01)

  rows <- grob$layout[grepl("^panel(?:-[0-9]+)?$", grob$layout$name, perl = TRUE), , drop = FALSE]
  rows <- rows[order(rows$t, rows$l), , drop = FALSE]
  stopifnot(nrow(rows) == length(panel_ids))
  x_edges <- c(0, cumsum(widths_pt)); top_edges <- c(0, cumsum(heights_pt))
  total_h <- sum(heights_pt)
  panels <- lapply(seq_len(nrow(rows)), function(i) {
    r <- rows[i, ]
    list(id = panel_ids[i],
         bbox_pt = unname(c(x_edges[r$l], total_h - top_edges[r$b + 1],
                            x_edges[r$r + 1], total_h - top_edges[r$t])),
         grid_id = "patchwork-grid-1",
         row_start = as.integer(r$t - 1), row_stop = as.integer(r$b),
         col_start = as.integer(r$l - 1), col_stop = as.integer(r$r))
  })
  manifest <- list(schema_version = 1L, backend = "r-patchwork",
                   figure = list(width_pt = width_in * 72, height_pt = height_in * 72),
                   panels = panels,
                   row_groups = list(list("A", "B")))
  jsonlite::write_json(manifest, manifest_path, auto_unbox = TRUE, pretty = TRUE, digits = NA)
  invisible(manifest)
}
write_fig3_manifest(fig, "main_figures/Figure3.panel-manifest.json",
                    W_MM / 25.4, H_MM / 25.4, c("A", "B", "C"))
cat("manifest written\n")

## ---------- PNG 像素/DPI 复核 ----------
print(system("sips -g pixelWidth -g pixelHeight -g dpiWidth main_figures/Figure3.png", intern = TRUE))
