# Figure 2 — 转移灶器官生态位转录印记（跨队列随机效应 meta 分析）
# Panel A: meta 效应量火山图（PM vs PT）
# Panel B: 腹膜趋向 top 基因森林图（GSE190609 + GSE225182 + meta 菱形）
# Panel C: 肝趋向 top 基因森林图（4 个 LM 队列 + meta 菱形）
set.seed(123)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

dir.create("main_figures", showWarnings = FALSE)

COL_PERI  <- "#B4464B"
COL_LIVER <- "#3E6E8E"
GREY_PT   <- "grey55"
GREY_MID  <- "grey75"
FONT <- "Helvetica"

meta_dir <- "data/processed/crc_organotropism/meta"
de_dir   <- "data/processed/crc_organotropism/de"

meta_pm     <- fread(file.path(meta_dir, "meta_PM_vs_PT_RNAseq.tsv.gz"))
meta_lm     <- fread(file.path(meta_dir, "meta_LM_vs_PT_mixed.tsv.gz"))
peri_genes  <- fread(file.path(meta_dir, "peritoneal_tropism_genes.tsv"))
liver_genes <- fread(file.path(meta_dir, "liver_tropism_genes.tsv"))

read_de <- function(f) {
  d <- fread(file.path(de_dir, f))
  d[, se := abs(logFC / t)]
  d[, .(symbol, logFC, se)]
}

# ---- sanity counts（核对预期数值） ---------------------------------------------
pm_sig <- meta_pm[fdr < 0.05 & abs(beta) >= 1]
lm_sig <- meta_lm[fdr < 0.05 & abs(beta) >= 1]
cat(sprintf("PM program: %d genes FDR<0.05 & |beta|>=1 (%d sign-concordant)\n",
            nrow(pm_sig), sum(pm_sig$sign_concord == 1)))
cat(sprintf("LM program: %d genes FDR<0.05 & |beta|>=1 (%d sign-concordant)\n",
            nrow(lm_sig), sum(lm_sig$sign_concord == 1)))
cat(sprintf("Tropism lists: peritoneal n=%d, liver n=%d\n",
            nrow(peri_genes), nrow(liver_genes)))

# ---- Panel A: volcano ----------------------------------------------------------
va <- copy(meta_pm)
va[, group := "Other"]
va[gene %in% liver_genes$gene, group := "Liver-tropism (n=57)"]
va[gene %in% peri_genes$gene,  group := "Peritoneal-tropism (n=49)"]
va[, group := factor(group, levels = c("Other", "Liver-tropism (n=57)", "Peritoneal-tropism (n=49)"))]
va[, neglogfdr := -log10(pmax(fdr, 1e-300))]

peri_lab <- c("FABP4", "CD36", "PLAC9", "RBP7", "ADAMTS15", "PI16", "SRPX", "LPL")
liver_lab <- va[group == "Liver-tropism (n=57)" & fdr < 0.05][order(fdr)]$gene
cat("Liver-tropism genes FDR<0.05 in PM meta (labelled in Panel A):",
    paste(liver_lab, collapse = ", "), "\n")
# ggrepel 避让全部前景色点（腹膜49+肝57）：未标注的点给空标签，仅作避让锚点
repel_dt <- va[group != "Other"]
repel_dt[, lab := ifelse(gene %in% c(peri_lab, liver_lab), gene, "")]
lab_dt <- repel_dt[lab != ""]

pA <- ggplot(va, aes(beta, neglogfdr)) +
  geom_point(data = va[group == "Other"], colour = GREY_MID, size = 0.4, alpha = 0.5, stroke = 0) +
  geom_point(data = va[group == "Liver-tropism (n=57)"], colour = COL_LIVER, size = 1.1, alpha = 0.9, stroke = 0) +
  geom_point(data = va[group == "Peritoneal-tropism (n=49)"], colour = COL_PERI, size = 1.1, alpha = 0.95, stroke = 0) +
  geom_hline(yintercept = -log10(0.05), colour = GREY_PT, linewidth = 0.25, linetype = "dashed") +
  geom_text_repel(
    data = repel_dt, aes(label = lab, colour = group),
    fontface = "italic", size = 6 / .pt, seed = 123, family = FONT,
    box.padding = 0.4, point.padding = 0.45, min.segment.length = 0,
    force_pull = 1.1, max.time = 2, max.iter = 20000,
    segment.size = 0.2, segment.colour = "grey60", max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = c("Peritoneal-tropism (n=49)" = COL_PERI,
                                 "Liver-tropism (n=57)" = COL_LIVER)) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.5,
           label = "282 genes FDR < 0.05, |β| ≥ 1\n(280 sign-concordant)",
           size = 6 / .pt, colour = "black", lineheight = 0.95, family = FONT) +
  labs(x = "Meta β (peritoneal metastasis vs primary tumour)",
       y = expression(-log[10] * "(FDR)")) +
  theme_bw(base_size = 7, base_family = FONT) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 0.4),
    axis.text = element_text(colour = "black", size = 6),
    axis.title = element_text(colour = "black", size = 7),
    plot.margin = margin(4, 4, 2, 2)
  )

# ---- forest builder（手动 y 布局，单 panel，兼容 alignment 审计） ----------------
# meta SE 由 pval 反推：z = qnorm(1 - p/2), se = |beta| / z
meta_se <- function(beta, pval) abs(beta) / qnorm(1 - pval / 2)

build_forest_dt <- function(genes, meta_tab, cohort_list) {
  rows <- list()
  for (g in genes) {
    for (nm in names(cohort_list)) {
      d <- cohort_list[[nm]][symbol == g]
      if (nrow(d) == 0) {
        cat(sprintf("  [warn] %s missing in %s\n", g, nm))
        next
      }
      rows[[length(rows) + 1]] <- data.table(
        gene = g, study = nm, est = d$logFC[1],
        lo = d$logFC[1] - 1.96 * d$se[1], hi = d$logFC[1] + 1.96 * d$se[1],
        kind = "cohort"
      )
    }
    m <- meta_tab[gene == g]
    se_m <- meta_se(m$beta, m$pval)
    fdr_lab <- ifelse(m$fdr < 0.001, "FDR < 0.001", sprintf("FDR = %.3f", m$fdr))
    rows[[length(rows) + 1]] <- data.table(
      gene = g, study = "Meta (RE)", est = m$beta,
      lo = m$beta - 1.96 * se_m, hi = m$beta + 1.96 * se_m,
      kind = "meta", fdr_text = sprintf("β = %.2f, %s", m$beta, fdr_lab)
    )
  }
  rbindlist(rows, fill = TRUE)
}

forest_plot <- function(dt, gene_order, study_levels, diamond_col, title) {
  dt <- copy(dt)
  dt[, gene := factor(gene, levels = gene_order)]
  dt[, study := factor(study, levels = study_levels)]
  setorder(dt, gene, study)
  m_per_gene <- length(study_levels)
  # 每个基因块 = 1 个基因名标签行（无几何元素）+ m 个研究行（末行为 meta）
  dt[, pos := (as.integer(gene) - 1) * (m_per_gene + 1) + 1 + as.integer(study)]
  dt[, y := -pos]

  # meta 菱形多边形（宽度 = 95% CI）
  dtm <- dt[kind == "meta"][, id := .I]
  dia <- dtm[, .(x = c(lo, est, hi, est), y = c(y, y + 0.32, y, y - 0.32)),
             by = .(gene, id)][, .(x, y), by = .(gene, id)]

  glab <- dt[, .(y_lab = max(y) + 1), by = gene]   # 标签行 = 块内第一行

  xmin <- min(dt$lo, na.rm = TRUE)
  xmax <- max(dt$hi, na.rm = TRUE)
  span <- xmax - xmin
  x_lab <- max(xmin, 0) + 0.025 * span   # 基因名：标签行左侧，避开 x=0 虚线
  x_txt <- xmax + 0.33 * span            # meta 统计文本右对齐位置
  x_hi  <- xmax + 0.35 * span

  ggplot(dt, aes(y = y)) +
    geom_vline(xintercept = 0, colour = GREY_PT, linewidth = 0.25, linetype = "dashed") +
    geom_errorbarh(data = dt[kind == "cohort"], aes(xmin = lo, xmax = hi),
                   height = 0, colour = GREY_PT, linewidth = 0.35) +
    geom_point(data = dt[kind == "cohort"], aes(x = est),
               colour = GREY_PT, size = 1.0, stroke = 0) +
    geom_polygon(data = dia, aes(x = x, y = y, group = interaction(gene, id)),
                 fill = diamond_col, colour = diamond_col, linewidth = 0.2) +
    geom_text(data = dt[kind == "meta"], aes(x = x_txt, label = fdr_text),
              hjust = 1, size = 5 / .pt, colour = "black", family = FONT) +
    geom_text(data = glab, aes(x = x_lab, y = y_lab, label = gene),
              hjust = 0, fontface = "italic", size = 5.5 / .pt, colour = "black", family = FONT) +
    scale_y_continuous(breaks = dt$y, labels = dt$study,
                       expand = expansion(mult = c(0.008, 0.04))) +
    coord_cartesian(xlim = c(xmin - 0.02 * span, x_hi), clip = "off") +
    labs(x = "log2 fold change (95% CI)", y = NULL, title = title) +
    theme_bw(base_size = 7, base_family = FONT) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", linewidth = 0.4),
      axis.text.y = element_text(colour = "black", size = 5),
      axis.text.x = element_text(colour = "black", size = 6),
      axis.title = element_text(colour = "black", size = 6.5),
      plot.title = element_text(size = 7, face = "bold", colour = "black", margin = margin(b = 2)),
      plot.margin = margin(2, 2, 2, 2)
    )
}

# ---- Panel B: peritoneal forest -------------------------------------------------
# 7 个基因（6–8 区间内）：按 meta FDR 舍弃最弱的 SRPX（FDR=2.96e-4），保留 LPL（脂肪谱系标志）
peri_top <- c("FABP4", "CD36", "PLAC9", "RBP7", "ADAMTS15", "PI16", "LPL")
peri_top <- peri_top[peri_top %in% peri_genes$gene]
cat("Panel B genes:", paste(peri_top, collapse = ", "), "\n")
dtB <- build_forest_dt(
  peri_top, meta_pm,
  list(GSE190609 = read_de("GSE190609_PM_vs_PT.tsv.gz"),
       GSE225182 = read_de("GSE225182_PM_vs_PT.tsv.gz"))
)
pB <- forest_plot(dtB, peri_top, c("GSE190609", "GSE225182", "Meta (RE)"),
                  COL_PERI, "Peritoneal-tropism genes (PM vs PT)")

# ---- Panel C: liver forest ------------------------------------------------------
liver_top <- liver_genes[order(-beta)]$gene
liver_top <- liver_top[liver_top %in% c("FGG", "HP", "APOC3", "GC", "APOA1", "F9", "SERPINA1")]
cat("Panel C genes:", paste(liver_top, collapse = ", "), "\n")
dtC <- build_forest_dt(
  liver_top, meta_lm,
  list(GSE190609 = read_de("GSE190609_LM_vs_PT.tsv.gz"),
       GSE50760  = read_de("GSE50760_LM_vs_PT.tsv.gz"),
       GSE41258  = read_de("GSE41258_LM_vs_PT.tsv.gz"),
       GSE41568  = read_de("GSE41568_LM_vs_PT.tsv.gz"))
)
pC <- forest_plot(dtC, liver_top,
                  c("GSE190609", "GSE50760", "GSE41258", "GSE41568", "Meta (RE)"),
                  COL_LIVER, "Liver-tropism genes (LM vs PT)")

# ---- assemble（扁平 wrap_plots 布局，panel viewport 命名为 panel-1..3） ------------
fig <- wrap_plots(A = pA, B = pB, C = pC, design = "AB\nAC") +
  plot_layout(widths = c(1, 1.12), heights = c(1, 1.5)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 8, colour = "black", family = FONT))

W_IN <- 180 / 25.4
H_IN <- 170 / 25.4

ragg::agg_png("main_figures/Figure2.png", width = W_IN, height = H_IN,
              units = "in", res = 600, background = "white")
print(fig)
dev.off()

# cairo_pdf 在本机 R 4.2 因缺 XQuartz 动态库不可用 → quartz(type="pdf")（同为矢量、文字可选可编辑）
grDevices::quartz(type = "pdf", file = "main_figures/Figure2.pdf",
                  width = W_IN, height = H_IN, family = FONT)
print(fig)
dev.off()

# ---- alignment manifest ---------------------------------------------------------
# skill 助手（panel_alignment.R）在本机 R 4.2/grid 下 convertWidth(null) 归零导致 bbox 塌缩；
# 改用与 scripts/23_M4_fig4 相同的 deviceLoc 方案从绝对设备坐标写同 schema 的 JSON。
source("/Users/wangtao/.agents/skills/nature-figure/scripts/panel_alignment.R")
library(grid)
write_panel_manifest <- function(plot, manifest_path, width_in, height_in,
                                 panel_ids, column_groups = NULL) {
  probe <- tempfile(fileext = ".pdf")
  pdf(probe, width = width_in, height = height_in, useDingbats = FALSE)
  grid.newpage()
  grid.draw(patchwork::patchworkGrob(plot))
  grid.force()
  vps <- grid.ls(viewports = TRUE, grobs = FALSE, print = FALSE)
  vp_names <- vps$name[grepl("^panel-[0-9]+[.]", vps$name)]
  if (length(vp_names) != length(panel_ids)) {
    dev.off(); unlink(probe)
    stop("panel viewport count does not match panel_ids")
  }
  panels <- lapply(seq_along(vp_names), function(i) {
    downViewport(vp_names[i])
    l00 <- deviceLoc(unit(0, "npc"), unit(0, "npc"))
    l11 <- deviceLoc(unit(1, "npc"), unit(1, "npc"))
    upViewport(0)
    list(id = panel_ids[i],
         bbox_pt = c(convertX(l00$x, "pt", valueOnly = TRUE),
                     convertY(l00$y, "pt", valueOnly = TRUE),
                     convertX(l11$x, "pt", valueOnly = TRUE),
                     convertY(l11$y, "pt", valueOnly = TRUE)),
         grid_id = "patchwork-grid-1")
  })
  dev.off(); unlink(probe)
  manifest <- list(schema_version = 1L, backend = "r-patchwork",
                   figure = list(width_pt = width_in * 72, height_pt = height_in * 72),
                   panels = panels, exemptions = list())
  if (!is.null(column_groups)) manifest$column_groups <- lapply(column_groups, as.character)
  jsonlite::write_json(manifest, manifest_path, auto_unbox = TRUE, pretty = TRUE, digits = NA)
  invisible(manifest)
}

write_panel_manifest(fig, "main_figures/Figure2.panel-layout.json",
                     W_IN, H_IN, panel_ids = c("A", "B", "C"),
                     column_groups = list(c("B", "C")))

cat("Done. Files written to main_figures/\n")
