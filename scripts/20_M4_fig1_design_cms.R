# Figure 1 — 研究设计（队列×样本类型构成）+ CMS4 腹膜互证
# Panel A: 6 队列 × 6 样本类型瓷砖图（n 标注）
# Panel B: GSE190609 腹膜 vs 肝转移灶 CMS 亚型构成（nearestCMS, n=113）
# Panel C: CMS4×腹膜关联森林图（patient-aware GLMM 主估计 + sample-level Fisher 描述性参照）
set.seed(123)
suppressMessages({
  library(ggplot2); library(patchwork); library(scales)
})

proc <- "data/processed/crc_organotropism"
outdir <- "main_figures"
dir.create(outdir, showWarnings = FALSE)

# ---- 调色板（跨图统一） ----
col_pm   <- "#B4464B"   # 腹膜 哑红
col_lm   <- "#3E6E8E"   # 肝 哑蓝
col_cms4 <- "#C08552"   # CMS4 高亮
g55 <- "grey55"; g75 <- "grey75"; g90 <- "grey90"

FIG_W_MM <- 180; FIG_H_MM <- 105
base_family <- "Helvetica"

theme_fig <- theme_classic(base_size = 6.5, base_family = base_family) +
  theme(
    axis.title  = element_text(size = 7, color = "black"),
    axis.text   = element_text(size = 6.5, color = "black"),
    plot.title  = element_text(size = 7.5, face = "bold", color = "black"),
    plot.caption = element_text(size = 5, color = "grey30", hjust = 0),
    line = element_line(color = "black")
  )

# ================= Panel A: 队列 × 样本类型 =================
# 数据来源：GSE190609/GSE50760/GSE225182 由 processed meta 表重算；
# GSE41258/GSE41568 由 series_matrix 的 tissue/metastatic site 注释重算；
# MSK-MET 用已知数（1,147 例转移灶样本，OS Cox n=1,136）。
m609  <- read.delim(file.path(proc, "GSE190609_meta.tsv"))
m760  <- read.delim(file.path(proc, "GSE50760_meta.tsv"))
m182  <- read.delim(file.path(proc, "GSE225182_meta.tsv"))
stopifnot(nrow(m609) == 113, nrow(m760) == 54, nrow(m182) == 33)

pa <- data.frame(
  cohort = c(rep("GSE190609", 4), rep("GSE225182", 3), rep("GSE50760", 3),
             rep("GSE41258", 4), rep("GSE41568", 4), "MSK-MET"),
  type   = c("PT", "PM", "LM", "Met",
             "PT", "PM", "Normal",
             "PT", "LM", "Normal",
             "PT", "LM", "LungM", "Normal",
             "PM", "LM", "LungM", "Met",
             "Met"),
  n      = c(35, 59, 6, 13,                 # GSE190609: LN 11 + ovary 2 = Met 13
             20, 7, 6,
             18, 18, 18,
             186, 47, 20, 54,               # GSE41258（不含息肉/细胞系/正常肝肺）
             1, 80, 8, 44,                  # GSE41568: omentum 4 + abdom. wall 1 + NA 39
             1147),
  stringsAsFactors = FALSE
)
cohort_lv <- c("GSE190609", "GSE225182", "GSE50760", "GSE41258", "GSE41568", "MSK-MET")
type_lv   <- c("PT", "PM", "LM", "LungM", "Normal", "Met")
pa$cohort <- factor(pa$cohort, levels = rev(cohort_lv))
pa$type   <- factor(pa$type,   levels = type_lv)
type_cols <- c(PT = g55, PM = col_pm, LM = col_lm, LungM = g55, Normal = g75, Met = g55)
pa$txt_col <- ifelse(pa$n >= 146, "white", "black")  # alpha(sqrt) 阈值，深色瓷砖用白字

pA <- ggplot(pa, aes(x = type, y = cohort)) +
  geom_tile(aes(fill = type, alpha = n), color = "white", linewidth = 0.8) +
  geom_text(aes(label = n, color = txt_col), size = 2.1, family = base_family, fontface = "bold") +
  scale_fill_manual(values = type_cols, guide = "none") +
  scale_alpha_continuous(trans = "sqrt", range = c(0.30, 1), guide = "none") +
  scale_color_identity(guide = "none") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL,
       title = "Cohort x sample-type composition") +
  theme_fig +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x.top = element_text(size = 6.5, color = "black", face = "bold"),
    axis.text.y = element_text(size = 6.5, color = "black"),
    panel.background = element_rect(fill = g90, color = NA)
  )

# ================= Panel B: CMS 构成（腹膜 vs 肝） =================
cms <- read.delim(file.path(proc, "GSE190609_CMS_nearest.tsv"), na.strings = c("", "NA"))
stopifnot(nrow(cms) == 113)
cms <- cms[cms$organ %in% c("peritoneal metastasis", "liver metastasis"), ]
stopifnot(all(!grepl(",", cms$cms)))   # 腹膜/肝样本均为单一分型，无歧义标签

pb_tab <- as.data.frame(table(organ = cms$organ, cms = cms$cms))
pb_tab$prop <- pb_tab$Freq / ave(pb_tab$Freq, pb_tab$organ, FUN = sum)
pb_tab$organ <- factor(pb_tab$organ,
  levels = c("peritoneal metastasis", "liver metastasis"),
  labels = c("Peritoneal\n(n = 59)", "Liver\n(n = 6)"))
pb_tab$cms <- factor(pb_tab$cms, levels = c("CMS4", "CMS3", "CMS2", "CMS1"))
cms_cols <- c(CMS1 = g90, CMS2 = g75, CMS3 = g55, CMS4 = col_cms4)
pb_lab <- subset(pb_tab, Freq > 0 & prop >= 0.08)
pb_lab$txt_col <- ifelse(pb_lab$cms %in% c("CMS1", "CMS2"), "black", "white")

pB <- ggplot(pb_tab, aes(x = organ, y = prop, fill = cms)) +
  geom_col(width = 0.62, color = "white", linewidth = 0.4) +
  geom_text(data = pb_lab,
            aes(label = sprintf("%d%%", round(100 * prop)), color = txt_col,
                group = cms),
            position = position_stack(vjust = 0.5),
            size = 2.1, family = base_family, fontface = "bold") +
  scale_fill_manual(values = cms_cols, breaks = c("CMS1", "CMS2", "CMS3", "CMS4"),
                    name = NULL) +
  scale_color_identity(guide = "none") +
  scale_y_continuous(labels = label_percent(accuracy = 1), expand = c(0, 0)) +
  labs(x = NULL, y = "% of metastases",
       title = "CMS subtype by metastatic site") +
  theme_fig +
  theme(legend.position = "right",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6, color = "black"),
        axis.text.x = element_text(size = 6.5, color = "black", lineheight = 0.9))

# ================= Panel C: CMS4 × 腹膜 森林图（GLMM 主估计 + Fisher 描述性） =================
run_fisher <- function(mm) {
  isPM <- mm$organ == "peritoneal metastasis"
  is4  <- mm$cms == "CMS4"
  t <- matrix(c(sum(isPM & is4), sum(isPM & !is4),
                sum(!isPM & is4), sum(!isPM & !is4)), nrow = 2, byrow = TRUE)
  f <- fisher.test(t)
  list(n = nrow(mm), or = unname(f$estimate), lo = f$conf.int[1],
       hi = f$conf.int[2], p = f$p.value)
}
f_near   <- run_fisher(read.delim(file.path(proc, "GSE190609_CMS_nearest.tsv"), na.strings = c("", "NA")))
cms_strict <- read.delim(file.path(proc, "GSE190609_CMS.tsv"), na.strings = c("", "NA"))
f_strict <- run_fisher(cms_strict[!is.na(cms_strict$cms), ])

# 主估计：patient-aware GLMM（随机截距 = 患者），读取 scripts/24 产出，不硬编码
pl <- read.delim("results/revision/cms_patient_level.tsv", quote = "")
glmm <- pl[pl$section == "glmm", ]
g_pt <- glmm[glmm$analysis == "PM_vs_PT", ]
g_lm <- glmm[glmm$analysis == "PM_vs_LM", ]
stopifnot(nrow(g_pt) == 1, nrow(g_lm) == 1)

fmt_p <- function(p) if (p < 0.001) "P < 0.001" else sprintf("P = %.3f", p)

XMAX <- 15   # log 轴上限；PM vs LM 的 CI 上限 (~1268.6) 远超坐标，截断 + 右向箭头
pc <- data.frame(
  label = c("GLMM: PM vs PT",
            "GLMM: PM vs LM (unstable; n LM = 6)",
            sprintf("Fisher: nearestCMS (n = %d)", f_near$n),
            sprintf("Fisher: minPosterior >= 0.5 (n = %d)", f_strict$n)),
  or  = c(g_pt$estimate, g_lm$estimate, f_near$or, f_strict$or),
  lo  = c(g_pt$ci_lo,    g_lm$ci_lo,    f_near$lo, f_strict$lo),
  hi  = c(g_pt$ci_hi,    g_lm$ci_hi,    f_near$hi, f_strict$hi),
  p   = c(g_pt$p_value,  g_lm$p_value,  f_near$p,  f_strict$p),
  grp = factor(c("glmm", "glmm", "fisher", "fisher"), levels = c("glmm", "fisher")),
  y   = c(5.0, 4.0, 2.0, 1),
  stringsAsFactors = FALSE
)
pc$off_scale <- pc$hi > XMAX
pc$ann <- ifelse(pc$off_scale,
                 sprintf("%.1f (%.2f-%.1f); %s", pc$or, pc$lo, pc$hi, sapply(pc$p, fmt_p)),
                 sprintf("%.2f (%.2f-%.2f); %s", pc$or, pc$lo, pc$hi, sapply(pc$p, fmt_p)))
pc_in  <- pc[!pc$off_scale, ]
pc_off <- pc[pc$off_scale, ]
pc_off$x_dot <- 11   # 真实 OR (~32.6) 超出坐标轴：点绘于界内右侧，CI 截断为箭头，实际值见文字

# 分组小标题作为 y 轴标签（ggtext markdown：加粗 + 着色区分 primary / descriptive）
hdr <- data.frame(
  y = c(6.1, 3.1),
  label = c(sprintf("<span style='color:%s'>**Patient-aware GLMM (primary)**</span>", col_pm),
            "<span style='color:grey40'>**Sample-level Fisher (descriptive)**</span>"),
  stringsAsFactors = FALSE
)
axis_br <- c(hdr$y, pc$y)
axis_lb <- c(hdr$label, pc$label)

pC <- ggplot(pc, aes(y = y)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = g55, linewidth = 0.4) +
  geom_errorbarh(data = pc_in, aes(xmin = lo, xmax = hi, color = grp),
                 height = 0.14, linewidth = 0.7) +
  geom_segment(data = pc_off,
               aes(x = lo, xend = XMAX - 0.8, y = y, yend = y, color = grp),
               linewidth = 0.7,
               arrow = grid::arrow(length = unit(1.2, "mm"), type = "closed")) +
  geom_point(data = pc_in,  aes(x = or, color = grp, shape = grp), size = 2.2) +
  geom_point(data = pc_off, aes(x = x_dot, color = grp, shape = grp), size = 2.2) +
  # 注释放置在 CI 线上方、OR=1 参考线右侧，避免与线段碰撞（collision 审计）
  geom_text(aes(x = 1.15, y = y + 0.4, label = ann), hjust = 0, size = 1.9,
            family = base_family, color = "black") +
  scale_color_manual(values = c(glmm = col_pm, fisher = g55), guide = "none") +
  scale_shape_manual(values = c(glmm = 16, fisher = 1), guide = "none") +
  scale_x_log10(breaks = c(0.5, 1, 2, 4, 8, 15),
                labels = c("0.5", "1", "2", "4", "8", "15"),
                limits = c(0.5, XMAX), expand = c(0, 0)) +
  scale_y_continuous(breaks = axis_br, labels = axis_lb, limits = c(0.55, 6.5),
                     expand = c(0, 0)) +
  labs(x = "OR for CMS4, peritoneal vs comparator (log scale)", y = NULL,
       title = "CMS4 enrichment in peritoneal metastases") +
  theme_fig +
  theme(axis.text.y = ggtext::element_markdown(size = 6, color = "black"),
        axis.title.x = element_text(size = 6.5))

# ================= 组装与导出 =================
# 扁平 layout（嵌套 patchwork 会让 alignment QA 找不到 panel 单元格）
lay_design <- c(
  area(t = 1, l = 1, b = 2, r = 1),   # A 纵向跨两行
  area(t = 1, l = 2, b = 1, r = 2),   # B 右上
  area(t = 2, l = 2, b = 2, r = 2)    # C 右下
)
fig <- pA + pB + pC +
  plot_layout(design = lay_design, widths = c(1.12, 0.88)) +
  plot_annotation(
    tag_levels = "A",
    caption = paste0(
      "PT = primary tumor; PM = peritoneal metastasis; LM = liver metastasis; LungM = lung metastasis; ",
      "Met = other/unspecified metastases\n(GSE190609: lymph node 11 + ovary 2; GSE41568: omentum 4, ",
      "abdominal wall 1, unspecified 39; MSK-MET: 1,147 multi-organ metastasis samples, OS Cox n = 1,136)."
    ),
    theme = theme(
      plot.tag = element_text(face = "bold", size = 9, color = "black", family = base_family),
      plot.caption = element_text(size = 5, color = "grey30", hjust = 0, family = base_family)
    )
  )

png_path <- file.path(outdir, "Figure1.png")
pdf_path <- file.path(outdir, "Figure1.pdf")
ragg::agg_png(png_path, width = FIG_W_MM, height = FIG_H_MM, units = "mm",
              res = 600, background = "white")
print(fig); dev.off()
# cairo_pdf 在本机缺 X11 依赖（libXrender）不可用；改用基础 pdf() 设备，
# 同为矢量输出且 Helvetica 以 Type1 嵌入、文字可编辑（已用 PyMuPDF 验证字体）
grDevices::pdf(pdf_path, width = FIG_W_MM / 25.4, height = FIG_H_MM / 25.4,
               family = base_family, bg = "white", useDingbats = FALSE)
print(fig); dev.off()

cat("Exported:", png_path, "and", pdf_path, "\n")
print(system(paste("sips -g pixelWidth -g pixelHeight -g dpiWidth", png_path), intern = TRUE))

# ---- Panel alignment manifest（最终导出尺寸） ----
# skill 自带 write_patchwork_panel_layout 在本机 grid 版本下 null 单位不解析
# （grid.force() 返回 NULL，面板列宽测为 0）。此处按同一 manifest schema，
# 用 deviceLoc() 实测各 panel 视口的设备坐标，供 audit_panel_alignment.py 消费。
source("/Users/wangtao/.agents/skills/nature-figure/scripts/panel_alignment.R")  # 保留接口一致性
write_manifest_measured <- function(plot, manifest_path, width_in, height_in,
                                    panel_ids, column_groups, exemptions) {
  grob <- patchwork::patchworkGrob(plot)
  probe <- tempfile(fileext = ".pdf")
  grDevices::pdf(probe, width = width_in, height = height_in, useDingbats = FALSE)
  on.exit({ grDevices::dev.off(); unlink(probe) }, add = TRUE)
  grid::grid.newpage(); grid::grid.draw(grob); grid::grid.force()
  vp_names <- grid::grid.ls(viewports = TRUE, grobs = FALSE, print = FALSE)$name
  panel_vps <- grep("^panel-[0-9]+[.]", vp_names, value = TRUE)
  lay <- grob$layout[grob$layout$name %in% sub("[.].*$", "", panel_vps), ]
  lay$vp <- panel_vps[match(lay$name, sub("[.].*$", "", panel_vps))]
  lay <- lay[order(lay$t, lay$l), ]
  stopifnot(nrow(lay) == length(panel_ids))
  panels <- lapply(seq_len(nrow(lay)), function(i) {
    grid::downViewport(lay$vp[i])
    bl <- grid::deviceLoc(grid::unit(0, "npc"), grid::unit(0, "npc"))
    tr <- grid::deviceLoc(grid::unit(1, "npc"), grid::unit(1, "npc"))
    grid::upViewport(0)
    list(
      id = panel_ids[i],
      bbox_pt = unname(c(grid::convertUnit(bl$x, "pt", valueOnly = TRUE),
                         grid::convertUnit(bl$y, "pt", valueOnly = TRUE),
                         grid::convertUnit(tr$x, "pt", valueOnly = TRUE),
                         grid::convertUnit(tr$y, "pt", valueOnly = TRUE))),
      grid_id = "patchwork-grid-1",
      row_start = as.integer(lay$t[i] - 1), row_stop = as.integer(lay$b[i]),
      col_start = as.integer(lay$l[i] - 1), col_stop = as.integer(lay$r[i])
    )
  })
  manifest <- list(
    schema_version = 1L, backend = "r-patchwork",
    figure = list(width_pt = width_in * 72, height_pt = height_in * 72),
    panels = panels,
    column_groups = column_groups,
    exemptions = exemptions
  )
  jsonlite::write_json(manifest, manifest_path, auto_unbox = TRUE, pretty = TRUE, digits = NA)
  invisible(manifest)
}
write_manifest_measured(
  fig,
  manifest_path = file.path(outdir, "Figure1.panel-layout.json"),
  width_in = FIG_W_MM / 25.4, height_in = FIG_H_MM / 25.4,
  panel_ids = c("a", "b", "c"),
  column_groups = list(c("b", "c")),
  exemptions = list(list(panels = I("a"), checks = I("all"),
                         reason = "Hero panel A intentionally spans the full left column (both grid rows); asymmetric design, excluded from equal-size comparisons."))
)
cat("Manifest written.\n")

# ---- 重算数值回报 ----
cat(sprintf("PANEL_C GLMM PM vs PT : OR=%.2f [%.2f, %.2f] %s (n patients=%d, n samples=%d)\n",
            g_pt$estimate, g_pt$ci_lo, g_pt$ci_hi, fmt_p(g_pt$p_value),
            g_pt$n_patients, g_pt$n_samples))
cat(sprintf("PANEL_C GLMM PM vs LM : OR=%.2f [%.2f, %.2f] %s (CI off-scale, truncated with arrow)\n",
            g_lm$estimate, g_lm$ci_lo, g_lm$ci_hi, fmt_p(g_lm$p_value)))
cat(sprintf("PANEL_C nearestCMS : n=%d OR=%.2f [%.2f, %.2f] %s\n",
            f_near$n, f_near$or, f_near$lo, f_near$hi, fmt_p(f_near$p)))
cat(sprintf("PANEL_C strict>=0.5: n=%d OR=%.2f [%.2f, %.2f] %s\n",
            f_strict$n, f_strict$or, f_strict$lo, f_strict$hi, fmt_p(f_strict$p)))
cat(sprintf("PANEL_B PM CMS4 = %.1f%% ; LM CMS4 = %.1f%%\n",
            100 * mean(read.delim(file.path(proc, "GSE190609_CMS_nearest.tsv"))$cms[
              read.delim(file.path(proc, "GSE190609_CMS_nearest.tsv"))$organ == "peritoneal metastasis"] == "CMS4"),
            100 * mean(read.delim(file.path(proc, "GSE190609_CMS_nearest.tsv"))$cms[
              read.delim(file.path(proc, "GSE190609_CMS_nearest.tsv"))$organ == "liver metastasis"] == "CMS4")))
