# Figure 4 - Honest negative result and hepatocyte-contamination warning
# Panel A: five AUC configurations with bootstrap 95% CI, Wilcoxon P and permutation P
# Panel B: gene retention after epithelial-purity adjustment
# Panel C: hepatocyte-marker log2FC (LM vs PT) across three cohorts
set.seed(123)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
})

base_dir <- "data/processed/crc_organotropism"
out_dir  <- "main_figures"
dir.create(out_dir, showWarnings = FALSE)

col_pm  <- "#B4464B"  # peritoneum, muted red
col_lm  <- "#3E6E8E"  # liver, muted blue
col_mid <- "grey55"

font_family <- "Helvetica"
base_size   <- 6.5

theme_fig <- theme_classic(base_size = base_size, base_family = font_family) +
  theme(
    axis.text        = element_text(colour = "black", size = base_size),
    axis.title       = element_text(colour = "black", size = base_size + 0.5),
    plot.title       = element_text(colour = "black", size = base_size + 1, face = "bold"),
    plot.subtitle    = element_text(colour = "grey25", size = base_size - 0.5),
    legend.text      = element_text(colour = "black", size = base_size - 0.5),
    legend.title     = element_text(colour = "black", size = base_size - 0.5),
    legend.key.size  = unit(3.2, "mm"),
    axis.line        = element_line(linewidth = 0.3),
    axis.ticks       = element_line(linewidth = 0.3)
  )


## ---------------------------------------------------------------- Panel A
# Five configurations with bootstrap 95% CI error bars (2,000 resamples) and two
# clearly separated P annotations per bar: Wilcoxon rank-sum P (group separation)
# and permutation empirical P (1,000 random gene sets, specificity null).
auc_main <- read.delim(file.path(base_dir, "ssgsea/layer3_auc_summary.tsv"),
                       stringsAsFactors = FALSE)
auc_loco <- read.delim(file.path(base_dir, "ssgsea/layer3_LOCO_auc.tsv"),
                       stringsAsFactors = FALSE)
auc_perm <- read.delim(file.path(base_dir, "ssgsea/revision_permutation_auc.tsv"),
                       stringsAsFactors = FALSE)

pm_app <- auc_main[grepl("主检验1", auc_main$tag), ]
lm_app <- auc_main[grepl("主检验2", auc_main$tag), ]
lm_loco_main <- auc_loco[grepl("^T1:", auc_loco$tag), ]       # LOCO main row
lm_sen <- auc_loco[grepl("T1-敏感性", auc_loco$tag), ]        # restricted subset
pm_sen <- auc_loco[grepl("^T2:", auc_loco$tag), ]             # cross-cohort

perm_main <- auc_perm[auc_perm$contrast == "main", ]
get_perm <- function(cfg) perm_main[perm_main$config_id == cfg, ]
perm_lm_app  <- get_perm("A_apparent_LMcore")
perm_lm_loco <- get_perm("C_LOCO_main_LMcore")
perm_lm_sen  <- get_perm("D_LOCO_sens_LMcore")
perm_pm_app  <- get_perm("B_apparent_PMcore")
perm_pm_sen  <- get_perm("E_LOCO_T2_PMcore")

panel_a_df <- data.frame(
  group  = rep(c("LM_core", "PM_core"), times = c(3, 2)),
  config = c("Apparent", "LOCO\nmain", "LOCO\nrestricted\nsubset",
             "Apparent", "Cross-\ncohort"),
  mode   = c("Apparent", "LOCO / sensitivity", "LOCO / sensitivity",
             "Apparent", "LOCO / sensitivity"),
  auc      = c(lm_app$auc, lm_loco_main$auc, lm_sen$auc, pm_app$auc, pm_sen$auc),
  p_wilcox = c(lm_app$p, lm_loco_main$p, lm_sen$p, pm_app$p, pm_sen$p),
  ci_lo  = c(perm_lm_app$boot_ci_lo, perm_lm_loco$boot_ci_lo, perm_lm_sen$boot_ci_lo,
             perm_pm_app$boot_ci_lo, perm_pm_sen$boot_ci_lo),
  ci_hi  = c(perm_lm_app$boot_ci_hi, perm_lm_loco$boot_ci_hi, perm_lm_sen$boot_ci_hi,
             perm_pm_app$boot_ci_hi, perm_pm_sen$boot_ci_hi),
  p_perm = c(perm_lm_app$emp_p_greater, perm_lm_loco$emp_p_greater,
             perm_lm_sen$emp_p_greater, perm_pm_app$emp_p_greater,
             perm_pm_sen$emp_p_greater),
  fill   = rep(c(col_lm, col_pm), times = c(3, 2)),
  stringsAsFactors = FALSE
)
# guard: permutation-table AUC must match the source TSV AUC for every config
stopifnot(isTRUE(all.equal(
  panel_a_df$auc,
  c(perm_lm_app$auc_real, perm_lm_loco$auc_real, perm_lm_sen$auc_real,
    perm_pm_app$auc_real, perm_pm_sen$auc_real),
  tolerance = 1e-9
)))
panel_a_df$group <- factor(panel_a_df$group, levels = c("LM_core", "PM_core"))
panel_a_df$config <- factor(panel_a_df$config,
                            levels = c("Apparent", "LOCO\nmain",
                                       "LOCO\nrestricted\nsubset", "Cross-\ncohort"))
panel_a_df$mode <- factor(panel_a_df$mode,
                          levels = c("Apparent", "LOCO / sensitivity"))
fmt_p_plain <- function(p) if (p < 0.001) "P < 0.001" else sprintf("P = %.3f", p)
panel_a_df$plab <- sprintf("Wilcoxon %s\nperm %s",
                           vapply(panel_a_df$p_wilcox, fmt_p_plain, character(1)),
                           vapply(panel_a_df$p_perm, fmt_p_plain, character(1)))

strip_lab <- c("LM_core" = "LM_core\nliver-met PT vs peritoneal-met PT",
               "PM_core" = "PM_core\nperitoneal-met PT vs liver-met PT")

pA <- ggplot(panel_a_df, aes(x = config, y = auc)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = col_mid, linewidth = 0.35) +
  geom_col(aes(fill = fill, alpha = mode), width = 0.62,
           colour = "grey20", linewidth = 0.25) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.22, linewidth = 0.3, colour = "grey15") +
  geom_text(aes(label = sprintf("%.2f", auc), y = ci_hi + 0.04), vjust = 0,
            size = 1.9, family = font_family, colour = "black") +
  geom_text(aes(label = plab, y = ci_hi + 0.115), vjust = 0,
            size = 1.6, lineheight = 0.85, family = font_family, colour = "black") +
  geom_text(data = data.frame(group = factor("LM_core", levels = c("LM_core", "PM_core"))),
            aes(x = 0.55, y = 1.29), inherit.aes = FALSE,
            label = "dashed line = chance (AUC = 0.5)",
            hjust = 0, size = 1.6, colour = "grey30", family = font_family) +
  facet_grid(cols = vars(group), scales = "free_x", space = "free_x",
             labeller = as_labeller(strip_lab)) +
  scale_fill_identity() +
  scale_alpha_manual(values = c("Apparent" = 1, "LOCO / sensitivity" = 0.38)) +
  scale_y_continuous(limits = c(0, 1.36), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(x = NULL, y = "AUC (primary tumour -> metastatic site)",
       title = "Primary-tumour transcriptome does not predict metastatic site",
       subtitle = "Error bars: bootstrap 95% CI (2,000 resamples); perm P: 1,000 random gene sets",
       alpha = NULL) +
  theme_fig +
  theme(legend.position = c(0.985, 0.985),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", colour = NA),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black", size = base_size - 0.5, lineheight = 0.9),
        axis.text.x = element_text(size = base_size - 1, lineheight = 0.85),
        panel.spacing.x = unit(4, "mm"))

## ---------------------------------------------------------------- Panel B
pm_total <- nrow(read.delim(file.path(base_dir, "meta/peritoneal_tropism_genes.tsv")))
lm_total <- nrow(read.delim(file.path(base_dir, "meta/liver_tropism_genes.tsv")))
pm_ret   <- nrow(read.delim(file.path(base_dir, "meta/PM_core_purity_retained.tsv")))
lm_ret   <- nrow(read.delim(file.path(base_dir, "meta/LM_core_purity_retained.tsv")))

panel_b_df <- data.frame(
  set  = factor(c("PM_core", "LM_core"), levels = c("PM_core", "LM_core")),
  retained = c(pm_ret, lm_ret),
  total    = c(pm_total, lm_total),
  fill     = c(col_pm, col_lm)
)
panel_b_df$frac <- panel_b_df$retained / panel_b_df$total

pB <- ggplot(panel_b_df, aes(x = set, y = frac * 100, fill = fill)) +
  geom_col(width = 0.58, colour = "grey20", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%d/%d", retained, total), y = frac * 100 - 7),
            size = 2.1, family = font_family, colour = "black") +
  geom_text(aes(label = sprintf("%.0f%%", frac * 100), y = frac * 100 + 5),
            size = 2.1, family = font_family, colour = "black", fontface = "bold") +
  scale_fill_identity() +
  scale_y_continuous(limits = c(0, 112), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(x = NULL, y = "Genes retained after\npurity adjustment (%)",
       title = "Liver-tropism signal survives\npurity correction; peritoneal partly") +
  theme_fig

## ---------------------------------------------------------------- Panel C
cohorts <- c("GSE50760", "GSE190609", "GSE41258")
markers <- c("ALB", "APOB", "APOA1", "HP", "SERPINA1")

de_list <- lapply(cohorts, function(cc) {
  d <- read.delim(gzfile(file.path(base_dir, "de", paste0(cc, "_LM_vs_PT.tsv.gz"))),
                  stringsAsFactors = FALSE)
  d$cohort <- cc
  d[d$symbol %in% markers, c("symbol", "logFC", "adj.P.Val", "cohort")]
})
panel_c_df <- do.call(rbind, de_list)
panel_c_df$symbol <- factor(panel_c_df$symbol,
                            levels = markers[order(vapply(markers, function(g) {
                              max(panel_c_df$logFC[panel_c_df$symbol == g])
                            }, numeric(1)))])
panel_c_df$cohort <- factor(panel_c_df$cohort, levels = cohorts)

coh_shapes <- c("GSE50760" = 16, "GSE190609" = 17, "GSE41258" = 15)

pC <- ggplot(panel_c_df, aes(x = logFC, y = symbol)) +
  geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.3) +
  geom_vline(xintercept = 2, linetype = "dashed", colour = col_mid, linewidth = 0.35) +
  geom_point(aes(shape = cohort), colour = col_lm, size = 2.2, stroke = 0.6) +
  annotate("text", x = 2.2, y = 0.52, label = "typical DE magnitude (log2FC = 2)",
           hjust = 0, size = 1.7, colour = "grey30", family = font_family) +
  scale_shape_manual(values = coh_shapes) +
  scale_x_continuous(limits = c(-0.5, 12.5), breaks = seq(0, 12, 2), expand = c(0.01, 0)) +
  labs(x = "log2 fold change (liver metastasis vs primary tumour)",
       y = NULL, shape = "Cohort",
       title = "Hepatocyte contamination: liver markers enriched in liver metastases vs primary tumours") +
  theme_fig +
  theme(legend.position = c(0.88, 0.30),
        legend.background = element_rect(fill = "white", colour = "grey85", linewidth = 0.2),
        legend.spacing.y = unit(0.5, "mm"))

## ---------------------------------------------------------------- assemble
# flat design keeps panel cells as "panel-1..3" so the alignment helper can see them
fig <- wrap_plots(A = pA, B = pB, C = pC, design = "AB\nCC") +
  plot_layout(widths = c(1.55, 1), heights = c(1, 0.88)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", colour = "black",
                                size = base_size + 2, family = font_family))

fig_w_in <- 180 / 25.4
fig_h_in <- 150 / 25.4

png_path <- file.path(out_dir, "Figure4.png")
pdf_path <- file.path(out_dir, "Figure4.pdf")

ragg::agg_png(png_path, width = fig_w_in, height = fig_h_in,
              units = "in", res = 600, background = "white")
print(fig)
dev.off()

# cairo_pdf unavailable on this machine (R cairo module needs missing X11 libs);
# quartz(type="pdf") gives the same contract: vector output with editable text.
quartz(type = "pdf", file = pdf_path, width = fig_w_in, height = fig_h_in,
       family = font_family)
print(fig)
dev.off()

## ---------------------------------------------------------------- QA: alignment manifest
# The skill helper (panel_alignment.R) collapses panel bboxes to zero under this
# R 4.2/grid build (convertWidth of null units returns 0). Same JSON schema is
# therefore written here from absolute device coordinates via grid::deviceLoc.
library(grid)
write_panel_manifest <- function(plot, manifest_path, width_in, height_in,
                                 panel_ids, row_groups = NULL) {
  probe <- tempfile(fileext = ".pdf")
  pdf(probe, width = width_in, height = height_in, useDingbats = FALSE)
  grid.newpage()
  grid.draw(patchwork::patchworkGrob(plot))
  grid.force()
  vps <- grid.ls(viewports = TRUE, grobs = FALSE, print = FALSE)
  # facetted patchwork cells yield one viewport per facet ("panel-1-1.",
  # "panel-1-2.", ...); group them by patchwork cell id (first number) and
  # merge each cell into a single union bbox
  vp_names <- vps$name[grepl("^panel-[0-9]+(-[0-9]+)?[.]", vps$name)]
  cell_id  <- as.integer(sub("^panel-([0-9]+).*", "\\1", vp_names))
  cells    <- sort(unique(cell_id))
  if (length(cells) != length(panel_ids)) {
    dev.off(); unlink(probe)
    stop("panel viewport count does not match panel_ids")
  }
  panels <- lapply(seq_along(cells), function(i) {
    vns <- vp_names[cell_id == cells[i]]
    x0 <- y0 <- Inf; x1 <- y1 <- -Inf
    for (vn in vns) {
      downViewport(vn)
      l00 <- deviceLoc(unit(0, "npc"), unit(0, "npc"))
      l11 <- deviceLoc(unit(1, "npc"), unit(1, "npc"))
      upViewport(0)
      x0 <- min(x0, convertX(l00$x, "pt", valueOnly = TRUE))
      y0 <- min(y0, convertY(l00$y, "pt", valueOnly = TRUE))
      x1 <- max(x1, convertX(l11$x, "pt", valueOnly = TRUE))
      y1 <- max(y1, convertY(l11$y, "pt", valueOnly = TRUE))
    }
    list(id = panel_ids[i], bbox_pt = c(x0, y0, x1, y1),
         grid_id = "patchwork-grid-1")
  })
  dev.off(); unlink(probe)
  manifest <- list(schema_version = 1L, backend = "r-patchwork",
                   figure = list(width_pt = width_in * 72, height_pt = height_in * 72),
                   panels = panels, exemptions = list())
  if (!is.null(row_groups)) manifest$row_groups <- lapply(row_groups, as.character)
  jsonlite::write_json(manifest, manifest_path, auto_unbox = TRUE, pretty = TRUE, digits = NA)
  invisible(manifest)
}

manifest_path <- file.path(out_dir, "Figure4.panel-layout.json")
write_panel_manifest(fig, manifest_path, fig_w_in, fig_h_in,
                     panel_ids = c("A", "B", "C"), row_groups = list(c("A", "B")))

system(paste("sips -g pixelWidth -g pixelHeight -g dpiWidth", shQuote(png_path)))

cat("Panel A values:\n"); print(panel_a_df[, c("group", "config", "mode", "auc",
                                              "ci_lo", "ci_hi", "p_wilcox", "p_perm")])
cat("Panel B values:\n"); print(panel_b_df[, c("set", "retained", "total")])
cat("Panel C values:\n"); print(panel_c_df[order(panel_c_df$symbol, panel_c_df$cohort), ])
