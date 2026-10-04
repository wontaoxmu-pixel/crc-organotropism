# Revision: 修正患者级突变表（腹膜受累 vs 无），补齐 APC 并统一 9 基因 BH FDR 口径
# 背景: data/processed/crc_organotropism/mskmet/patient_mutation_peritoneal.tsv 为 scripts/15 旧产出，
#       仅 8 基因（漏 APC）且 FDR 为 8 基因 BH。本脚本从原始 JSON 重算 9 基因版并覆盖写回该文件。
#       重算逻辑与 scripts/22_M4_fig3_genomic_axis.R Panel B 完全一致（不改 scripts/15，避免牵动 M3 链路）。
# 输出: 覆盖 data/processed/crc_organotropism/mskmet/patient_mutation_peritoneal.tsv（9 基因全字段）
# 断言: APC OR∈[0.45,0.47]、P∈[3.0e-4,3.4e-4]；SMAD4 FDR∈[0.20,0.23]；RNF43 P∈[0.13,0.15]
set.seed(123)
suppressMessages({
  library(data.table); library(jsonlite)
})

raw  <- "data/raw/crc_organotropism"
proc <- "data/processed/crc_organotropism/mskmet"

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

genes <- c("RNF43","SMAD4","PTEN","KRAS","BRAF","NRAS","TP53","PIK3CA","APC")
res <- rbindlist(lapply(genes, function(g) {
  ids <- unique(mut_p[gene == g]$patientId)
  a <- sum(pw$hasPM & pw$patientId %in% ids);  b <- sum(pw$hasPM & !pw$patientId %in% ids)
  c1 <- sum(!pw$hasPM & pw$patientId %in% ids); d <- sum(!pw$hasPM & !pw$patientId %in% ids)
  ft <- fisher.test(matrix(c(a, b, c1, d), nrow = 2, byrow = TRUE))
  data.table(gene = g, PM_mut = a, PM_n = a + b, nonPM_mut = c1, nonPM_n = c1 + d,
             pct_PM = 100 * a / (a + b), pct_nonPM = 100 * c1 / (c1 + d),
             OR = as.numeric(ft$estimate), lo = ft$conf.int[1], hi = ft$conf.int[2],
             p = ft$p.value)
}))
res[, fdr := p.adjust(p, "BH")]   # 9 基因集合统一 BH
res[, `:=`(pct_PM = round(pct_PM, 1), pct_nonPM = round(pct_nonPM, 1))]

## ---------- 断言（对齐主稿 / Figure 3B 口径） ----------
stopifnot(nrow(res) == 9, setequal(res$gene, genes))
apc   <- res[gene == "APC"]
smad4 <- res[gene == "SMAD4"]
rnf43 <- res[gene == "RNF43"]
stopifnot(
  apc$OR   >= 0.45    && apc$OR   <= 0.47,
  apc$p    >= 3.0e-4  && apc$p    <= 3.4e-4,
  smad4$fdr >= 0.20   && smad4$fdr <= 0.23,
  rnf43$p  >= 0.13    && rnf43$p  <= 0.15
)
cat("断言全部通过：\n")
cat(sprintf("  APC   OR = %.3f [%.3f–%.3f]  P = %.2e  FDR = %.4f\n",
            apc$OR, apc$lo, apc$hi, apc$p, apc$fdr))
cat(sprintf("  SMAD4 OR = %.3f  P = %.4f  FDR = %.3f\n", smad4$OR, smad4$p, smad4$fdr))
cat(sprintf("  RNF43 OR = %.3f  P = %.3f  FDR = %.3f\n", rnf43$OR, rnf43$p, rnf43$fdr))

## ---------- 覆盖写回（修正 8 基因旧文件，预期行为） ----------
out <- file.path(proc, "patient_mutation_peritoneal.tsv")
fwrite(res, out, sep = "\t")
cat("written:", out, "(", nrow(res), "genes )\n")
print(res[, .(gene, PM_mut, PM_n, nonPM_mut, nonPM_n, pct_PM, pct_nonPM,
              OR = round(OR, 3), lo = round(lo, 3), hi = round(hi, 3),
              p = signif(p, 4), fdr = signif(fdr, 4))])
