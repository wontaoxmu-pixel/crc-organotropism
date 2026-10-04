# M3 层次 4：MSK-MET 2021 CRC 队列基因组 × 转移器官互证 + OS
# 数据：3,548 CRC 样本（API 筛选）/ 14,220 条驱动突变（18 基因，3,519 样本）
#       S_*.json 逐样本临床属性（转移部位等）+ patient.txt（OS）
# 分析：
#   1) 转移灶样本 突变 × 器官 Fisher（Intra-Abdominal/Ovary vs Liver vs Lung），BH 校正
#   2) 患者级器官受累 × 驱动突变频率
#   3) Cox OS ~ 腹膜受累 + BRAF + 年龄
set.seed(123)
suppressMessages({library(data.table); library(jsonlite); library(survival)})
raw <- "data/raw/crc_organotropism"
proc <- "data/processed/crc_organotropism"
dir.create(file.path(proc, "mskmet"), showWarnings = FALSE)

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
cat("CRC 样本（有部位标注）:", nrow(clin), " 患者:", uniqueN(clin$patientId), "\n")
print(table(clin$sample_type))

mut <- fromJSON(file.path(raw, "mskmet_crc_mutations.json"), simplifyVector = TRUE)
mut <- data.table(sampleId = mut$sampleId, patientId = mut$patientId,
                  gene = mut$gene$hugoGeneSymbol, mutType = mut$mutationType)
# 仅保留非沉默突变
mut <- mut[mutType != "Silent"]
genes <- sort(unique(mut$gene))

## ---------- 1) 转移灶样本 突变 × 器官 ----------
met <- clin[sample_type == "Metastasis" & met_site %in% c("Liver", "Lung", "Intra-Abdominal", "Ovary")]
met[, peritoneal := met_site %in% c("Intra-Abdominal", "Ovary")]
cat("\n转移灶样本器官分布:\n"); print(table(met$met_site))
mut_s <- mut[sampleId %in% met$sampleId, .(mutated = TRUE), by = .(sampleId, gene)]
res1 <- rbindlist(lapply(genes, function(g) {
  mgs <- unique(mut_s[gene == g]$sampleId)
  tab <- sapply(c("Liver", "Lung", "Intra-Abdominal", "Ovary"), function(org) {
    ids <- met[met_site == org]$sampleId
    c(mut = sum(ids %in% mgs), wt = sum(!ids %in% mgs))
  })
  # 腹膜（IA+Ovary）vs 肝 Fisher
  pm_ids <- met[peritoneal == TRUE]$sampleId; lv_ids <- met[met_site == "Liver"]$sampleId
  m2 <- matrix(c(sum(pm_ids %in% mgs), length(pm_ids) - sum(pm_ids %in% mgs),
                 sum(lv_ids %in% mgs), length(lv_ids) - sum(lv_ids %in% mgs)),
               nrow = 2, byrow = TRUE)
  ft <- fisher.test(m2)
  data.table(gene = g,
             liver_pct = round(100 * tab[1, "Liver"] / sum(tab[, "Liver"]), 1),
             lung_pct = round(100 * tab[1, "Lung"] / sum(tab[, "Lung"]), 1),
             intraab_pct = round(100 * tab[1, "Intra-Abdominal"] / sum(tab[, "Intra-Abdominal"]), 1),
             ovary_pct = round(100 * tab[1, "Ovary"] / sum(tab[, "Ovary"]), 1),
             OR_PM_vs_Liver = round(as.numeric(ft$estimate), 2),
             p_PM_vs_Liver = signif(ft$p.value, 3))
}))
res1[, fdr := p.adjust(p_PM_vs_Liver, "BH")]
setorder(res1, p_PM_vs_Liver)
cat("\n== 突变 × 器官（腹膜 vs 肝，按 P 排序）==\n"); print(res1)
fwrite(res1, file.path(proc, "mskmet", "mutation_by_organ.tsv"), sep = "\t")

## ---------- 2) 患者级器官受累 ----------
pw <- clin[sample_type == "Metastasis",
           .(hasPM = any(met_site %in% c("Intra-Abdominal", "Ovary")),
             hasLiver = any(met_site == "Liver"),
             hasLung = any(met_site == "Lung")), by = patientId]
cat("\n患者级: 腹膜受累", sum(pw$hasPM), "/", nrow(pw),
    " 肝", sum(pw$hasLiver), " 肺", sum(pw$hasLung), "\n")
# 关键基因患者级：腹膜受累 vs 无腹膜 的突变率
mut_p <- mut[, .(mutated = TRUE), by = .(patientId, gene)] |> unique()
res2 <- rbindlist(lapply(c("BRAF", "KRAS", "NRAS", "TP53", "PIK3CA", "RNF43", "PTEN", "SMAD4"), function(g) {
  ids <- unique(mut_p[gene == g]$patientId)
  a <- sum(pw$hasPM & pw$patientId %in% ids); b <- sum(pw$hasPM & !pw$patientId %in% ids)
  c1 <- sum(!pw$hasPM & pw$patientId %in% ids); d <- sum(!pw$hasPM & !pw$patientId %in% ids)
  ft <- fisher.test(matrix(c(a, b, c1, d), nrow = 2, byrow = TRUE))
  data.table(gene = g, PM_mut = a, PM_n = a + b, nonPM_mut = c1, nonPM_n = c1 + d,
             pct_PM = round(100 * a / (a + b), 1), pct_nonPM = round(100 * c1 / (c1 + d), 1),
             OR = round(as.numeric(ft$estimate), 2), p = signif(ft$p.value, 3))
}))
res2[, fdr := p.adjust(p, "BH")]
cat("\n== 患者级 突变 × 腹膜受累 ==\n"); print(res2)
fwrite(res2, file.path(proc, "mskmet", "patient_mutation_peritoneal.tsv"), sep = "\t")

## ---------- 3) Cox OS ----------
phdr <- readLines(file.path(raw, "mskmet_patient.txt"), n = 1)
pcn <- strsplit(sub("^#", "", phdr), "\t")[[1]]
pat <- fread(file.path(raw, "mskmet_patient.txt"), skip = 4, header = FALSE, col.names = pcn)
pat <- pat[, .(patientId = `Patient Identifier`, os_m = `Overall Survival (Months)`,
               os_status = `Overall Survival Status`, age = `Age at Sequencing`)]
pat <- pat[!is.na(os_m) & os_status %in% c("0:LIVING", "1:DECEASED")]
pat[, os_m := as.numeric(os_m)]
pat[, age := as.numeric(age)]
pat <- pat[!is.na(os_m)]
pat[, dead := as.integer(os_status == "1:DECEASED")]
df <- merge(pw, pat, by = "patientId")
braf_ids <- unique(mut_p[gene == "BRAF"]$patientId)
df[, BRAF_mut := patientId %in% braf_ids]
cat("\nCox 队列:", nrow(df), "患者, 死亡", sum(df$dead), "\n")
fit <- coxph(Surv(os_m, dead) ~ hasPM + hasLiver + hasLung + BRAF_mut + age, data = df)
print(summary(fit))
s <- as.data.table(summary(fit)$coefficients, keep.rownames = "term")
ci <- as.data.table(summary(fit)$conf.int)
fwrite(cbind(s, ci[, .(`lower .95`, `upper .95`)]),
       file.path(proc, "mskmet", "cox_os.tsv"), sep = "\t")

## ---------- 4) CNA 提示性分析（仅 404 样本，探索性） ----------
cna <- fromJSON(file.path(raw, "mskmet_crc_cna.json"), simplifyVector = TRUE)
cna <- data.table(sampleId = cna$sampleId, gene = cna$gene$hugoGeneSymbol, alt = cna$alteration)
cna <- cna[sampleId %in% clin$sampleId]
cat("\nCNA 覆盖样本:", uniqueN(cna$sampleId), "（量太少，仅探索）\n")
amp <- cna[alt == 2, .(amp = TRUE), by = .(sampleId, gene)]
cna_met <- clin[sampleId %in% cna$sampleId & sample_type == "Metastasis" &
                met_site %in% c("Liver", "Intra-Abdominal", "Ovary", "Lung")]
res4 <- rbindlist(lapply(unique(amp$gene), function(g) {
  ids <- unique(amp[gene == g]$sampleId)
  pm_ids <- cna_met[met_site %in% c("Intra-Abdominal", "Ovary")]$sampleId
  lv_ids <- cna_met[met_site == "Liver"]$sampleId
  m2 <- matrix(c(sum(pm_ids %in% ids), length(pm_ids) - sum(pm_ids %in% ids),
                 sum(lv_ids %in% ids), length(lv_ids) - sum(lv_ids %in% ids)),
               nrow = 2, byrow = TRUE)
  if (any(m2 < 0)) return(NULL)
  ft <- fisher.test(m2)
  data.table(gene = g, PM_amp_pct = round(100 * m2[1, 1] / sum(m2[1, ]), 1),
             Liver_amp_pct = round(100 * m2[2, 1] / sum(m2[2, ]), 1),
             OR = round(as.numeric(ft$estimate), 2), p = signif(ft$p.value, 3))
}))
if (!is.null(res4) && nrow(res4)) { setorder(res4, p); print(res4)
  fwrite(res4, file.path(proc, "mskmet", "cna_amp_by_organ.tsv"), sep = "\t") }
