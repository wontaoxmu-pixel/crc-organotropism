# 修订 A2：MSK-MET 深化（RA-M2 / RB-m4 / RA-M4）
#  1) 1,147 转移灶样本 vs 患者数：伪重复核实；多灶患者 → 每患者单样本重跑样本级 APC Fisher
#  2) 排除卵巢转移的敏感性分析：样本级 + 患者级 APC/SMAD4/RNF43 × 腹膜（仅 Intra-Abdominal）
#  3) 患者级 APC×腹膜多变量 logistic：年龄、性别（全队列）；MSI/TMB/原发侧别（sample.txt 792 子集）
# 输出：mskmet/revision_*.tsv
set.seed(123)
suppressMessages({library(data.table); library(jsonlite)})
raw <- "data/raw/crc_organotropism"
dir.create("mskmet", showWarnings = FALSE)

crc_ids <- fromJSON(file.path(raw, "mskmet_crc_sample_ids.json"), simplifyVector = TRUE)
read_attr <- function(f) {
  x <- fromJSON(file.path(raw, f), simplifyVector = FALSE)
  rbindlist(lapply(x, function(r) data.table(sampleId = r$sampleId, patientId = r$patientId,
                                             value = r$value)))
}
site <- read_attr("mskmet_S_METASTATIC_SITE.json");  setnames(site, "value", "met_site")
stype <- read_attr("mskmet_S_SAMPLE_TYPE.json");     setnames(stype, "value", "sample_type")
clin <- merge(site, stype[, .(sampleId, sample_type)], by = "sampleId", all.x = TRUE)
clin <- clin[sampleId %in% crc_ids]

mut <- fromJSON(file.path(raw, "mskmet_crc_mutations.json"), simplifyVector = TRUE)
mut <- data.table(sampleId = mut$sampleId, patientId = mut$patientId,
                  gene = mut$gene$hugoGeneSymbol, mutType = mut$mutationType)
mut <- mut[mutType != "Silent"]

## ========== 1) 转移灶样本 → 患者：伪重复核实 ==========
met <- clin[sample_type == "Metastasis" & met_site %in% c("Liver", "Lung", "Intra-Abdominal", "Ovary")]
spp <- met[, .N, by = patientId]
cat("== 转移灶样本（4 器官标注）==\n")
cat("样本数:", nrow(met), " 患者数:", uniqueN(met$patientId), "\n")
cat("每患者样本数分布:\n"); print(table(spp$N))
multi <- spp[N > 1]
cat("多灶患者数:", nrow(multi), "；涉及样本:", sum(multi$N), "\n")
fwrite(spp[order(-N)], "mskmet/revision_samples_per_patient.tsv", sep = "\t")

# 每患者保留 1 个转移灶样本（按 sampleId 排序取第一个，确定性去重）
setorder(met, patientId, sampleId)
met1 <- met[, .SD[1], by = patientId]
cat("去重后样本级分析集:", nrow(met1), "样本 =", uniqueN(met1$patientId), "患者\n")

fisher_or <- function(case_mut, case_n, ctrl_mut, ctrl_n) {
  m2 <- matrix(c(case_mut, case_n - case_mut, ctrl_mut, ctrl_n - ctrl_mut), nrow = 2, byrow = TRUE)
  ft <- fisher.test(m2)
  list(OR = as.numeric(ft$estimate), lo = ft$conf.int[1], hi = ft$conf.int[2], p = ft$p.value)
}

# 样本级 APC Fisher（PM=IA+Ovary vs Liver）：全部样本 vs 每患者单样本
run_sample_fisher <- function(dt, genes, pm_def = c("IAO", "IA")) {
  pm_sites <- if (pm_def == "IAO") c("Intra-Abdominal", "Ovary") else "Intra-Abdominal"
  pm_ids <- dt[met_site %in% pm_sites]$sampleId
  lv_ids <- dt[met_site == "Liver"]$sampleId
  rbindlist(lapply(genes, function(g) {
    mgs <- unique(mut[sampleId %in% dt$sampleId & gene == g]$sampleId)
    r <- fisher_or(sum(pm_ids %in% mgs), length(pm_ids), sum(lv_ids %in% mgs), length(lv_ids))
    data.table(gene = g, PM_mut = sum(pm_ids %in% mgs), PM_n = length(pm_ids),
               Liver_mut = sum(lv_ids %in% mgs), Liver_n = length(lv_ids),
               OR = r$OR, CI_lo = r$lo, CI_hi = r$hi, p = r$p)
  }))
}

cat("\n== 样本级 APC：全部样本 vs 每患者单样本（PM=IA+Ovary vs Liver）==\n")
a_all <- run_sample_fisher(met, "APC", "IAO")[, `:=`(level = "sample_all", pm_def = "IA+Ovary")]
a_dd  <- run_sample_fisher(met1, "APC", "IAO")[, `:=`(level = "sample_dedup_1perPt", pm_def = "IA+Ovary")]
res_apc <- rbind(a_all, a_dd)
print(res_apc)
fwrite(res_apc, "mskmet/revision_apc_samplelevel_dedup.tsv", sep = "\t")

## ========== 2) 排除卵巢的敏感性分析（APC / SMAD4 / RNF43）==========
genes3 <- c("APC", "SMAD4", "RNF43")
cat("\n== 样本级：IA+Ovary vs 仅 IA（对照均为 Liver）==\n")
s_iao <- run_sample_fisher(met, genes3, "IAO")[, `:=`(level = "sample_all", pm_def = "IA+Ovary")]
s_ia  <- run_sample_fisher(met, genes3, "IA")[, `:=`(level = "sample_all", pm_def = "IA_only")]
sens_s <- rbind(s_iao, s_ia)
sens_s[, fdr_within_def := p.adjust(p, "BH"), by = pm_def]
print(sens_s)

# 患者级：hasPM 定义前后
pw_iao <- met[, .(hasPM = any(met_site %in% c("Intra-Abdominal", "Ovary"))), by = patientId]
pw_ia  <- met[, .(hasPM = any(met_site == "Intra-Abdominal")), by = patientId]
mut_p <- mut[, .(patientId, gene)] |> unique()
run_patient_fisher <- function(pw, genes) {
  rbindlist(lapply(genes, function(g) {
    ids <- unique(mut_p[gene == g]$patientId)
    a <- sum(pw$hasPM & pw$patientId %in% ids); b <- sum(pw$hasPM & !pw$patientId %in% ids)
    c1 <- sum(!pw$hasPM & pw$patientId %in% ids); d <- sum(!pw$hasPM & !pw$patientId %in% ids)
    r <- fisher_or(a, a + b, c1, c1 + d)
    data.table(gene = g, PM_mut = a, PM_n = a + b, nonPM_mut = c1, nonPM_n = c1 + d,
               OR = r$OR, CI_lo = r$lo, CI_hi = r$hi, p = r$p)
  }))
}
cat("\n== 患者级：IA+Ovary vs 仅 IA ==\n")
p_iao <- run_patient_fisher(pw_iao, genes3)[, `:=`(level = "patient", pm_def = "IA+Ovary")]
p_ia  <- run_patient_fisher(pw_ia, genes3)[, `:=`(level = "patient", pm_def = "IA_only")]
sens_p <- rbind(p_iao, p_ia)
sens_p[, fdr_within_def := p.adjust(p, "BH"), by = pm_def]
print(sens_p)
sens_all <- rbind(sens_s, sens_p, fill = TRUE)
fwrite(sens_all, "mskmet/revision_ovary_excluded_sensitivity.tsv", sep = "\t")

## ========== 3) 患者级 APC×腹膜 多变量 logistic ==========
phdr <- readLines(file.path(raw, "mskmet_patient.txt"), n = 1)
pcn <- strsplit(sub("^#", "", phdr), "\t")[[1]]
pat <- fread(file.path(raw, "mskmet_patient.txt"), skip = 4, header = FALSE, col.names = pcn)
pat <- pat[, .(patientId = `Patient Identifier`, sex = Sex,
               age = as.numeric(`Age at Sequencing`))]
apc_ids <- unique(mut_p[gene == "APC"]$patientId)
df <- merge(pw_iao, pat, by = "patientId")
df[, APC_mut := patientId %in% apc_ids]
df <- df[!is.na(age) & sex %in% c("Male", "Female")]
df[, sex_male := as.integer(sex == "Male")]
cat("\n== logistic 全队列（有年龄/性别者）:", nrow(df), "患者；PM", sum(df$hasPM), "==\n")
fit1 <- glm(hasPM ~ APC_mut + age + sex_male, data = df, family = binomial)
s1 <- summary(fit1)$coefficients
ci1 <- suppressMessages(confint(fit1))
tab1 <- data.table(model = "full: hasPM ~ APC + age + sex", term = rownames(s1),
                   n = nrow(df), OR = exp(s1[, 1]), CI_lo = exp(ci1[, 1]),
                   CI_hi = exp(ci1[, 2]), p = s1[, 4])
print(tab1)

# 扩展协变量（仅 sample.txt 覆盖的 792 样本子集）：MSI、TMB、原发侧别
shdr <- readLines(file.path(raw, "mskmet_sample.txt"), n = 1)
scn <- strsplit(sub("^#", "", shdr), "\t")[[1]]
sm <- fread(file.path(raw, "mskmet_sample.txt"), skip = 4, header = FALSE, col.names = scn)
sm <- sm[`Sample Identifier` %in% crc_ids,
         .(patientId = `Patient Identifier`, primary = `Primary Tumor Site`,
           msi = `MSI Type`, tmb = as.numeric(`TMB (nonsynonymous)`))]
sm <- unique(sm, by = "patientId")  # sample.txt 每患者 1 行（已核实 6531 行=6531 患者）
right <- c("Cecum", "Ascending Colon", "Transverse Colon")
sm[, side := fifelse(primary %in% right, "Right",
             fifelse(primary %in% c("Descending Colon", "Sigmoid Colon",
                                    "Rectosigmoid Colon", "Rectum"), "Left", NA_character_))]
df2 <- merge(df, sm, by = "patientId")
df2[, msi_unstable := as.integer(msi == "Instable")]
df2 <- df2[!is.na(tmb) & !is.na(side) & msi %in% c("Stable", "Instable")]
cat("\n== logistic 子集（sample.txt 覆盖，MSI/TMB/侧别可用）:", nrow(df2),
    "患者；PM", sum(df2$hasPM), "==\n")
cat("子集侧别分布:", table(df2$side), " MSI-H:", sum(df2$msi_unstable), "\n")
fit2 <- glm(hasPM ~ APC_mut + age + sex_male + msi_unstable + tmb + side,
            data = df2, family = binomial)
s2 <- summary(fit2)$coefficients
ci2 <- suppressMessages(confint(fit2))
tab2 <- data.table(model = "subset: + MSI + TMB + side", term = rownames(s2),
                   n = nrow(df2), OR = exp(s2[, 1]), CI_lo = exp(ci2[, 1]),
                   CI_hi = exp(ci2[, 2]), p = s2[, 4])
print(tab2)
fwrite(rbind(tab1, tab2), "mskmet/revision_apc_logistic.tsv", sep = "\t")

cat("\n完成。输出: mskmet/revision_samples_per_patient.tsv / revision_apc_samplelevel_dedup.tsv /",
    "revision_ovary_excluded_sensitivity.tsv / revision_apc_logistic.tsv\n")
