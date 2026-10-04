# 修复 C：Cox OS 加入转移负荷（每位患者受累器官数）协变量
set.seed(123)
suppressMessages({library(data.table); library(jsonlite); library(survival)})
raw <- "data/raw/crc_organotropism"
proc <- "data/processed/crc_organotropism"

crc_ids <- fromJSON(file.path(raw, "mskmet_crc_sample_ids.json"), simplifyVector = TRUE)
read_attr <- function(f) {
  x <- fromJSON(file.path(raw, f), simplifyVector = FALSE)
  rbindlist(lapply(x, function(r) data.table(sampleId = r$sampleId, patientId = r$patientId,
                                             value = r$value)))
}
site <- read_attr("mskmet_S_METASTATIC_SITE.json")
setnames(site, "value", "met_site")
clin <- site[sampleId %in% crc_ids]

mut <- fromJSON(file.path(raw, "mskmet_crc_mutations.json"), simplifyVector = TRUE)
mut <- data.table(patientId = mut$patientId, gene = mut$gene$hugoGeneSymbol,
                  mutType = mut$mutationType)[mutType != "Silent"]
braf_ids <- unique(mut[gene == "BRAF"]$patientId)

pw <- clin[, .(hasPM = any(met_site %in% c("Intra-Abdominal", "Ovary")),
               hasLiver = any(met_site == "Liver"),
               hasLung = any(met_site == "Lung"),
               n_site = uniqueN(met_site[met_site != ""])), by = patientId]

phdr <- readLines(file.path(raw, "mskmet_patient.txt"), n = 1)
pcn <- strsplit(sub("^#", "", phdr), "\t")[[1]]
pat <- fread(file.path(raw, "mskmet_patient.txt"), skip = 4, header = FALSE, col.names = pcn)
pat <- pat[, .(patientId = `Patient Identifier`, os_m = `Overall Survival (Months)`,
               os_status = `Overall Survival Status`, age = `Age at Sequencing`)]
pat[, os_m := as.numeric(os_m)][, age := as.numeric(age)]
pat <- pat[!is.na(os_m) & os_status %in% c("0:LIVING", "1:DECEASED")]
pat[, dead := as.integer(os_status == "1:DECEASED")]

df <- merge(pw, pat, by = "patientId")
df[, BRAF_mut := patientId %in% braf_ids]
cat("Cox 队列:", nrow(df), "患者, 死亡", sum(df$dead), "\n")

fit0 <- coxph(Surv(os_m, dead) ~ hasPM + hasLiver + hasLung + BRAF_mut + age, data = df)
cat("\n== 主模型（hasLiver/hasLung 已部分控制共受累）==\n")
print(round(summary(fit0)$coefficients[, c("exp(coef)", "Pr(>|z|)")], 4))

# 负荷敏感性：sample.txt 快照覆盖的子集（792 样本）用 Met Count
shdr <- readLines(file.path(raw, "mskmet_sample.txt"), n = 1)
scn <- strsplit(sub("^#", "", shdr), "\t")[[1]]
sm <- fread(file.path(raw, "mskmet_sample.txt"), skip = 4, header = FALSE, col.names = scn)
sm <- sm[, .(patientId = `Patient Identifier`, met_count = as.numeric(`Met Count`))]
df2 <- merge(df, sm, by = "patientId")
cat("\n负荷敏感性子集:", nrow(df2), "患者\n")
fit1 <- coxph(Surv(os_m, dead) ~ hasPM + hasLiver + hasLung + BRAF_mut + age + met_count, data = df2)
cat("== 含 Met Count（子集敏感性）==\n")
print(round(summary(fit1)$coefficients[, c("exp(coef)", "Pr(>|z|)")], 4))
s <- as.data.table(summary(fit1)$coefficients, keep.rownames = "term")
ci <- as.data.table(summary(fit1)$conf.int)
fwrite(cbind(s, ci[, .(`lower .95`, `upper .95`)]),
       file.path(proc, "mskmet", "cox_os_with_burden.tsv"), sep = "\t")
