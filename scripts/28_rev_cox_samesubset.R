# 任务 A5：生存分析修订（RA-M5 / RB-M4）
# 1) 在负荷（Met Count）标注子集上，同一子集拟合 M0（不含 met_count）与 M1（含 met_count），
#    对比 hasPM HR 变化 —— 区分「负荷解释」与「功效损失（子集缩小）」
# 2) 确认 OS 时间起点定义（patient.txt 字段注释行），并量化测序距首次转移诊断的时间间隔
# 3) 探查肺转移保护性 HR：hasLung 与 met_count / 年龄 / 共受累模式的关系
# 输出: data/processed/crc_organotropism/mskmet/revision_cox_samesubset.tsv
#       data/processed/crc_organotropism/mskmet/revision_lung_exploration.tsv
set.seed(123)
suppressMessages({library(data.table); library(jsonlite); library(survival)})
raw <- "data/raw/crc_organotropism"
outdir <- file.path("data/processed/crc_organotropism", "mskmet")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

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

phdr <- readLines(file.path(raw, "mskmet_patient.txt"), n = 2)
pcn <- strsplit(sub("^#", "", phdr[1]), "\t")[[1]]
pdef <- strsplit(sub("^#", "", phdr[2]), "\t")[[1]]
pat <- fread(file.path(raw, "mskmet_patient.txt"), skip = 4, header = FALSE, col.names = pcn)
pat <- pat[, .(patientId = `Patient Identifier`, os_m = `Overall Survival (Months)`,
               os_status = `Overall Survival Status`, age = `Age at Sequencing`,
               age_mets = `Age at First Mets Dx`)]
pat[, os_m := as.numeric(os_m)][, age := as.numeric(age)][, age_mets := as.numeric(age_mets)]
pat <- pat[!is.na(os_m) & os_status %in% c("0:LIVING", "1:DECEASED")]
pat[, dead := as.integer(os_status == "1:DECEASED")]

df <- merge(pw, pat, by = "patientId")
df[, BRAF_mut := patientId %in% braf_ids]
cat("全队列 Cox 患者数:", nrow(df), " 死亡:", sum(df$dead), "\n")

## ---- OS 起点定义（任务③）----
cat("\n== patient.txt 字段官方定义（第 2 行注释） ==\n")
for (i in seq_along(pcn)) cat(sprintf("  %-28s : %s\n", pcn[i], pdef[i]))

# 测序距首次转移诊断的间隔（age 单位若为月则直接为月）
df[, seq_after_mets := age - age_mets]
cat("\n== 测序日期 - 首次转移诊断日期 间隔（age 字段差值，全队列有值者） ==\n")
gap <- df[!is.na(seq_after_mets)]$seq_after_mets
cat("n =", length(gap), " 中位 =", median(gap), " IQR =", quantile(gap, .25), "-", quantile(gap, .75),
    " 范围 =", min(gap), "to", max(gap), "\n")
cat("age 字段量级参考: age 中位 =", median(df$age, na.rm = TRUE), " age_mets 中位 =", median(df$age_mets, na.rm = TRUE), "\n")

## ---- Met Count 整理：核对患者内一致性并去重（19 号脚本直接 merge 会因多样本复制患者行）----
shdr <- readLines(file.path(raw, "mskmet_sample.txt"), n = 1)
scn <- strsplit(sub("^#", "", shdr), "\t")[[1]]
sm <- fread(file.path(raw, "mskmet_sample.txt"), skip = 4, header = FALSE, col.names = scn)
sm <- sm[, .(patientId = `Patient Identifier`, met_count = as.numeric(`Met Count`))]
sm <- sm[!is.na(met_count)]
chk <- sm[, .(n_val = uniqueN(met_count), n_smp = .N, mc = max(met_count)), by = patientId]
cat("\n== Met Count 患者内一致性 ==\n")
cat("有 Met Count 的患者:", nrow(chk), "; 多样本患者:", sum(chk$n_smp > 1),
    "; 患者内取值不一致者:", sum(chk$n_val > 1), "\n")
smu <- chk[, .(patientId, met_count = mc)]

df2 <- merge(df, smu, by = "patientId")
cat("负荷标注子集（去重后）:", nrow(df2), "患者, 死亡", sum(df2$dead), "\n")

## ---- 同子集双模型（任务②）----
f_M0 <- Surv(os_m, dead) ~ hasPM + hasLiver + hasLung + BRAF_mut + age
m0_sub <- coxph(f_M0, data = df2)                                   # 子集, 不含负荷
m1_sub <- coxph(update(f_M0, . ~ . + met_count), data = df2)          # 子集, 含负荷
m0_full <- coxph(f_M0, data = df)                                   # 全队列, 不含负荷（参照）

getrow <- function(fit, model) {
  s <- summary(fit)
  ci <- s$conf.int
  data.table(model = model, term = rownames(s$coefficients),
             coef = s$coefficients[, "coef"], HR = s$coefficients[, "exp(coef)"],
             se = s$coefficients[, "se(coef)"], z = s$coefficients[, "z"],
             p = s$coefficients[, "Pr(>|z|)"],
             HR_lo95 = ci[, "lower .95"], HR_hi95 = ci[, "upper .95"],
             n = s$n, events = s$nevent)
}
res <- rbind(getrow(m0_full, "M0_fullcohort_no_burden"),
             getrow(m0_sub, "M0_subset_no_burden"),
             getrow(m1_sub, "M1_subset_with_met_count"))

lrt <- anova(m0_sub, m1_sub)
cat("\n== 同子集 M0 vs M1：hasPM 对照 ==\n")
pm <- res[term == "hasPMTRUE"]
print(pm[, .(model, HR, HR_lo95, HR_hi95, se, p, n, events)])
cat(sprintf("hasPM logHR 变化: M0子集 %.4f -> M1子集 %.4f (delta = %.4f, 相对衰减 %.1f%%)\n",
            pm$coef[pm$model == "M0_subset_no_burden"], pm$coef[pm$model == "M1_subset_with_met_count"],
            pm$coef[pm$model == "M1_subset_with_met_count"] - pm$coef[pm$model == "M0_subset_no_burden"],
            100 * (1 - pm$coef[pm$model == "M1_subset_with_met_count"] / pm$coef[pm$model == "M0_subset_no_burden"])))
cat("LRT (M1 vs M0, 同子集): Chisq =", round(lrt$Chisq[2], 3), " df =", lrt$Df[2],
    " P =", format.pval(lrt$`Pr(>|Chi|)`[2], digits = 3), "\n")
cat("\n全队列 M0（参照, n =", m0_full$n, "）hasPM HR =",
    round(summary(m0_full)$coefficients["hasPMTRUE", "exp(coef)"], 4),
    " P =", format.pval(summary(m0_full)$coefficients["hasPMTRUE", "Pr(>|z|)"], digits = 3), "\n")

fwrite(res, file.path(outdir, "revision_cox_samesubset.tsv"), sep = "\t")

## ---- 肺转移保护性 HR 探查（任务④）----
cat("\n== hasLung 特征对比（负荷子集 df2） ==\n")
cmp <- df2[, .(n = .N,
               deaths = as.numeric(sum(dead)),
               met_count_med = as.numeric(median(met_count)), met_count_mean = round(mean(met_count), 2),
               n_site_med = as.numeric(median(n_site)),
               age_med = round(median(age), 1),
               pct_hasLiver = round(100 * mean(hasLiver), 1),
               pct_hasPM = round(100 * mean(hasPM), 1),
               pct_BRAF = round(100 * mean(BRAF_mut), 1)), by = hasLung]
print(cmp)
cat("Wilcoxon met_count ~ hasLung: P =",
    format.pval(wilcox.test(met_count ~ hasLung, data = df2)$p.value, digits = 3), "\n")
cat("Wilcoxon age ~ hasLung: P =",
    format.pval(wilcox.test(age ~ hasLung, data = df2)$p.value, digits = 3), "\n")
cat("Fisher hasLung x hasLiver: P =",
    format.pval(fisher.test(df2$hasLung, df2$hasLiver)$p.value, digits = 3),
    " OR =", round(fisher.test(df2$hasLung, df2$hasLiver)$estimate, 3), "\n")
cat("Fisher hasLung x hasPM: P =",
    format.pval(fisher.test(df2$hasLung, df2$hasPM)$p.value, digits = 3),
    " OR =", round(fisher.test(df2$hasLung, df2$hasPM)$estimate, 3), "\n")

# 共受累模式
df2[, lung_pattern := fifelse(!hasLung, "no_lung",
                      fifelse(!hasLiver & !hasPM & n_site == 1, "lung_only",
                      fifelse(!hasLiver & !hasPM, "lung_other_no_liverPM",
                      "lung_with_liver_or_PM")))]
pat_tab <- df2[, .(n = .N, deaths = as.numeric(sum(dead)), os_med = round(median(os_m), 1),
                   met_count_med = as.numeric(median(met_count))), by = lung_pattern]
cat("\n== 肺受累模式分层 ==\n"); print(pat_tab)

# 分层/替代参数化：lung_only vs 其他, 及低负荷层内 hasLung 效应
m_lungonly <- coxph(Surv(os_m, dead) ~ I(lung_pattern == "lung_only") + hasPM + hasLiver +
                      BRAF_mut + age + met_count, data = df2)
cat("\n== 替代模型：lung_only（孤立肺转移）指示变量 ==\n")
print(round(summary(m_lungonly)$coefficients[, c("exp(coef)", "Pr(>|z|)")], 4))

df2_low <- df2[met_count <= median(met_count)]
m_low <- coxph(Surv(os_m, dead) ~ hasPM + hasLiver + hasLung + BRAF_mut + age, data = df2_low)
cat("\n== 低负荷层（met_count <= 中位", median(df2$met_count), ", n =", nrow(df2_low), "）==\n")
print(round(summary(m_low)$coefficients[, c("exp(coef)", "Pr(>|z|)")], 4))

lung_out <- list(cmp = cmp, pattern = pat_tab)
fwrite(cmp, file.path(outdir, "revision_lung_exploration.tsv"), sep = "\t")
fwrite(pat_tab, file.path(outdir, "revision_lung_patterns.tsv"), sep = "\t")

# 全队列 hasLung 共受累交叉表（验证子集中 OR=0 是否为子集现象）
cat("\n== 全队列（df, n =", nrow(df), "）hasLung 共受累 ==\n")
print(df[, .(n = .N, deaths = as.numeric(sum(dead)), os_med = round(median(os_m), 1),
             pct_hasLiver = round(100 * mean(hasLiver), 1),
             pct_hasPM = round(100 * mean(hasPM), 1),
             n_site_med = as.numeric(median(n_site))), by = hasLung])
cat("Fisher hasLung x hasLiver (全队列): P =",
    format.pval(fisher.test(df$hasLung, df$hasLiver)$p.value, digits = 3),
    " OR =", round(fisher.test(df$hasLung, df$hasLiver)$estimate, 3), "\n")
df[, lung_pattern_full := fifelse(!hasLung, "no_lung",
                          fifelse(!hasLiver & !hasPM & n_site == 1, "lung_only",
                          fifelse(!hasLiver & !hasPM, "lung_other_no_liverPM",
                          "lung_with_liver_or_PM")))]
print(df[, .(n = .N, deaths = as.numeric(sum(dead)), os_med = round(median(os_m), 1)),
         by = lung_pattern_full])

cat("\n输出已写入:", file.path(outdir, "revision_cox_samesubset.tsv"), "\n")
