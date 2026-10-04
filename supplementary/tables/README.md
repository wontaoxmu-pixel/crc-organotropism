# Supplementary Tables S1–S10（生成日期 2026-10-04）

全部 CSV 为 UTF-8 编码、首行表头，由 `scripts/33_build_supplementary_tables.R` 从已落盘的结果文件整合生成，未新造任何数字。生成脚本运行方式：在项目根目录 `Rscript scripts/33_build_supplementary_tables.R`。

## 逐表说明与来源

| 表 | 文件 | 内容（一句话） | 来源文件 |
|---|---|---|---|
| S1 | `TableS1_PM_core_49genes.csv`（49 行） | PM_core 49 个腹膜趋向基因：per-cohort log2FC（GSE190609 PM vs PT、GSE225182 PM vs PT、GSE190609 PM vs LM）、meta β/SE/P/FDR、τ²/I²/Q、FE 敏感性、GSE190609 纯度调整后保留状态 | `data/processed/crc_organotropism/meta/peritoneal_tropism_genes.tsv`、`meta/revision_heterogeneity_per_gene_PM.tsv.gz`、`de/GSE190609_PM_vs_PT.tsv.gz`、`de/GSE225182_PM_vs_PT.tsv.gz`、`de/GSE190609_PM_vs_LM.tsv.gz`、`meta/PM_core_purity_retained.tsv` |
| S2 | `TableS2_LM_core_57genes.csv`（57 行） | LM_core 57 个肝趋向基因，字段同 S1（per-cohort log2FC 覆盖 GSE190609/GSE50760/GSE41258/GSE41568 四个 LM vs PT 队列），另含 GSE50760 纯度调整保留状态 | `meta/liver_tropism_genes.tsv`、`meta/revision_heterogeneity_per_gene_LM.tsv.gz`、`de/*_LM_vs_PT.tsv.gz`（4 个）、`de/GSE190609_PM_vs_LM.tsv.gz`、`meta/LM_core_purity_retained.tsv`、`results/revision/lmcore_gse50760_purity.tsv` |
| S3 | `TableS3_full_meta_output.csv`（23,418 行 = PM 14,085 + LM 9,333） | 全基因组 meta 输出：REML β/P/FDR、τ²、I²、Cochran Q 及 Q 的 P、固定效应 β/P/FDR、方向一致性、`program` 列区分 PM_vs_PT / LM_vs_PT | `meta/revision_heterogeneity_per_gene_PM.tsv.gz`、`meta/revision_heterogeneity_per_gene_LM.tsv.gz` |
| S4 | `TableS4_MSKMET_samplelevel_mutation_by_organ.csv`（21 行 = 18 基因主分析 + 3 基因敏感性） | MSK-MET 样本级突变×器官 Fisher 全表（腹膜 IA+Ovary 混合 vs 肝，18 驱动基因，BH FDR），附卵巢排除（仅 IA，n=81）敏感性分析（APC/SMAD4/RNF43） | `data/processed/crc_organotropism/mskmet/mutation_by_organ.tsv`、`mskmet/revision_ovary_excluded_sensitivity.tsv`（项目根 `mskmet/` 目录） |
| S5 | `TableS5_MSKMET_patientlevel_fisher.csv`（9 行） | 患者级 Fisher（124 例腹膜受累 vs 1,023 例无），9 基因（含 APC），含 OR、P、9 基因统一 BH FDR 及 95% CI | `data/processed/crc_organotropism/mskmet/patient_mutation_peritoneal.tsv`（2026-10-04 已由 `scripts/34_rev_patient_mutation_table.R` 从原始 JSON 重算为 9 基因版）；CI 由该文件 2×2 计数以 `fisher.test` 精确重算 |
| S6 | `TableS6_cox_models_full_coefficients.csv`（16 行） | 三个 Cox 模型全系数表：同子集（n=405，226 事件）不含/含转移负荷项（M0/M1）+ 全队列 M0（n=1,136，524 事件），含 HR、95% CI、z、P | `data/processed/crc_organotropism/mskmet/revision_cox_samesubset.tsv` |
| S7 | `TableS7_CMS4_GLMM_and_duplicateCorrelation.csv`（26 行） | CMS4×器官 GLMM（lme4，患者随机截距）全输出 + duplicateCorrelation 阻断敏感性 + 各对比样本/患者数 + 患者级检验 | `results/revision/cms_patient_level.tsv` |
| S8 | `TableS8_permutation_baseline_1000random_sets.csv`（13 行） | 置换基线：5 个评分配置（表观/LOCO 主行/LOCO 可检出基因敏感子集/跨队列）× main+分层对比，含真实 AUC、bootstrap 95% CI（2,000 次）、1,000 个随机基因集的零分布统计（均值/SD/分位数/极值）与经验 P | `data/processed/crc_organotropism/ssgsea/revision_permutation_auc.tsv` |
| S9 | `TableS9_MSKMET_logistic_models1_2.csv`（11 行） | 多变量 logistic 回归全系数：模型 1（腹膜 vs 肝转移 ~ APC + 年龄 + 性别，n=887）与模型 2（+MSI+TMB+原发侧别，n=302） | `mskmet/revision_apc_logistic.tsv`（项目根目录） |
| S10 | `TableS10_external_signature_overlap.csv`（14 行） | 外部签名重叠超几何检验：PM_core/LM_core × Laoukili 2022、Lenos 2022、Gelli 2023、Yu 2025，宇宙 14,085/9,333 主分析 + 20,000 基因敏感性 | `results/revision/signature_overlap.tsv` |

## 缺口与注意事项（如实声明）

1. **S1/S2 的 `se_wald_derived` 列**：中间结果文件未直接保存 meta SE；该列由 metafor 正态 Wald 关系 z = β/SE 从全精度存储的 `beta_re`、`pval_re` 代数反推（SE = |β| / |qnorm(1−P/2)|），非估计新值。
2. **S4 列空缺**：主分析 18 基因行中仅 APC/SMAD4/RNF43 有突变计数与 95% CI（来自敏感性文件回补），其余 15 基因源文件只存百分比与 OR/P/FDR；敏感性行（IA_only）无分器官百分比列。空白即源文件未保存该字段。
3. **S8 的 `tag`/对比命名沿用分析脚本原始命名**（含中文 tag 与缩写），列含义见 `ssgsea/revision_permutation_auc_report.md`。

## 修订记录

- 2026-10-04：S5 源文件 `patient_mutation_peritoneal.tsv` 原为 `scripts/15_M3_layer4_mskmet.R` 的 8 基因旧产出（缺 APC，FDR 为 8 基因 BH 口径）。已由 `scripts/34_rev_patient_mutation_table.R` 从原始 JSON 重算并覆盖写回为 9 基因版（含 APC，OR = 0.46, 95% CI 0.31–0.71, P = 3.2×10⁻⁴, FDR = 0.003），FDR 统一为 9 基因 BH（与主稿及 Figure 3B 口径一致，如 SMAD4 FDR = 0.214）。
