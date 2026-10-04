# Supplementary display items 状态总览（2026-10-04）

对应 `manuscript/main_manuscript.md` 的 "Supplementary display items (proposal)" 一节（第 184–203 行）。

## Supplementary Figures S1–S5：全部已生产（`supplementary/figures/`，PNG 600 dpi + 矢量 PDF）

生成脚本：`scripts/35_supplementary_figures.R`（S1 另需先跑 `scripts/13_M2_CMS.R` 落盘后验概率）。碰撞审计记录：`supplementary/figures/SupplementaryFigureS{1-5}.collision-audit.json`（S1/S3/S4/S5 PASS；S2 为 0 fail + 8 warn，warn 均为 Panel A 的 EMT 基因标签与散点填充部分重叠，属火山图标注的固有情况）。

| 图 | 内容 | 数据依赖 | 关键数字（脚本实算复核） |
|---|---|---|---|
| S1 | 严格 CMS 分型（minPosterior ≥ 0.5, n = 74）的 CMS×器官分布 + Fisher 完整结果 | `data/processed/crc_organotropism/GSE190609_CMS_posteriors.tsv`（新增，由修改后的 `scripts/13_M2_CMS.R` 生成）、`GSE190609_meta.tsv` | n = 74；OR = 2.63，95% CI 0.92–7.88，P = 0.060（fisher.test 实算复核一致）；CMS4：PM 22/41 (53.7%) vs 非 PM 10/33 (30.3%) |
| S2 | GSE190609 三联配对 DE 火山图（PM vs PT / LM vs PT / PM vs LM）+ EMT sanity-check 基因标注 | `data/processed/crc_organotropism/de/GSE190609_{PM_vs_PT,LM_vs_PT,PM_vs_LM}.tsv.gz` | 显著基因数（FDR < 0.05 且 \|log2FC\| ≥ 1）：861 / 534 / 1,053（487 上 + 566 下），脚本内 stopifnot 断言通过与正文一致；EMT 基因方向全部符合预期（SPARC/CDH2/FN1/FAP/COL3A1 上调，EPCAM/CDH1/MUC1 下调） |
| S3 | 平台分层 meta（PM 分层即单队列；LM 分 RNA-seq 层与芯片层）vs 合并 β 散点 | 合并 meta（只读）+ 各队列 DE 表；新落盘 `data/processed/crc_organotropism/meta/revision_platform_stratified_{PM,LM}.tsv`（不覆盖现有文件） | Spearman ρ：PM 层 GSE190609 = 0.88 / GSE225182 = 0.39；LM 层 RNA-seq = 0.99 / 芯片 = 0.95 |
| S4 | 各评分配置的 1,000 随机基因集 AUC 零分布 + 观察值 + 经验 P | `ssgsea/revision_permutation_auc.tsv`、`ssgsea/revision_permutation_null_distributions.tsv.gz`（5 主配置 × 1000 行） | 经验 P：表观 LM_core 0.031；表观 PM_core 0.302；LOCO 主行 0.223；受限子集 0.367；跨队列 0.637（观察 0.532 < 零均值 0.600） |
| S5 | MSK-MET CNA×器官（探索性，n = 177，9 基因扩增频率 + Fisher P） | `mskmet/cna_amp_by_organ.tsv` | 全部 P > 0.05，无可信器官特异信号 |

## Supplementary Tables S1–S10：全部已生成（`supplementary/tables/`）

| 表 | 文件 | 行数 | 状态 |
|---|---|---|---|
| S1 | `TableS1_PM_core_49genes.csv` | 49 | ✅ 就绪（SE 为 Wald 反推列，见 tables/README 注 3） |
| S2 | `TableS2_LM_core_57genes.csv` | 57 | ✅ 就绪（同上） |
| S3 | `TableS3_full_meta_output.csv` | 23,418 | ✅ 就绪 |
| S4 | `TableS4_MSKMET_samplelevel_mutation_by_organ.csv` | 21 | ✅ 就绪（部分列空白为源文件未保存，见 tables/README 注 4） |
| S5 | `TableS5_MSKMET_patientlevel_fisher.csv` | 8 | ⚠️ 缺 APC 行（源文件基因列表遗漏，且 FDR 为 8 基因 BH 而非主稿的 9 基因 BH）；投稿前须补跑 `scripts/15_M3_layer4_mskmet.R` |
| S6 | `TableS6_cox_models_full_coefficients.csv` | 16 | ✅ 就绪 |
| S7 | `TableS7_CMS4_GLMM_and_duplicateCorrelation.csv` | 26 | ✅ 就绪 |
| S8 | `TableS8_permutation_baseline_1000random_sets.csv` | 13 | ✅ 就绪 |
| S9 | `TableS9_MSKMET_logistic_models1_2.csv` | 11 | ✅ 就绪 |
| S10 | `TableS10_external_signature_overlap.csv` | 14 | ✅ 就绪 |

各表的逐字段来源与口径说明见 `supplementary/tables/README.md`。生成脚本：`scripts/33_build_supplementary_tables.R`。
