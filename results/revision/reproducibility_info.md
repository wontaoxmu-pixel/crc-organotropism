# A11 — 可重复性信息收集（RA-m5 / RC-m4）

生成日期：2026-10-04。所有版本号均来自本机实测（`Rscript -e 'sessionInfo()'` 与 `packageVersion()`），非凭记忆填写。

---

## ① 软件环境：sessionInfo 摘要与 Methods 版本句

### 实测 sessionInfo 摘要

| 项目 | 实测值 |
|---|---|
| R 版本 | R version 4.2.0 (2022-04-22) |
| 平台 | x86_64-apple-darwin17.0 (64-bit) |
| 操作系统 | macOS 26.6.2 |
| Bioconductor | 3.15 |
| 随机种子 | set.seed(123)，全部脚本首行固定 |

### 实际用到的包版本（以 scripts/10–23 的 library() 调用为准，逐个用 packageVersion() 核实）

| 包 | 版本 | 用途（对应脚本） |
|---|---|---|
| data.table | 1.14.4 | 全脚本数据读写（10–19, 21, 22） |
| readxl | 1.4.3 | GEO 补充表解析（10） |
| limma | 3.52.2 | 配对差异表达（11, 12, 18） |
| metafor | 5.2.1 | 随机效应 meta 分析 REML（12, 17） |
| org.Hs.eg.db | 3.15.0 | Ensembl↔Symbol/Entrez 映射（11, 13, 16, 17, 18） |
| AnnotationDbi | 1.58.0 | 同上（mapIds） |
| CMSclassifier | 1.0.0 | CMS 分型 RF 模型（13） |
| hgu133a.db | 3.13.0 | GSE41258 芯片探针注释（14） |
| hgu133plus2.db | 3.13.0 | GSE41568 芯片探针注释（14） |
| jsonlite | 1.8.8 | MSK-MET cBioPortal JSON 解析（15, 19, 22） |
| survival | 3.5.7 | Cox 比例风险模型（15, 19） |
| ggplot2 | 3.4.4 | 图 1–4（20–23） |
| patchwork | 1.1.3 | 多面板拼版（20–23） |
| scales | 1.2.1 | 坐标轴格式化（20, 22, 23） |
| ggrepel | 0.9.4 | 图 2 标签防重叠（21） |
| grid | 4.2.0（R 自带） | 面板定位（21, 23） |

### 关于任务提示中点名、但脚本实际未调用的包（重要，防止 Methods 误写）

- **GSVA 1.44.5**：已安装但**未使用**。脚本 14/16/17 注释已写明：GSVA 1.44.5 与 matrixStats（实测 1.2.0）存在已知不兼容，故 ssGSEA 采用手工实现的 Barbie 2009 秩次 KS 游走统计（scripts/14 的 `ssgsea_manual()`）。Methods 第 55 行已按此表述，版本句中不应出现 GSVA。
- **GEOquery 2.64.2 / edgeR 3.38.4 / sva 3.44.0**：均已安装，但 scripts/10–23 的 library() 调用中未出现，下游分析未使用（数据下载为 GEO series matrix 文件 + cBioPortal JSON，批次效应用平台分层 meta 而非 ComBat）。版本句不列入，以免审稿人按"声明即使用"追问。

### Methods 可直接粘贴的英文版本句（替换 main_manuscript.md 第 67 行占位符）

> Analyses were performed in R 4.2.0 (2022-04-22; x86_64-apple-darwin17.0, Bioconductor 3.15) using data.table 1.14.4, readxl 1.4.3, limma 3.52.2, metafor 5.2.1, survival 3.5.7, jsonlite 1.8.8, CMSclassifier 1.0.0, org.Hs.eg.db 3.15.0, AnnotationDbi 1.58.0, hgu133a.db 3.13.0, and hgu133plus2.db 3.13.0; figures were produced with ggplot2 3.4.4, patchwork 1.1.3, scales 1.2.1, and ggrepel 0.9.4. Single-sample gene-set scoring used a custom implementation of the rank-based ssGSEA statistic of Barbie et al. because of GSVA version-compatibility constraints. All tests were two-sided, and false-discovery rates were controlled by the BH method unless otherwise stated. The random seed was fixed (set.seed(123)) in all stochastic steps.

（第 67 行第二处占位符 [Evidence needed: repository URL] 需用户提供仓库地址，见 ③。）

---

## ② scripts/ 目录全量清单（Code availability 底稿）

| 脚本 | 一句话功能 |
|---|---|
| `00_get_manifest.py` | 查询 GDC API 生成 TCGA-THCA 诊断切片 manifest（注：属并行课题的病理切片管线，与本 CRC 论文无关） |
| `01_download_slides.py` | 按 manifest 断点续传下载 TCGA-THCA 全切片图像（同上，非本论文） |
| `02_tile_qc.py` | 切片 20× 下 256×256 patch 切分 + 组织检出 QC（同上，非本论文） |
| `10_M1_crc_data_qc.R` | CRC 器官趋向性 GEO 队列元数据解析、表达矩阵对齐与 QC，产出 data/processed/crc_organotropism/ |
| `11_M2_layer1_GSE190609_de.R` | GSE190609 同患者多器官配对 limma-trend 差异分析（PM/LM/PT 三组对比） |
| `12_M2_layer2_meta.R` | 各队列（GSE50760/225182/41258/41568）limma 差异分析 + metafor REML 随机效应整合，产出 PM_core（49 基因）/ LM_core（57 基因） |
| `13_M2_CMS.R` | GSE190609 全 113 样本 CMSclassifier RF 分型 + CMS4×腹膜 Fisher 关联（含 minPosterior≥0.5 严格版） |
| `14_M3_layer3_ssgsea.R` | 五队列原发灶手工 ssGSEA 打分（Barbie 2009 秩次法）+ 表观 AUC/ROC（PM_core、LM_core 回推验证） |
| `15_M3_layer4_mskmet.R` | MSK-MET 3,548 样本基因组×转移器官 Fisher（样本级+患者级）与 Cox OS 分析 |
| `16_audit_purity_confound.R` | 趋向基因集与肿瘤纯度/TME 混淆量化（EPCAM+KRT8/18 上皮代理相关 + 基因细胞来源标注） |
| `17_fix_B_LOCO.R` | 留一队列（LOCO）无循环敏感性验证：芯片队列重推 LM_core、GSE190609 重推 PM_core，跨队列打分 |
| `18_fix_A_purity.R` | GSE190609 配对 DE 加入上皮含量协变量，量化纯度控制后器官差异保留率 |
| `19_fix_C_cox_burden.R` | Cox OS 加入转移负荷（受累器官数）协变量的敏感性模型 |
| `20_M4_fig1_design_cms.R` | 图 1：队列×样本类型瓷砖图 + CMS4 腹膜富集（600-dpi PNG + 矢量 PDF） |
| `21_M4_fig2_niche_meta.R` | 图 2：meta 效应量火山图 + PM/LM 趋向 top 基因森林图 |
| `22_M4_fig3_genomic_axis.R` | 图 3：MSK-MET 突变×器官 dumbbell、患者级富集森林图、Cox OS |
| `23_M4_fig4_negative_contamination.R` | 图 4：LOCO 阴性结果 + 纯度调整保留率 + 肝细胞标志物跨队列 log2FC |

### Code availability 建议文本（待用户补仓库 URL 后可用）

> All analysis code is available at [repository URL]. The pipeline comprises 14 R scripts (scripts/10–23) covering data QC (10), paired differential expression (11, 18), cross-cohort random-effects meta-analysis (12), CMS classification (13), single-sample gene-set scoring and LOCO sensitivity analysis (14, 17), MSK-MET genomic and survival analyses (15, 19), purity-confounding audits (16), and figure generation (20–23).

注意：scripts/00–02（Python, TCGA-THCA 切片管线）不属于本论文，Code availability 不应包含；如仓库为整个项目目录，建议论文中明确 "scripts relevant to this study are scripts/10–23"，或单独建仓。

---

## ③ manuscript/main_manuscript.md 占位符扫描结果

全文（198 行）共 8 行、9 个占位符实例（第 67 行含 2 处），逐条如下：

| # | 行号 | 占位符原文 | 需要用户提供的信息 |
|---|---|---|---|
| 1 | 7 | `[Author list to be completed]` | 作者列表（姓名、排序） |
| 2 | 9 | `[Affiliations to be completed]` | 各单位全称与编号对应 |
| 3 | 11 | `[Name, address, email to be completed]` | 通讯作者姓名、地址、邮箱 |
| 4 | 67（第 1 处） | `[Evidence needed: R version and key package versions]` | **本报告 ① 已解决**：直接粘贴上方英文版本句替换 |
| 5 | 67（第 2 处） | `[Evidence needed: repository URL]` | 代码仓库地址（GitHub/Zenodo 等；需用户先建仓或给 DOI） |
| 6 | 117 | `[Evidence needed: repository URL]`（Code availability 段） | 同上，与第 5 处共用同一 URL |
| 7 | 121 | `[To be completed]`（Funding） | 基金号与资助机构全称 |
| 8 | 123 | `[To be completed]`（Conflict of interest） | 利益冲突声明（无则写 "The authors declare no competing interests."） |
| 9 | 125 | `[To be completed]`（Author contributions） | 各作者贡献（建议按 CRediT 角色填写） |

另：未发现其他 `[Evidence needed: ...]`、TODO、FIXME、XXX、CITATION 等模式残留。Data availability（115 行）与 Acknowledgments（119 行）已成稿，无占位符。

---

## 验证记录

- `Rscript -e 'sessionInfo()'` 实测输出（全文存 /tmp/sessionInfo_A11.txt）：R 4.2.0 / x86_64-apple-darwin17.0 / macOS 26.6.2 / Bioc 3.15。
- 16 个包版本逐个经 `packageVersion()` 实测（见 ① 表）；GSVA/GEOquery/edgeR/sva 的"已装未用"判定基于 `grep library()/require()` 对 scripts/10–23 的全文扫描。
- 脚本功能说明全部取自各脚本首行注释（实测 head 读取），非推测。
- 占位符扫描用 ripgrep 对 main_manuscript.md 全模式匹配（`\[Evidence needed:`、`\[To be completed\]`、`\[Author/Affiliations/Name.*to be completed\]`、TODO、FIXME、XXX、CITATION），共命中 8 行 9 处。
