# Title page

**Title:** Metastatic organotropism in colorectal cancer: a peritoneal stromal–adipocytic imprint, an APC-depleted genomic background, and the limits of primary-tumor prediction

**Running title:** Niche imprint and APC depletion in CRC metastasis

**Authors:** [Author list to be completed]

**Affiliations:** [Affiliations to be completed]

**Corresponding author:** [Name, address, email to be completed]

# Abstract

**Background:** Colorectal cancer (CRC) metastasizes preferentially to the liver and peritoneum, yet the molecular basis of this organotropism remains incompletely defined, and whether the primary-tumor transcriptome encodes metastatic destination is unresolved.

**Methods:** We pooled patient-matched comparisons across five public gene-expression cohorts to derive a peritoneal-tropism gene set (PM_core) and a liver-tropism gene set (LM_core). Consensus molecular subtypes (CMS) were assigned in a multi-organ RNA-seq cohort, with per-patient modeling of subtype enrichment. Genomic associations of peritoneal versus liver spread were tested in 3,548 CRC samples from MSK-MET. Whether primary-tumor expression predicts metastatic site was evaluated by single-sample gene-set scoring, restricting validation to cohorts not used for gene-set derivation and comparing performance against size-matched random gene sets.

**Results:** Peritoneal metastases showed a stromal/adipocytic niche imprint (PM_core, 49 genes, including *FABP4* and *PI16*), with CMS4 enrichment confirmed as a positive control under patient-aware modeling (OR = 3.72, 95% CI 1.23–11.29, *P* = 0.020 versus primary tumors). Liver metastases carried a hepatocyte-associated program (LM_core, 57 genes). Peritoneal spread tracked an APC-depleted genomic background: *APC* mutations were depleted at the metastasis-sample level (64.5% vs 82.2%; OR = 0.39, FDR = 5.1 × 10⁻⁴), with a nominal *SMAD4* enrichment that did not survive multiple-testing correction. Primary-tumor expression of these gene sets could not be distinguished from size-matched random gene sets under leave-one-cohort-out validation at the current sample sizes. The liver-tropism signal was quantitatively consistent with substantial hepatocyte admixture.

**Conclusions:** Metastatic lesions retain organ-specific transcriptomic niche imprints, and peritoneal spread tracks an APC-depleted genomic background. The primary-tumor transcriptome was not reliably supported as a predictor of metastatic destination, and bulk liver-metastasis transcriptomes require routine control for hepatocyte contamination.

**Keywords:** colorectal cancer; peritoneal metastasis; liver metastasis; organotropism; consensus molecular subtypes; APC; hepatocyte contamination; transcriptomics

---
## Internal notes

### Section outline
1. Title：预审修订后定稿（蓝图 §2 候选 3 修订版）——「APC–SMAD4 genomic axis」已废弃，改为「APC-depleted genomic background」，保留"诚实阴性结果进标题"的防守性设计；running title 相应改为 "Niche imprint and APC depletion in CRC metastasis"。
2. Background：两句——器官趋向未明 + 原发灶能否预测转移部位未决（对应 P1/P3 的问题）。
3. Methods：patient-matched 比较汇集五队列 → 基因集；CMS 增加 per-patient modeling；MSK-MET 基因组；ssGSEA 限定非推导队列验证 + size-matched random gene sets 基线。只写回答问题所需的设计逻辑。
4. Results：P1（生态位印记 49/57 + CMS4 GLMM OR = 3.72 阳性对照）→ P2（APC 样本级 64.5% vs 82.2%，SMAD4 仅 nominal 且不过校正）→ P3（cannot be distinguished from size-matched random gene sets，克制句式）→ P4（quantitatively consistent with substantial hepatocyte admixture 一句）。
5. Conclusions：生态位印记 + APC-depleted 基因组背景为候选框架；原发灶预测 not reliably supported；肝转移 bulk 研究须控污染。
6. Keywords：8 个——SMAD4 移出、hepatocyte contamination 纳入（与主稿一致）。

### Assumptions or missing inputs
- SMAD4 口径升级：主稿摘要现已明确写出 "a nominal *SMAD4* enrichment that did not survive multiple-testing correction"（患者级 OR = 1.62, *P* = 0.048, FDR = 0.214；样本级 24.7% vs 17.0%, ns）——不再使用旧稿方向性表述 "enriched in the same peritoneal context"。
- APC 主口径改为样本级 64.5% vs 82.2%（腹膜 = 腹腔内 + 卵巢合并，*n* = 124 vs 肝 624）；58.0%（47/81）仅为排除卵巢的敏感性分析，不进摘要。
- CMS4 摘要口径改为患者感知 GLMM 主分析（OR = 3.72, 95% CI 1.23–11.29, *P* = 0.020 vs PT），并明确标注 "confirmed as a positive control"；旧稿样本级 Fisher（OR = 3.13）降为正文描述性结果，不进摘要。
- P3 摘要措辞锁定 "could not be distinguished from size-matched random gene sets ... at the current sample sizes"；P4 锁定 "quantitatively consistent with substantial hepatocyte admixture"。
- 49/57 基因集名称按蓝图 §8 锁定形式（PM_core / LM_core）；摘要中首次出现给全称。
- "five public gene-expression cohorts" 未逐一列 accession（摘要篇幅限制）；正文 Methods 详列。
- 作者列表/单位/通讯信息仍为占位——投稿包阶段补充。

### Claim-evidence map

| Claim | Evidence | Status |
|---|---|---|
| PM 生态位印记（stromal/adipocytic, PM_core 49 基因） | 蓝图 §4.2（meta 282 显著 ∩ 方向一致 → 49 基因；FABP4 β=3.43 等） | 已核实 |
| CMS4 富集为患者感知 GLMM 阳性对照（OR = 3.72, 95% CI 1.23–11.29, *P* = 0.020 vs PT） | 蓝图 §4.3 | 已核实 |
| LM_core 57 基因肝细胞程序 | 蓝图 §4.2（FGG β=6.93 等） | 已核实 |
| APC 缺失（64.5% vs 82.2%; OR = 0.39, FDR = 5.1×10⁻⁴，样本级主口径） | 蓝图 §4.5 | 已核实 |
| SMAD4 仅 nominal 富集、不过多重校正（*P* = 0.048, FDR = 0.214） | 蓝图 §4.5 | 已核实（措辞与主稿一致） |
| 原发灶预测不可与随机基因集区分（LOCO + 1,000 排列基线） | 蓝图 §4.4 | 已核实 |
| 肝趋向信号 quantitatively consistent with substantial hepatocyte admixture | 蓝图 §4.6 | 已核实 |
