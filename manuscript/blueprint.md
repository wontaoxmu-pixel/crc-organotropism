# Manuscript Blueprint — CRC 转移器官趋向性（选题 10）

> 本文件是全稿唯一事实来源（single source of truth）。所有 section 起草者只能使用本文件中的数字与结论；缺事实写占位符 `[Evidence needed: ...]`，禁止编造。内部工作文件，中文注释+英文术语。
> 依据规范：`~/.claude/rules/03_科研写作分析投稿项目级规范.md`（3.2.0 攻防平衡、3.2.4 参考文献、第二章统计格式）。
> **2026-10-04 预审修订**：主稿已按预审稿报告全量修订（GLMM 主口径、APC 主口径改 64.5% vs 82.2%、SMAD4 降为 nominal、同子集 Cox、随机基因集排列基线、污染归因分析、28 条参考文献新编号）；本蓝图已同步，「APC–SMAD4 genomic axis」表述全稿废弃。

## 0. 一句话论点（one-sentence argument）

In colorectal cancer (CRC), we show that metastatic lesions retain organ-specific transcriptomic niche imprints — a stromal/adipocytic program in peritoneal metastases and a hepatocyte-derived signal in liver metastases — and that peritoneal spread tracks an APC-depleted genomic background (with a nominal SMAD4 enrichment that did not survive multiple-testing correction), using paired multi-organ RNA-seq meta-analysis (5 GEO cohorts) and a 3,548-sample clinical-genomic cohort (MSK-MET), while demonstrating that primary-tumor prediction of metastatic destination cannot be distinguished from size-matched random gene sets at the current sample sizes.

## 1. 目标期刊与篇幅

- 首投候选：Clinical Cancer Research（CCR, 2025 JIF 10.9, Q1）或 British Journal of Cancer（7.8, Q1）；保底 International Journal of Cancer（4.7, Q1）。
- journal 轴 = generic（非 Nature 系）。
- Abstract ≤ 250 词（structured: Background / Methods / Results / Conclusions）。
- 主文（Intro+Methods+Results+Discussion）目标 3,800–4,500 词。
- 语言：美式英语（"tumor"）。注意：四张主图图内文字目前用 "tumour"，投稿包阶段统一重导，本稿文字一律美式。
- 语态/时态：Methods/Results 被动+过去时；已确立事实用现在时。

## 2. 题目候选（禁止 robust/novel/comprehensive/superior）

1. Organ-specific niche imprints and the genomic background of peritoneal versus liver metastasis in colorectal cancer
2. Transcriptomic niche imprint and an APC-depleted genomic background distinguish peritoneal from liver metastasis in colorectal cancer
3. Metastatic organotropism in colorectal cancer: a peritoneal stromal–adipocytic imprint, an APC-depleted genomic background, and the limits of primary-tumor prediction

**已定稿：候选 3**（主稿采用；running title: Niche imprint and APC depletion in CRC metastasis）。预审修订后 SMAD4 降为 nominal 趋势，候选 2/3 中的「APC–SMAD4 genomic axis」已同步改为「APC-depleted genomic background」；候选 3 保留"诚实阴性结果进标题"的防守性设计。

## 3. 四支柱 claim 链 ↔ 图 ↔ Results 小节

| 支柱 | Claim（校准后动词） | 证据 | 图 | Results 小节 |
|---|---|---|---|---|
| P1 | Peritoneal metastases **show** a stromal/adipocytic niche imprint, corroborated by CMS4 enrichment under patient-aware modeling (GLMM) — **treated as a positive control**, not a novel finding | 配对 DE + meta + CMS GLMM | Fig 1, Fig 2 | §1 队列与CMS、§2 生态位 |
| P2 | Peritoneal spread **is associated with** an APC-depleted genomic background (nominal SMAD4 enrichment not surviving multiple-testing correction); the serrated/WNT-ligand framework is retained **only as literature context** — not supported by our genomic data | MSK-MET Fisher + logistic + Cox | Fig 3 | §3 基因组背景 |
| P3 | Primary-tumor prediction of metastatic site **is not reliably supported** — no scoring configuration **can be distinguished from size-matched random gene sets at the current sample sizes** (beyond a marginal apparent LM_core signal) | ssGSEA 表观 vs LOCO + 1,000 随机基因集排列基线 | Fig 4A | §4 阴性结果 |
| P4 | Liver-met "liver-tropism" signal **is quantitatively consistent with** substantial hepatocyte admixture | 纯度校正 + 肝标志物 + 定量归因 | Fig 2C, Fig 4B/C | §4 污染警示 |

## 4. 已核实数字库（只准用这些）

### 4.1 队列（Fig 1A 数值以此为准）
- GSE190609（RNA-seq RPKM）：113 样本 / 12 患者：35 原发区域（PT）、59 腹膜（PM）、6 肝（LM）、13 其他转移（淋巴结 11 + 卵巢 2）。
- GSE50760（RNA-seq FPKM）：54 样本 / 18 患者 ×（正常 18 / PT 18 / LM 18）。
- GSE225182（RNA-seq log2CPM）：33 样本（PT 20 / PM 7 / 正常 6）；配对 meta 用 7 对。
- GSE41258（芯片）：PT 186 / LM 47 / 肺转移 20 / 正常 54。
- GSE41568（芯片）：133 样本：LM 80 / 肺 8 / PM 1 / 网膜 4 / 腹壁 1 / 未标注 39。（注意：旧文档写 79 肝，实为 80，已核实闭合）
- MSK-MET CRC 亚组：3,548 样本；**1,147 器官标注转移灶样本 = 1,147 个独立患者（严格一样本一患者，样本级基因组分析无患者内伪重复）**；Cox 可评估 n=1,136（524 events）；负荷标注子集 n=405（226 events）；患者级腹膜受累 n=124 vs 未受累 1,023；样本级腹膜（腹腔内+卵巢合并）n=124 vs 肝 n=624；测序距首次转移诊断中位 0.9 年（IQR 0.2–2.5；n=1,125）；146 例肺受累患者均为孤立肺转移（抽样选择，勿解释为生物学保护）。
- GSE41258/41568 因配对结构不成立仅作参考队列；GSE190609 全部 12 患者均有腹膜转移（无纯肝对照）。

### 4.2 层次 1–2：配对 DE + meta（Fig 2）
- GSE190609 配对 limma（患者阻断）：PM vs PT 861、LM vs PT 534、PM vs LM 1,053（上调 487/下调 566），FDR<0.05 且 |log2FC|≥1（仅 FDR 时为 3,046/553/1,233，S2 图口径为双阈值）。EMT sanity：PM 中 *SPARC/CDH2/FN1/FAP/COL3A1* 上调、*EPCAM/CDH1/MUC1* 下调。
- 随机效应 meta（metafor REML）：
  - PM 程序（GSE190609+GSE225182，14,085 共同基因）：282 个 FDR<0.05 且 |β|≥1，280（99.3%）方向一致。
  - LM 程序（GSE190609/50760/41258/41568，9,333 共同基因）：146 显著，145 方向一致。
- 器官趋向基因集（meta 显著且方向一致 ∩ GSE190609 PM_vs_LM 同向）：
  - **PM_core 49 个**（间质/脂肪谱系）：*FABP4* β=3.43 FDR<0.001；*CD36* 1.90；*PLAC9* 1.55；*RBP7* 1.52；*ADAMTS15* 1.55；*PI16* 2.35；*LPL* 1.28（均 FDR<0.001）；其余 *CADM3/SRPX/MEOX2* 等。
  - **LM_core 57 个**（肝功能/急性期蛋白）：*FGG* β=6.93 FDR=0.003；*HP* 5.84/0.001；*APOC3* 5.79/0.004；*GC* 5.75/0.020；*APOA1* 5.17/0.050；*F9* 3.71/0.045；*SERPINA1* 2.22/0.014；其余 *FGA/CPB2/ORM1/C9* 等。
- 肝程序在 PM vs PT 对比中基本不显著（57 中仅少数沉底）→ 两器官程序互不串扰，特异性证据。
- **异质性（预审修订新增）**：PM 程序 median τ² = 0（64.2% 基因）、median I² = 0%、Q 显著 7.1%；PM_core 内 median τ² = 0、median I² = 0%、Q 仅 *KCNA1* 1 个显著。LM 程序 median I² = 96.9%，57 个 LM_core 中 54 个 Q 显著。固定 vs 随机效应 β 相关：PM 0.988（近一致）/ LM 0.766（背离）。
- **k = 2 限定**：腹膜 meta 仅两个 RNA-seq 队列，REML τ² 估计不精确 → 异质性统计（τ²/I²/Q）仅作描述性解释；对策 = 固定效应敏感性 + LOCO（k=2 时退化为单队列重分析）+ 符号一致准入（单队列不能独自把基因送进 PM_core）。
- **PM_core LOCO（层次 2 内）**：方向 100% 跨队列一致；GSE225182 单队列 38/49 基因 |log2FC| ≥ 1（median 1.264 vs GSE190609 1.426）；0/49 达 FDR 显著 = 7 对功效所限。
- **外部签名重叠（hypergeometric；universe 14,085/9,333；20,000 基因敏感性结论不变）**：Laoukili et al.（[6]，366 腹膜上调基因）∩ PM_core = 24/49（25.1-fold，*P* = 1.4×10⁻²⁸，FDR = 2.0×10⁻²⁷）；Gelli et al.（[12]，61 基因签名）仅 14/40 在可检出 universe 内 → 不可评估（功效不足，非不重叠证据）；Lenos et al.（[5]，855 基因 PM_high）重叠 0；Yu et al.（[28]，805 基因单细胞程序）重叠 1（*FABP4*）。→ PM_core 为同一腹膜信号的跨队列 meta 汇聚核心（25/49 未被 Laoukili 单队列分析捕获），且为首个由三重汇聚（配对跨队列 meta + CMS4 阳性对照 + 外部签名富集）支持的腹膜核心集。

### 4.3 CMS 阳性对照（Fig 1B/C）
- CMSclassifier RF（272/273 模型基因覆盖）；nearestCMS 全 113 样本。
- **主口径（预审修订升级为患者感知 GLMM）**：GLMM（binomial，per-patient random intercept (1|patient)，lme4 1.1.30，收敛）：PM vs PT OR = 3.72（95% CI 1.23–11.29），*P* = 0.020（94 样本 / 12 患者）；PM vs LM OR = 32.57（0.84–1,268.63），*P* = 0.062（65 样本 / 12 患者；仅 6 个 LM 样本、来自 6 个不同患者 → 估计不稳定，须注明）；LM vs PT 对比 41 样本 / 12 患者。
- 患者级：9/12 患者至少有一个 CMS4 PM 样本；患者级平均 CMS4 比例 PM 52.8% vs LM 16.7%。样本级：CMS4 占 PM 47.5%（28/59）vs LM 16.7%（1/6）。
- 样本级 Fisher（忽略患者内相关）**仅作 descriptive**：OR = 3.13（95% CI 1.30–7.89），*P* = 0.006。
- 敏感性：duplicateCorrelation 分析与主分析方向一致（Jaccard = 0.904）；严格版（minPosterior≥0.5，n=74）方向保留：OR = 2.63（0.92–7.88），*P* = 0.060。
- **定位：阳性对照（positive control）而非新发现**——在患者感知建模下重导出同一队列既往报道的 CMS4–腹膜关联（Laoukili 2022 [6]，即 GSE190609 源文献）及独立外部系列（Lenos 2022 [5]），证明流程能恢复已知信号。（旧稿 "Saris 2025" 引用已证伪删除，见 §9。）

### 4.4 层次 3：原发灶回推 → 诚实阴性（Fig 4A）
- 表观（五队列 293 PT，手工 ssGSEA 秩次法）：LM_core 区分肝转移患者 vs 腹膜患者原发灶 AUC = 0.933（bootstrap 95% CI 0.845–0.996，*P* = 4.5×10⁻⁸）；PM_core 反向 AUC = 0.719（95% CI 0.560–0.859，*P* = 0.006）。
- LOCO 敏感性：LM_core 纯芯片推导（21 基因）→ RNA-seq PT：AUC = 0.736（95% CI 0.598–0.863，*P* = 0.003）仍 nominal 显著，但**限制为两组均可检出基因后** AUC = 0.596（*P* = 0.229）；PM_core 仅 GSE190609 推导（91 基因）→ 跨队列验证 AUC = 0.532（*P* = 0.753）。
- **随机基因集排列基线（预审修订新增）**：每个配置评分 1,000 个 size- 与 detection-rate–matched 随机基因集（bootstrap 95% CI，2,000 resamples）。仅表观 LM_core 勉强超过零假设（随机集均值 AUC = 0.489；empirical *P* = 0.031）；其余配置均未超过：表观 PM_core *P* = 0.302；LOCO 主行 *P* = 0.223；LOCO 限制子集 *P* = 0.367；cross-cohort 配置的随机集均值本身 = 0.600（队列间批次结构抬高零基线）。
- 结论表述（精确，不可夸大）：表观区分度由循环论证（基因集含同队列信息）与跨平台秩次偏移驱动；**no configuration can be distinguished from size-matched random gene sets at the current sample sizes（除表观 LM_core 行勉强显著外）→ primary-tumor prediction is not reliably supported — rather than disproven — at the current sample sizes**。**注意 LOCO 主行 0.736 仍 nominal 显著，写作必须写明结论依赖可检出基因限制子集 + 排列基线这一限定**。
- 方法学备注（预审修订新增）：cross-sample range normalization 为保序仿射变换，不影响秩次判别（严格单样本 vs 跨样本归一实现：Spearman ρ = 1.000、ΔAUC = 0，六种评分配置一致）；独立 GSEApy 1.3.1（Python）实现方向一致（LM_core AUC 0.903 / PM_core 0.716 vs 主实现 0.933 / 0.719）。

### 4.5 层次 4：MSK-MET 基因组背景 + 生存（Fig 3）
- **样本级主口径（预审修订改）**：腹膜 = 腹腔内 + 卵巢**合并**（n = 124）vs 肝（n = 624），18 driver 基因 BH；每个转移灶样本对应独立患者（无伪重复）。***APC* 突变 64.5% vs 82.2%，OR = 0.394，*P* = 2.86×10⁻⁵，FDR = 5.1×10⁻⁴**。*BRAF* 7.4% 不富集（FDR = 0.298）；*SMAD4* 24.7% vs 17.0%（ns）；*KRAS* 56.8% vs 41.3%（ns，样本级）。
- **卵巢排除敏感性（intra-abdominal only，n = 81）**：*APC* 47/81（**58.0% —— 仅限此敏感性分析，不作主口径**），OR = 0.300（95% CI 0.179–0.504），*P* = 2.98×10⁻⁶，FDR = 8.94×10⁻⁶（关联更强）。同一分析中 *SMAD4* 不显著（OR = 1.601，*P* = 0.092）；*RNF43* 富集由卵巢样本驱动（排除后 *P* = 0.121，ns）。
- **多变量 logistic（结局 = 腹膜 vs 肝转移；预审修订新增）**：model 1（校正 age+sex；n = 887，122 腹膜）：*APC* OR = 0.379（95% CI 0.249–0.584），*P* = 8.17×10⁻⁶（独立显著、与单变量基本不变）。model 2（小子集再校正 MSI/TMB/原发部位；n = 302）：方向保留但功效不足（*APC* OR = 0.544，*P* = 0.140）；右半结肠独立富集腹膜扩散（OR = 2.56，95% CI 1.18–5.49，*P* = 0.016）；MSI-high 仅 2 例不可评估。MSI/TMB/sidedness 注释来自覆盖不全的样本级快照，故全队列不可得。
- 患者级（腹膜受累 n=124 vs 未受累 n=1,023，9 基因 BH）：*APC* OR = 0.46（0.31–0.71），*P* = 3.2×10⁻⁴，FDR = 0.003；***SMAD4* 仅 nominal**：OR = 1.62（0.98–2.60），*P* = 0.048，**FDR = 0.214 不过多重校正**；*RNF43* 趋势 8.1% vs 5.0%，OR = 1.67（0.74–3.45），*P* = 0.142（ns）。
- **口径警示**：写作必须区分样本级（OR = 0.394）与患者级（OR = 0.46）两个 APC 数值；58.0% 仅出现在卵巢排除敏感性语境；SMAD4 一律写 "nominal enrichment that did not survive multiple-testing correction"；禁用「APC–SMAD4 axis」「SMAD4-enriched」表述。
- Cox OS：**OS_MONTHS 自测序日期起算；测序距首次转移诊断中位 0.9 年（IQR 0.2–2.5；n=1,125）→ 左截断（immortal-time bias），不可外推至诊断时起算**。全队列多变量模型（不含负荷项；n=1,136；524 events；协变量 = 腹膜/肝/肺受累 + *BRAF* + 测序年龄）：腹膜受累 HR = 1.42（1.05–1.92），*P* = 0.022；*BRAF* HR = 1.72（1.29–2.31），*P* < 0.001；肺受累 HR = 0.60（0.43–0.83），*P* = 0.002。
- **同子集比较（预审修订新口径；n=405；226 events；含/不含负荷项拟合同一批样本）**：无负荷项腹膜 HR = 1.38（0.86–2.22），*P* = 0.181（子集中已不显著）→ 加转移部位数后 HR = 1.49（0.92–2.39），*P* = 0.102。**点估计未降（log HR 反升）→ 失显著反映小子集功效损失，而非负荷解释预后信号**。met_count HR = 1.08/部位（1.05–1.11），*P* = 4.0×10⁻⁸；*BRAF* 负荷校正后 HR = 1.76（1.16–2.68），*P* = 0.008（稳健）。
- 结论口径：腹膜–OS 关联在全队列显著、在负荷标注子集中因功效损失失显著（同子集证据）；BRAF 预后作用稳健。**旧「部分由转移负荷解释」口径废弃**（见 §6）。
- CNA 分析仅 177 样本，无可信信号 → 仅探索性，进 SI。
- MSK-MET 腹膜 BRAF 7.4% vs GSE190609 32% 不一致 → Discussion 诚实处理（队列选择差异），*BRAF* 不纳入腹膜基因组背景。

### 4.6 纯度控制 + 污染（Fig 4B/C）
- GSE190609 加上皮分数协变量（~organ+patient+epi）：PM_vs_PT 51% 保留、PM_vs_LM 73% 保留；PM_core 34/49（69%）校正后仍显著；LM_core 57/57（100%）保留，且 **GSE50760 独立复现 57/57（100%）** → EPCAM 协变量管不住器官驻留细胞污染。**限定：GSE50760 参与 LM meta、非完全 held-out；EPCAM 代理不含肝细胞成分，故保留与 admixture 相容而非反证。**
- **定量归因（预审修订新增）**：56/57 LM_core 基因为已知肝功能基因；LM_core 基因级效应量与肝特异性相关（Spearman ρ = 0.635，*P* = 1.1×10⁻⁷）；12 基因肝标志物评分归因 → **quantitatively consistent with substantial hepatocyte admixture**。
- **PM_core 不对称对照（预审修订新增）**：10/49 注释为脂肪细胞/脂质处理、16/49 为成纤维/ECM 身份、0/49 间皮标志物；脂肪细胞评分归因约 97% 的 PM_core 信号、但归属混合生态位成分 → PM_core 不反映单一驻留细胞类型污染。
- 肝细胞标志物 LM vs PT log2FC：*ALB*：GSE190609 10.98 / GSE50760 6.11 / GSE41258 0.19；*APOB* 8.24/4.95/0.69；*APOA1* 9.03/7.78/0.57；*HP* 8.88/7.98/2.85；*SERPINA1* 2.60/4.08/0.99。（GSE41258 芯片倍数低是平台压缩效应，不要解释为"无污染"）
- **措辞锁定**：肝信号一律写 "quantitatively consistent with substantial hepatocyte admixture"（旧稿 "largely explained by hepatocyte contamination" 废弃）。

## 5. Results 证据分配（main-text discipline）

- **核心发现（主文）**：PM_core/LM_core meta 基因集（P1）、CMS4 GLMM 阳性对照（P1）、APC 缺失 + SMAD4 nominal（P2，样本级+患者级+logistic）、LOCO + 排列基线阴性（P3）、污染定量归因（P4）。
- **必要支撑（主文一句+SI）**：EMT sanity、方向一致率（280/282、145/146）、异质性（PM median I² = 0% / LM 96.9%）、固定 vs 随机效应（β 相关 0.988/0.766）、严格版 CMS（*P* = 0.060）、duplicateCorrelation（Jaccard = 0.904）、ssGSEA 保序性验证与 GSEApy 基准。
- **限定/敏感性（主文克制呈现）**：同子集 Cox（*P* = 0.181 → 0.102，功效损失框架）、LOCO 主行 0.736 nominal 显著的限定、model 2（n = 302）功效不足、左截断声明、肺 HR 抽样选择说明、卵巢排除敏感性。
- **SI**：CNA 177 例探索性、GSE41258/41568 参考队列细节、芯片 vs RNA-seq 分层 meta 细节、全部逐队列 DE 表。
- 禁止把 LOCO 阴性藏进 SI——这是卖点之一（有意义的阴性结果）。

## 6. 攻防平衡写作规则（03 规范 3.2.0，强制执行）

- 亮点前置：Abstract/Intro 末段/Results 开头先给 P1+P2 最强证据链。
- 不利结果：完整报告但降叙事中心性。**负荷校正 Cox 推荐句式（2026-10-04 修订，旧 "partly carried by metastatic burden" 句式废弃）**："peritoneal involvement was already non-significant in the subset model without the burden term (HR = 1.38, 95% CI 0.86–2.22, *P* = 0.181) and remained so after adding the number of metastatic sites (HR = 1.49, 95% CI 0.92–2.39, *P* = 0.102); the point estimate did not decrease — the loss of significance reflects reduced power in the smaller subset rather than metastatic burden explaining the prognostic signal." 禁写 "failed validation / performed poorly / lack of robustness"。
- 禁止措辞：robust clinical prediction model、novel signature、markedly improved、superior、first ever、comprehensive、mechanism confirmed；**以及预审修订新增禁用：「APC–SMAD4 (genomic) axis」、「SMAD4-enriched」、「partly carried by (metastatic burden)」、"largely explained by hepatocyte contamination"**。
- 推荐措辞：predefined gene set、transcriptomic context、was associated with、direction preserved、supports ... as a candidate requiring further validation；**新增锁定：「APC-depleted genomic background」、「nominal *SMAD4* enrichment that did not survive multiple-testing correction」、「positive control」（CMS4）、「not reliably supported — rather than disproven — at the current sample sizes」（P3）、「quantitatively consistent with substantial hepatocyte admixture」（P4）**。
- 结论强度校准：单队列相关→"associated with"；多队列方向一致→"consistently associated with"；无前瞻性验证→不写 clinically implementable；无实验→不写因果。

## 7. Discussion 五段结构（03 规范 3.2.0 §7）——预审修订后实为七段

1. 总结主要发现（P1+P2 最强链；CMS4 标注阳性对照；SMAD4 标注 nominal 不过校正）。
2. 生物学合理性与文献定位：腹膜间皮/脂肪微环境（[25,26] + FABP4 耦合 [27]）vs CMS4 间质型（[5,6]）互证；**基因组背景——serrated/WNT-ligand（*RNF43*–*RSPO*）框架仅作 literature context 保留，明确本数据不支持（*BRAF* 不富集、*RNF43* 趋势 ns 且卵巢驱动、*SMAD4* 属 TGF-β 与配体依赖 WNT 无机制链接、*APC* 缺失指向经典 WNT 输出减弱）**；同子集 Cox 功效损失框架；与 Golas 2025（CNA 层 [8]）和 Yu 2025（scRNA 层 [28]）互补；肝污染方法学信息（56/57 肝功能基因、ρ = 0.635；[9,10,11]）。
3. **相对既往工作的增量（新增段）**：CMS4 = 阳性对照非新发现；外部签名系统重叠分析；PM_core = 跨队列 meta 汇聚核心，首个三重汇聚腹膜核心集。
4. 优势：预定义生物学假设（stromal peritoneal niche + APC-depleted genomic background）、同患者多器官配对设计消除批次混淆 + 患者感知混合建模、跨平台 meta 方向一致 + 低腹膜异质性、独立大队列互证、诚实阴性结果（LOCO + 排列基线）完整报告。
5. 限制（克制）：回顾性公共数据；GSE190609 无纯肝对照患者；腹膜 meta 仅 k = 2（REML 异质性不精确、仅描述性）；肝转移污染估计依赖转录标志物；BRAF 跨队列不一致（选择差异）；MSI/TMB/sidedness 覆盖不全（model 2 功效不足、MSI-high 2 例不可评估）；OS 左截断；肺 HR = 0.60 为抽样选择；CNA 功效不足；LOCO 结论依赖敏感子集 + 排列基线；无蛋白/空间/实验验证。
6. **下一步（新增段）**：无即时临床可行性；前瞻多中心队列 + 匹配原发–转移采样（含腹膜）+ 空间/单细胞 + 预注册分析计划（检验原发灶 *APC* 缺失 → 后续腹膜复发；sidedness 与 MSI 预设分层）。
7. 结论（不过界）：生态位印记 + APC-depleted 基因组背景作为候选框架值得前瞻与实验验证；原发灶转录组预测转移部位 not reliably supported；肝转移 bulk 转录组研究须常规控制肝细胞污染。

## 8. Terminology Ledger（锁定形式，全稿统一）

- 样本类型缩写：PT（primary tumor）、PM（peritoneal metastasis）、LM（liver metastasis）、LungM（lung metastasis）——首次出现给全称。
- 基因集：peritoneal-tropism gene set (PM_core, 49 genes)、liver-tropism gene set (LM_core, 57 genes)。
- 方法：paired limma with patient blocking；random-effects meta-analysis (REML, metafor)；ssGSEA (rank-based)；leave-one-cohort-out (LOCO) sensitivity analysis；CMS（consensus molecular subtypes, CMSclassifier random-forest）；**patient-aware generalized linear mixed model（GLMM；binomial family，per-patient random intercept (1|patient)，lme4 1.1.30）**；**size- and detection-rate–matched random gene sets（1,000 permutations；empirical *P*）**。
- 队列：MSK-MET（Memorial Sloan Kettering–Metastatic Events and Tropisms）。
- 基因斜体（*APC*、*SMAD4*、*BRAF*）、蛋白正体（APC、SMAD4）；统计：*P*、*n*、*I*² 斜体；HR/OR 带 95% CI；*P*<0.001 或三小数。
- 通路表述：serrated-type / WNT-ligand–dependent (*RNF43*–*RSPO*) context **仅作 literature context**——本数据不支持该框架，禁用作本结果的解释；描述文献时用 "consistent with" 而非 "demonstrates"。
- 基因组背景锁定表述：「APC-depleted genomic background」；SMAD4 只能以 "nominal enrichment that did not survive multiple-testing correction (OR = 1.62, *P* = 0.048, FDR = 0.214)" 出现。
- APC 样本级主口径：64.5% vs 82.2%（腹腔内+卵巢合并）；**58.0%（47/81）仅允许出现在卵巢排除敏感性分析语境**。

## 9. 参考文献（2026-10-04 预审修订：28 条新编号，全部核验落位）

全部引用已落实为 **28 条编号引用**（按正文首次出现编号），唯一事实源：`manuscript/references_verified.tsv`（含 PMID/DOI/verdict/claim 对应）。正文与 sections 已同步替换，无残留 `[REF:`。

**旧（26 条）→ 新（28 条）编号映射**：3→4, 4→5, 5→6, 6→7, 7→8, 8 删除→由新 9/10 替代, 9→11, 10→12, 11 删除→由新 13 替代, 12→14, 13→15, 14→16, 15→17, 16→18, 17→19, 18→20, 19→21, 20→22, 21→23, 22→24, 23→25, 24→26, 25→27, 26→28。
**新增 4 条**：[3] Quenet 2021（PRODIGE 7，HIPEC 增量获益质疑）、[9] Ki 2007（肝转移签名须排除器官特异基因）、[10] Moosavi 2021（肝微环境浸润塑造转移灶转录组）、[13] Watanabe 2010（原发灶预测肝转移支持性证据）。
**删除 2 条（审计记录，勿再引用）**：旧 [8] Vasudevan 2020（仅脑膜瘤个案证据，证据强度不足，由新 [9]/[10] 替代）；旧 [11] Tie J（由新 [13] Watanabe 2010 替代）。

核验关键结论（承继旧蓝图，编号已更新）：
- **Saris2025 证伪**：真实存在的 Saris J 2025 论文为腹膜巨噬细胞研究，与 CMS4–腹膜关联无关；正文改引 [5] Lenos 2022 Nat Commun 13:4443（PMID 35927254）+ [6] Laoukili 2022 Br J Cancer 126:1824-33（PMID 35194192；即 GSE190609 源文献，一引两用）。
- **MSK-MET 更正**：正确出处为 Nguyen et al., **Cell** 2022;185(3):563-75.e11，PMID **35120664**（本节原记 Cancer Cell / PMID 35120658 有误，后者为无关论文）。
- **MSK-IMPACT-50K（Bandlamudi, Cancer Cell 2026）弃用**：摘要不含样本级器官标注，正文未引用。
- Guinney 2015 [4] 仅支持 CMS 框架本身，不支持腹膜关联；正文已拆句（CMS 框架引 [4]，腹膜-CMS4 引 [5,6]）。
- 方法学引用全部落实：limma [22]、metafor [23]、ssGSEA（"as described by Barbie et al."）[24]、CMSclassifier [4]、cBioPortal [20,21]。
- 数据集源文献：GSE190609→[6]、GSE50760→[15]、GSE225182→[16]、GSE41258→[17]、GSE41568→[18]。
- 讨论引用：腹膜/网膜生态位 [25,26]、FABP4 脂肪细胞耦合 [27]、Golas CNA [8]、Yu scRNA [28]、原发灶预测支持 [12,13] / 质疑 [14]。
- 年代结构：2023+ 文献 6/28（Lund-Andersen 2024、Gelli 2023、Ha 2024、Golas 2025、Dai 2025、Yu 2025），其余为方法学与数据集经典（符合 03 规范对奠基性方法学的豁免）。

## 10. 交付结构

- 各 section 草稿：`manuscript/sections/00_title_abstract.md` … `05_figure_legends.md`
- 汇总稿（主 agent 拼装）：`manuscript/main_manuscript.md`（2026-10-04 预审修订后为主稿唯一权威版本，sections 与本蓝图已同步）
- 补充声明段（主 agent 写）：Data availability（GEO accession + MSK-MET cBioPortal）、Code availability（scripts/ 列表；repository URL 已落实：https://github.com/wontaoxmu-pixel/crc-organotropism）、Acknowledgments（数据库致谢，格式：The authors gratefully acknowledge ... for making their data publicly available.）
