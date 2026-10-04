# 预审稿报告（nature-reviewer 模拟）

- 生成日期：2026-10-04
- 审稿件：`manuscript/main_manuscript.md`（含 26 条已核验参考文献）
- 方法：3 名互盲审稿人（各自独立子代理上下文，仅见同一不可变审稿包 + 各自侧重点简报）+ 1 份冻结后综合
- 侧重点简报：A = 技术可靠性/统计与生信严谨性；B = 临床意义/转化相关性；C = 原创性/广泛兴趣/可读性
- 互盲声明：三份报告由三个互不通信的独立上下文生成，综合前未向任何审稿人展示其他报告；重叠与分歧均为独立评审的真实证据，未做事后改写。

## Review setup

- **Input scope**：全文稿件（摘要、正文、方法、参考文献、图注）+ 证据蓝本（blueprint.md）+ 参考文献核验台账（references_verified.tsv）
- **Assessment boundary**：无图件实物（仅图注）；无代码/原始数据可复跑；参考文献真伪已另行核验，预审只评估引用-claim 匹配
- **Shared manuscript claim summary**：5 个 GEO 配对队列 meta + MSK-MET 3,548 例 → ①PM 间质/脂肪生态位印记（PM_core, 49 基因）+ CMS4 富集；②腹膜播散伴 APC 缺失/SMAD4 富集基因组背景；③原发灶转录组不能可靠预测转移部位（LOCO）；④肝转移 bulk 信号大部分由肝细胞污染解释
- **Visible evidence base**：正文全部数字可追溯至 blueprint；26 条文献均有 PMID/DOI 核验记录
- **Missing materials affecting confidence**：图件实物、分析代码、R/包版本、仓库 URL（正文仍残留 `[Evidence needed: ...]` 占位）

---

## Reviewer A 报告（侧重点：技术可靠性 / 统计与生信严谨性）

- **Overall assessment**：稿件把配对多器官转录组 meta 分析、CMS 分型、MSK-MET 基因组关联和 LOCO 阴性结果组织成一条论点链，写作克制、阴性结果前置，态度可取。但从统计与生信技术角度看，四条支柱中至少三条存在方法学缺口：(i) 原发灶预测分析中"转移部位标签与队列/平台几乎完全混杂"，使表观 AUC 与 LOCO AUC 都无法作生物学解释，且自定义 ssGSEA 引入跨样本归一化，破坏了 single-sample 性质；(ii) CMS4 富集与 MSK-MET 样本级突变比较均把同一患者的多个样本当作独立观测，存在伪重复；(iii) 上皮分数校正的逻辑不能支撑"肝细胞污染解释 LM 信号"这一结论；(iv) "APC–SMAD4 基因组轴"写进了标题与摘要，但 SMAD4 在任何一层均未通过多重校正（患者级 FDR = 0.214）。核心框架（PM_core 生态位印记 + APC 缺失）在两队列 REML 和患者级分析下有立足点，但多个定量主张在当前证据下不能成立。倾向 major revision。
- **Who would be interested**：结直肠癌转移生物学与器官趋向性研究者；做 bulk 转移灶转录组的方法学人群（肝细胞污染警示）；做原发灶预测模型的转化研究者（LOCO 阴性结果是证伪型约束）。
- **Major strengths**：①同患者多器官配对设计 + 患者阻断 limma；②阴性结果完整前置并声明依赖敏感性子集；③转录组发现与独立大队列患者级分析互证；④负荷校正 Cox、严格阈值 CMS、平台分层 meta 等敏感性分析均如实报告边界。

### Major Concerns（A）

**RA-M1**（Blocking: Yes；technical soundness）
- Claim pointer：Results §4 "The primary-tumor transcriptome does not reliably encode metastatic destination"（:93-97）及 Methods §Single-sample scoring（:55）
- Evidence pointer：GSE190609 全部 12 例均有腹膜转移，GSE50760 全部 18 例均为肝转移三联（:39, :73），标签与队列/平台近乎一一对应；Methods 写明自定义打分使用 cross-sample range normalization（:55）
- Concern：(a) 标签-队列混杂：表观 AUC = 0.933 本质上可由"区分 GSE50760 与其余队列"的批次结构驱动，LOCO 主行（AUC = 0.736, P = 0.003）同此；没有随机基因集置换基线或队列内交叉验证，"表观信号来自循环论证"与"LOCO 下失败"都无法归因于生物学。(b) 自定义 ssGSEA 的跨样本极差归一化使单样本分数依赖同批次其他样本构成，破坏 single-sample 性质，与 Barbie 2009 定义冲突；以 GSVA 版本兼容为由弃用已验证实现却无基准一致性对比。
- Why it matters：P3（原发灶不可预测）是标题三大主张之一，当前设计既不能证真也不能证伪生物学假设。
- Resolution test：①随机基因集置换（1,000 次，匹配基因数与检出率）的 AUC 零分布；②队列内验证或按队列分层报告 AUC；③与 GSVA 官方实现的基准对比（或改纯 per-sample 归一化重跑）；④每个 AUC 附 95% CI。

**RA-M2**（Blocking: No；technical soundness）
- Claim pointer：CMS4 富集检验（:75, OR = 3.13, P = 0.006）与 MSK-MET 样本级突变比较（:87）
- Evidence pointer：CMS 检验基于 113 样本但仅 12 例患者（:39）；MSK-MET 样本级基于 1,147 转移灶样本，未说明每患者单样本
- Concern：两处 Fisher 精确检验均存在伪重复（同一患者多样本当独立观测）；59 个 PM 样本来自 12 例患者，CMS 分型患者内高度相关，P = 0.006 反保守；6 个 LM 样本中 CMS4 仅 1 个，OR 估计极不稳定。
- Why it matters：CMS4 富集是 P1 的阳性对照支柱；样本级 APC 是摘要直接引用数字。
- Resolution test：CMS4 改患者级（或混合效应/聚类稳健检验）或明确声明探索性；MSK-MET 说明每患者样本数分布，以患者级为主报告。

**RA-M3**（Blocking: Yes；technical soundness）
- Claim pointer：Results §4 "largely explained by admixed hepatocyte transcripts"（:99）及摘要同句（:19）
- Evidence pointer：纯度校正仅用 EPCAM/KRT8/KRT19 三基因均值作上皮分数协变量（:63）；保留率 PM_core 34/49、LM_core 57/57（:99）
- Concern：非 sequitur 推理链。三基因上皮代理与肝细胞含量无已建立协变关系：若捕捉不到肝细胞比例，保留率结果对污染假说无信息量；若污染伴随低上皮分数，100% 保留反而与"污染解释"张力并存。真正支持污染的只有 5 个标志基因 log2FC，不能量化 57 个基因中多大比例来自肝细胞。对称地，PM_core 的脂肪/间质信号同样可能来自网膜组织混杂，双侧解读标准不一致。
- Why it matters：P4（污染警示）是标题与结论三大主张之一，证据形式不足以支撑 "largely explained by" 的强度。
- Resolution test：用去卷积方法（CIBERSORTx/xCell 或 GTEx 肝特异性基因集）定量估计器官驻留细胞比例，检验 LM_core 效应量与肝特异性程度相关；PM_core 用同一框架报告脂肪/间质贡献；若仅能做标志物论证，措辞降为 "consistent with substantial hepatocyte admixture"。

**RA-M4**（Blocking: Yes；technical soundness / claim calibration）
- Claim pointer：标题 "an APC–SMAD4 genomic axis"（:3）；摘要 "with SMAD4 mutations enriched"（:19）
- Evidence pointer：SMAD4 样本级 24.7% vs 17.0% ns（:87）；患者级 OR = 1.62, P = 0.048, FDR = 0.214（:89）
- Concern：SMAD4 在两个分析层级均未通过多重校正，标题与摘要却把它提升为与 APC 并列的"轴"，与自身数据矛盾。患者级比较未校正原发灶侧别、MSI 等已知混杂；BRAF 频率低不能替代直接分层分析。
- Why it matters：标题级主张与数据不匹配是审稿阶段必须解决的校准问题。
- Resolution test：标题摘要改为数据支撑的形式（如 "APC-depleted genomic background, with a nominal SMAD4 trend"）；补患者级侧别/MSI 分布与多变量 logistic 回归敏感性分析。

**RA-M5**（Blocking: No；technical soundness）
- Claim pointer："attenuated after adjustment ... suggesting the prognostic signal is partly carried by metastatic burden"（:91, :105）
- Evidence pointer：未校正模型 n = 1,136 时 HR = 1.42；负荷校正模型 n = 405 时 HR = 1.49
- Concern：HR 从 1.42 变为 1.49，点估计并未衰减；显著性消失更可能反映样本量收缩（1,136→405）的功效损失，且两模型拟合在不同子集，"attenuated" 的因果解读不成立。"unadjusted" 模型实为含肝/肺受累、BRAF、年龄的多变量模型，命名不当。
- Why it matters：该句是 Discussion 对腹膜预后意义的定性结论，方向性误读会误导读者。
- Resolution test：n = 405 子集上同时拟合含与不含 met_count 的模型，报告 HR 变化；若点估计基本不变，改写为"在更小的负荷标注子集中不再显著，无法区分负荷校正还是功效损失"。

**RA-M6**（Blocking: No；technical soundness）
- Claim pointer：PM 程序随机效应 meta（282 个显著基因，:81）与 Methods（:47）
- Evidence pointer：全文未报告 τ²、I² 或 Q；PM meta 仅 k = 2 队列（其一仅 7 对）；LM meta 4 队列含两个非配对芯片参考队列
- Concern：k = 2 时 REML 对 τ² 估计极不稳定，随机效应权重近乎任意；旗舰产物 PM_core 建立在两队列 REML 之上且无任何异质性统计量。方向一致率 280/282 提供部分安慰但不能替代正式异质性报告。
- Why it matters：PM_core 构成直接决定 P1 及全部下游分析。
- Resolution test：报告各 meta 的 τ²/I²/Q 并声明 k = 2 局限；补 leave-one-cohort-out 的 meta 敏感性；考虑附固定效应模型对照。

### Minor Comments（A）

- **RA-m1**（technical soundness；Results §3 与 Fig 3C 图注）："unadjusted Cox" 实为多变量模型 → 改称 "multivariable model without metastatic-burden term"，图注同步。
- **RA-m2**（technical soundness；:59, :87）：卵巢转移并入腹膜 compartment 有争议 → 补剔除卵巢样本的敏感性分析或论证合并理由。
- **RA-m3**（technical soundness；:99）：LM_core 定义与保留率评估用同一队列同一对比，轻度循环，100% 保留率有乐观偏倚 → 在 GSE50760/芯片队列重复，或文中明示。
- **RA-m4**（引用匹配；:31, :105 引 [8]）：单病例脑膜瘤甲基化研究支撑 CRC 肝转移污染论断偏弱 → 补直接文献或降级为方法学背景。
- **RA-m5**（可重复性；:67, :117）：`[Evidence needed: R version/package versions/repository URL]` 占位残留，"scripts are publicly available" 无法验证 → 投稿前填实。
- **RA-m6**（technical soundness；:95-97）：全部 AUC 无 95% CI，多配置并行检验未说明多重性 → 补 DeLong/bootstrap CI 并声明策略。
- **RA-m7**（technical soundness；:43）：LM 对比患者自由度极小且未报告 6 个 LM 样本的患者来源；多区域采样建模细节未说明 → 报告各对比实际患者数与样本数。

- **Technical failings（A）**：RA-M1、RA-M3、RA-M4 解决前核心论点不成立；另需 RA-M2 修正与 RA-m5 占位符填实。
- **五轴评价**：originality——框架组合有新意，组件均有在先文献；scientific importance——若成立对领域有约束性价值，当前可靠性不足；interdisciplinary readership——三个社群各取所需，达标；technical soundness——最弱一轴，需系统性返工；readability——清晰克制，但标题/摘要 SMAD4 表述误导非专业读者。
- **Recommendation posture**：major revision 倾向；若无法提供去卷积定量与批次置换分析（RA-M1、RA-M3），则 reject-and-resubmit。

---

## Reviewer B 报告（侧重点：临床意义 / 转化相关性）

- **Overall assessment**：稿件论证腹膜转移的间质/脂肪生态位印记与 APC 缺失基因组背景，并诚实报告两个阴性/警示结果。稿件与 blueprint 数字一致，写作克制，局限段坦诚。主要问题集中在临床转化层面：标题摘要的「APC–SMAD4 基因组轴」中 SMAD4 未通过多重校正；摘要对 LOCO 阴性结果选择性呈现；「锯齿状/WNT 配体依赖」框架与自身 BRAF、RNF43 数据方向相悖；生存分析时间原点与不朽时间偏倚未交代。总体倾向 major revision。
- **Who would be interested**：消化道肿瘤外科与内科医生、CRC 转移机制研究者、转录组方法学使用者、CMS/WNT 通路研究者。外科医生读者会注意到缺少可操作临床路径，兴趣在机制与阴性警示而非实践改变。
- **Major strengths**：①同患者多器官配对设计罕见且优于非配对比较；②阴性结果呈现在主文与摘要且写明敏感性子集限定；③肝细胞污染警示可立即采纳；④转录组层与基因组层 orthogonal 互证；⑤引用-claim 对应经逐条核验。

### Major Concerns（B）

**RB-M1**（Blocking: Yes；technical soundness / 临床可行动性）
- Claim pointer：标题（:3）、摘要 Results/Conclusions（:19, :21）、结论段（:111）
- Evidence pointer：SMAD4 样本级 ns（:87）；患者级 OR = 1.62, P = 0.048, FDR = 0.214（:89）；RNF43 P = 0.142
- Concern：双基因「轴」的两端中仅 APC 缺失有统计支持。对临床读者，「APC–SMAD4 轴」暗示两个可检测、潜在可指导分层的标志物，实际证据只支持 APC。
- Why it matters：标题摘要是临床读者主要接触面；以未显著信号命名轴会误导后续研究设计与引用。
- Resolution test：标题摘要降格为 APC 单信号（SMAD4 作方向性趋势并明确 FDR 阴性），或独立队列复现 SMAD4 富集（校正后显著）。

**RB-M2**（Blocking: No；阴性结果的临床含义）
- Claim pointer：摘要（:19）"did not reliably predict ... (AUC = 0.596 and 0.532, both nonsignificant)"
- Evidence pointer：正文（:97）与 limitations（:109）写明主 LOCO 配置仍显著（AUC = 0.736, P = 0.003），阴性依赖检出基因限制子集（0.596, P = 0.229）；blueprint 4.4 明确要求写明该限定
- Concern：正文遵守了限定，摘要没有。只读摘要的临床读者会把「证据不一、严格子集下失效」误读为「已被证伪」，对「是否值得继续做原发灶预测」的含义截然不同。
- Why it matters：阴性结论是三大卖点之一，不能从一个方向的过度声称换成另一个方向。
- Resolution test：摘要补入主 LOCO 行仍显著的事实与敏感性限定，措辞改为 "not reliably supported" 量级。

**RB-M3**（Blocking: No；scientific importance / 临床框定准确性）
- Claim pointer：摘要与 Discussion 的 "serrated-like, WNT-ligand–dependent context"（:33, :105）
- Evidence pointer：MSK-MET 腹膜病灶 BRAF 仅 7.4% 且不富集（FDR = 0.298，:87），与 GSE190609 的 32% 严重不一致；RNF43 ns；SMAD4 为 TGF-β 通路基因，与 RNF43–RSPO 配体依赖 WNT 无机制连接
- Concern：该框架由文献 [2,6] 外推，但本稿自身基因组数据（BRAF 不富集、RNF43 ns、SMAD4 ns）均不站在该框架一侧；唯一稳健的 APC 缺失恰恰指向经典 WNT 减弱而非配体依赖激活。Discussion 未正面讨论此张力。
- Why it matters：框架暗示特定临床亚群（锯齿状、BRAF、潜在 PORCN 抑制剂敏感），由不支持它的数据背书会把转化讨论引向错误方向。
- Resolution test：降格为文献背景假设，明确列出自身数据与该框架的不一致，或提供队列内锯齿状/MSI 分层直接检验。

**RB-M4**（Blocking: No；预后表述的临床准确性）
- Claim pointer：Results 小节标题句（:87）与 Cox 结果（:91）
- Evidence pointer：Methods（:59）含 "age at sequencing"，提示 OS 自测序日期起算；MSK-MET 为转移灶测序队列，存在不朽时间/超前时间偏倚风险，且无治疗信息（CRS+HIPEC、化疗线数）；肺受累 HR = 0.60（P = 0.002）的反直觉结果未获解释
- Concern：腹膜 HR = 1.42 作为小节核心句呈现，但关联在负荷校正后消失，且时间原点未定义。对以预后判断为日常工作的临床读者，这直接决定该数字可否被引用。
- Why it matters：采样时点差异可部分驱动部位-生存关联；不明确时间原点，HR 有效性无法评估。
- Resolution test：明确 OS 起点并讨论偏倚（或 landmark 敏感性分析）；解释或降格肺受累保护性 HR；小节标题句改以负荷校正结果为参照。

### Minor Comments（B）

- **RB-m1**（临床框定；:27）：CRS+HIPEC 适用性论断无直接文献支撑，且未提 PRODIGE 7 对 HIPEC 获益的质疑，框定滞后 → 补直接文献或删子句。
- **RB-m2**（引用匹配；:31 引 [8]）：单病例异癌种证据支撑器官特异性方法学论断力度不匹配 → 换 CRC/肝转移特异文献或改写为本研究观察。
- **RB-m3**（引用匹配；:31 引 [11]）：Tie 2011 为 KRAS 基因型-复发部位关联，与转录组预编程论证轴不同 → 改写明确证据层次或替换。
- **RB-m4**（technical soundness；:59 vs :39）：卵巢归类标准两个数据体系不一致 → 论证合并理由或补排除卵巢敏感性分析。
- **RB-m5**（originality；:75）：CMS4 关联已由 [4,5] 建立且 LM 对照 n = 6，实质为验证性阳性对照却以发现性语气呈现 → 明确 corroboration 定位并注明 n = 6 限制。
- **RB-m6**（转化落点；:27 vs :111）：动机（surveillance/治疗设计）与落点（"warrants prospective validation"）之间空档，外科读者 "so what" 缺位 → Discussion 补一小段明确当前无临床可行动性并给出一条具体前瞻性验证设计（队列构成、终点、分层变量）。
- **RB-m7**（可重复性；:55, :67）：自写 ssGSEA 无基准验证而 P3 完全建立在该打分器上；占位符残留 → 补与标准实现的一致性基准，填实占位符。

- **Technical failings（B）**：RB-M1、RB-M2、RB-M3、RB-M4 + RB-m7。
- **五轴评价**：originality 中等（增量在配对 meta + APC 缺失 + 诚实 LOCO，属渐进式）；scientific importance 中等偏上（取决于 RB-M1/M2 解决）；interdisciplinary readership 较好但转化路径缺位削弱临床获得感；technical soundness 骨架扎实但三处呈现超出统计支持、打分器未基准验证；readability 良好，摘要选择性呈现是主要扣分点。
- **Recommendation posture**：major revision 倾向；四条 Major 中三条不需新实验即可解决；RB-M1 若坚持保留轴命名需独立复现。若 SMAD4 轴与摘要选择性呈现维持现状，则 reject-and-resubmit。

---

## Reviewer C 报告（侧重点：原创性 / 广泛兴趣 / 可读性）

- **Overall assessment**：四层主张组合清晰，内部数字一致性良好，写作纪律值得肯定。主要问题：①相对 Lenos 2022、Laoukili 2022（即 GSE190609 源文献）、Yu 2025、Golas 2025 的增量从未显式划界，CMS4 富集实质是对同一队列已发表结论的重复推导；②标题摘要把未通过校正的 SMAD4 提为"轴"的一半支柱；③正文与图注存在多处图版错位。以当前形态更适合专科期刊（blueprint 自定的 CCR/BJC 定位现实）；投高影响综合期刊则原创性增量与广泛可读性均未达标。
- **Who would be interested**：CRC 转移生物学与转化研究者；计算肿瘤学群体（污染警示与"上皮协变量不足"这一反直觉发现）；原发灶预测/风险分层监测研究者。后两条方法学信息受众面其实比主线更宽，但叙事位置靠后削弱了广泛兴趣。
- **Major strengths**：①配对设计 + 跨平台 meta，方向一致率 280/282 与 145/146；②LOCO 阴性完整报告并写明限定；③对 BRAF 跨队列不一致、负荷校正失显著等不利结果处理克制透明；④文献台账可追溯。

### Major Concerns（C）

**RC-M1**（Blocking: Yes；originality / scientific importance）
- Claim pointer：引言（:29）与 Discussion（:105）；摘要把 CMS4 富集作为第一支柱（:19）
- Evidence pointer：[5] Laoukili 2022 即 GSE190609 源文献，标题已含 "Peritoneal metastases ... belong to CMS4"；[4] Lenos 2022 与 [26] Yu 2025 分别在 bulk 分型与单细胞层报道腹膜-间质趋向；稿件未提供 PM_core/LM_core 与任何已发表签名（Lenos 2022、Gelli 2023 [10]、Yu 2025）的重叠或一致性比较
- Concern：第一支柱在多大程度上是再确认而非新发现，读者无法判断；用与 Laoukili 2022 完全相同的队列重新得出 CMS4–腹膜关联并放进摘要主结果，使原创性增量模糊化。
- Why it matters：对综合期刊增量划界是生死问题；对专科期刊，缺失的签名重叠分析也让人无法定位贡献边界。
- Resolution test：三段式划界——(a) CMS4 改述为阳性对照/重复验证；(b) PM_core 与 Lenos 2022、Gelli 2023、Yu 2025 签名的重叠统计（超几何检验或方向一致性）；(c) 一句话说明独有增量（同患者配对 + 跨队列 meta + PM-vs-LM 交集定义 + LOCO 阴性 + 污染警示的组合）。

**RC-M2**（Blocking: Yes；readability / claim calibration）
- Claim pointer：标题（:3）、running title（:5）、摘要 Results（:19）与 Conclusions（:21）
- Evidence pointer：SMAD4 患者级 OR = 1.62, P = 0.048, FDR = 0.214；样本级 ns
- Concern：三处标题级位置把仅名义显著的信号与 FDR = 0.003 的 APC 缺失并列为"轴"，读者会得到两个等强度关联的错误印象。与稿件自定措辞纪律相悖。
- Why it matters：标题是稿件被引用检索时的永久声称；以未通过多重检验的结果命名轴可被统计审稿一票否掉。
- Resolution test：标题与 running title 移除或降级 SMAD4；摘要补明 "nominal enrichment that did not survive multiple-testing correction"；坚持保留则需独立验证。

**RC-M3**（Blocking: No；readability / 图注自明性）
- Claim pointer：:75 "(Fig. 1C)"；:81 "(Fig. 2B)"；:99 "(Fig. 2C, Fig. 4B)"；:79 "(Fig. 2A)"
- Evidence pointer：图注（:160, :162）显示 Fig. 1C 为 OR 森林图（47.5%/16.7% 属 Fig. 1B）；Fig. 2B/2C 均为逐基因森林图而正文用以指代"交集产生 PM_core/LM_core"与"上皮校正保留率"，图注中无对应面板；Fig. 2A 图注为 meta 火山图而正文引作单队列火山图
- Concern：至少四处正文–图注面板错位；无论图注滞后还是正文写错，图注自明性目前不成立。
- Why it matters：图注错位是核对证据链的直接障碍，提示图件版本管理失序。
- Resolution test：逐面板核对四张主图，统一正文引用与图注面板字母；补"正文引用句→面板"逐句核对记录。

**RC-M4**（Blocking: No，若统计审稿确认则升级；technical soundness 跨轴提示）
- Claim pointer：摘要与 :75 CMS4×腹膜 OR = 3.13, P = 0.006（Fisher，n = 113 样本）
- Evidence pointer：113 样本仅来自 12 名患者（59 个 PM 聚于 12 人），LM 对照仅 6 样本（1/6 CMS4）；MSK-MET 样本级检验未说明多灶去重
- Concern：Fisher 假设观测独立，伪重复风险；LM n = 6 使 OR 对单样本翻转极敏感；CMS4 头条关联无患者级稳健性分析。
- Why it matters：若患者级重分析后失稳，P1 的 CMS4 佐证需整体降调。
- Resolution test：患者级（代表样本或多数表决）CMS4×器官敏感性分析或患者聚类方法；MSK-MET 说明多灶处理规则。

### Minor Comments（C）

- **RC-m1**（readability；Fig 4A 图注 :166）：图注无 AUC 数值，"LOCO/sensitivity AUC" 合并了性质不同的两个配置 → 分行列出两组 AUC 与 P 值并注明推翻依赖受限子集。
- **RC-m2**（引用匹配；:31 引 [8]）：单病例脑膜瘤论文语境相距甚远 → 替换/补充直接文献或降为辅证。
- **RC-m3**（readability；摘要 Methods :17）：patient blocking / random-effects meta / ssGSEA / LOCO 术语无解释连发 → 摘要用功能性表述，缩写留正文。
- **RC-m4**（交付完整性；:67, :7-9, :117, :121-125）：占位符残留与 "scripts are publicly available" 并存 → 投稿前补齐；URL 未落实前弱化该句。
- **RC-m5**（一致性）：图内文字英式 tumour vs 正文美式 tumor → 投稿包阶段统一重导（已自记）。
- **RC-m6**（图注自明性；Fig 2A :162）："blue, liver-tropism genes detected in the PM meta-analysis" 易误读 → 改 "blue, LM_core genes overlaid for reference"。

- **Technical failings（C）**：RC-M2 与 RC-M1 解决前核心叙事不能按现标题成立；RC-M4 需患者级敏感性分析兜底。
- **五轴评价**：originality 中等偏低（真正新的是配对 meta 框架与两条方法学阴性信息）；scientific importance 中等（主线偏确认性，方法学纠偏有价值）；interdisciplinary readership 目前不达标（术语密度高、无概念图、方法学信息靠后）；technical soundness 内部一致且诚实但有三处需处理；readability 结构清晰但图注错位与摘要术语墙扣分。
- **Recommendation posture**：高影响综合期刊标准下 reject-and-resubmit 倾向；专科期刊（CCR/BJC）定位下 major revision——解决 RC-M1 至 RC-M4 后是扎实论文。

---

## Cross-review synthesis（post-review；未向审稿人展示）

### 共识强项（≥2 名审稿人独立指出）

1. 同患者多器官配对设计 + 患者阻断 limma，结构上优于同类研究（A、B、C）。
2. 阴性结果前置且写明敏感性子集限定，学术诚实度罕见（A、B、C）。
3. 转录组层与基因组层 orthogonal 互证（A、B）。
4. 肝细胞污染警示与原发灶预测 LOCO 约束具有即时方法学价值，受众面甚至可能宽于主线发现（A、B、C）。
5. 引用-claim 核验纪律良好（B、C）。

### 共识 Blocking concerns（≥2 名审稿人独立提出同一底层问题）

**S1. SMAD4「基因组轴」声称超出数据（RA-M4 + RB-M1 + RC-M2，三方共识，全部标 Blocking）**
标题、running title、摘要 Results 与 Conclusions 把 SMAD4（样本级 ns；患者级 P = 0.048, FDR = 0.214）与 APC 缺失（FDR = 0.003）并列。解决路径一致：标题摘要降格为 APC 单信号 + SMAD4 标注为名义趋势，或提供独立队列校正后复现。**这是修订的第一优先级，纯文字工作，不需新分析。**

**S2. 原创性增量未划界（RC-M1 + RB-m5）**
CMS4 富集是用 Laoukili 2022 [5] 同一队列重复推导其已发表结论，却以发现语气进入摘要第一支柱；PM_core 与任何已发表签名无重叠比较。解决：CMS4 改述阳性对照 + 签名重叠统计 + 一句话增量声明。

### 其他共识 Major concerns

**S3. 伪重复 / unit-of-analysis（RA-M2 + RC-M4）**
CMS4 Fisher（113 样本/12 患者）与 MSK-MET 样本级突变比较（1,147 样本）均把同一患者多病灶当独立观测。解决：患者级敏感性分析或聚类稳健方法；MSK-MET 说明多灶处理规则。

**S4. 摘要选择性呈现 LOCO 结果（RB-M2 + RC-m1）**
摘要只引两个非显著 AUC，略去主 LOCO 行 AUC = 0.736（P = 0.003）仍显著。解决：摘要补全限定，措辞 "not reliably supported"；Fig 4A 图注同步分行列出。

**S5. 生存分析表述链（RA-M5 + RB-M4）**
A 指出 "attenuated" 误读（HR 1.42→1.49 点估计未降，n 1,136→405 功效损失）；B 指出 OS 时间原点（自测序起算）与不朽时间偏倚未交代、肺保护性 HR 未解释。解决：同子集对照模型 + 时间原点声明 + 小节标题句改以负荷校正结果为参照。

**S6. 自定义 ssGSEA 未基准验证（RA-M1b + RB-m7）**
跨样本极差归一化破坏 single-sample 性质，P3 完全建立在该打分器上。解决：与 GSVA/GSEApy 基准对比或改 per-sample 归一化重跑。

### 单一审稿人提出的重要 Major（未形成共识但值得重视）

- **RA-M1a（标签-队列混杂）**：表观 AUC 与 LOCO AUC 都可能由批次结构驱动，缺随机基因集置换基线——这是三份报告中最锐利的一条，直接威胁 P3 的可解释性。
- **RA-M3（上皮协变量逻辑）**：三基因代理管不住肝细胞比例，"largely explained by" 需去卷积定量支撑或措辞降级。
- **RA-M6（k = 2 REML）**：全文无 τ²/I²/Q，旗舰产物 PM_core 的 meta 基础单薄。
- **RB-M3（锯齿/WNT 框架冲突）**：自身 BRAF/RNF43/SMAD4 数据均不支持该框架，需正面讨论而非一句带过。
- **RC-M3（图注错位 ≥4 处）**：投稿前必须逐面板核对。

### 审稿人分歧点

- **定位判断**：C 认为高影响综合期刊不达标、专科期刊（CCR/BJC）现实；A、B 未直接评期刊层级但均给 major revision 而非 reject，隐含认为问题可在专科期刊框架内解决。三方对「数据本身有价值、问题集中在声称校准与方法学归因」判断一致。
- **P3 阴性结果的价值**：A 视其为「当前设计既不能证真也不能证伪」（最严厉）；B 接受其临床含义但要求摘要完整呈现；C 视其为广泛兴趣点之一。
- **肝细胞污染结论**：A 要求去卷积定量否则措辞降级；B 认为现有标志物 + 校正实验已构成「可立即采纳的纠偏」。严厉程度不同，方向一致。

### Minor revision 清单（合并去重）

1. "unadjusted Cox" 改名（RA-m1）。
2. 卵巢并入腹膜 compartment 的敏感性分析或论证（RA-m2 + RB-m4，共识）。
3. LM_core 定义与保留率评估的队列循环性说明（RA-m3）。
4. 引文 [8] 支撑力弱（RA-m4 + RB-m2 + RC-m2，三方共识）：补直接文献或降级。
5. 填实全部 `[Evidence needed]` / `[To be completed]` 占位符（RA-m5 + RC-m4，共识）。
6. AUC 补 95% CI + 多重性声明（RA-m6）。
7. 报告 GSE190609 各对比实际患者数与 LM 样本来源（RA-m7）。
8. CRS+HIPEC 框定补文献或删子句（RB-m1）。
9. 引文 [11] 层次改写或替换（RB-m3）。
10. Discussion 补转化落点小段（RB-m6）。
11. 摘要 Methods 术语功能化改写（RC-m3）。
12. 图内 tumour/tumor 统一（RC-m5）。
13. Fig 2A 蓝色基因措辞改写（RC-m6）。

### Broad-interest / significance readout

主线生物学发现偏确认性（CMS4-腹膜、FABP4-脂肪细胞均有在先文献），真正的广泛价值在两条方法学信息：bulk 肝转移转录组的肝细胞污染警示、原发灶预测在 LOCO 下的失效。当前叙事把二者放在第三、四位支柱，若投专科期刊无需调整；若日后瞄准更高层级，应考虑把方法学警示前置。

### 投稿前必须解决的最重要问题（排序）

1. **S1（SMAD4 轴降级）**：标题/running title/摘要/结论四处，纯文字，当天可完成。
2. **S4（摘要 LOCO 完整呈现）**：纯文字，当天可完成。
3. **S2（增量划界三段式）**：需要一次签名重叠统计（超几何检验），约半天分析。
4. **S3（患者级 CMS4 + MSK-MET 多灶说明）**：需要重跑两个分析，约 1 天。
5. **RA-M1a + S6（置换基线 + ssGSEA 基准）**：需要新分析，约 1-2 天；若不做则 P3 措辞必须降为「当前数据无法区分生物学信号与批次结构」。
6. **RA-M3（去卷积定量或措辞降级）**：去卷积约 1 天；否则措辞降级为 "consistent with substantial hepatocyte admixture"。
7. **S5（生存分析同子集对照 + 时间原点声明）**：约半天。
8. **RC-M3（图注面板核对）**：对照 display_item_manifest.tsv 逐面板核对。

## Risk / unsupported claims

- 「APC–SMAD4 基因组轴」：SMAD4 端无校正后支持（三方共识 Blocking）。
- 「serrated-like, WNT-ligand–dependent context」：自身 BRAF/RNF43 数据不支持，SMAD4 属 TGF-β 通路与该框架无机制连接（RB-M3）。
- 「attenuated after burden adjustment → 预后信号部分由负荷承载」：点估计未衰减，子集不同，因果解读不成立（RA-M5）。
- 「largely explained by hepatocyte contamination」：三基因上皮代理无法定量污染比例（RA-M3）。
- 摘要 LOCO 表述：省略主配置仍显著的事实，构成选择性呈现（RB-M2）。
- 不可评估项：图件实物与正文引用对应（RC-M3 指出的错位只能核对图注推断）；代码可重复性（仓库 URL 缺失）。
