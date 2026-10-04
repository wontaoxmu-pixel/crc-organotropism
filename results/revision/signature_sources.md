# 已发表签名基因清单来源（A8 / RC-M1/S2）

提取日期 2026-10-04。四篇目标文献的签名均成功获得完整基因列表。

| 签名 | 文献 | 列表来源 | 获取途径 | 基因数 |
|---|---|---|---|---|
| Gelli2023_PM / Gelli2023_LM | Gelli et al., Cancers 2023;15:4418 (PMC10648258) | Supplementary Table S1（61 个 limma DEG，logFC>0=腹膜侧 40 / <0=肝侧 21） | Europe PMC supplementaryFiles → cancers-2534275-supplementary.pdf，pypdf 解析 | 40 / 21 |
| Laoukili2022_PM366 / PT138 | Laoukili et al., Br J Cancer 2022;126:1824-33 (PMC9174226，GSE190609 源文献) | Supplementary Table S2（PM 高于配对 PT 366 个 / PT 高于 PM 138 个，p<0.01 FDR） | Springer static-content MOESM2_ESM.docx 解析 | 366 / 138 |
| Lenos2022_PMhigh_up | Lenos et al., Nat Commun 2022;13:4443 (PMC9352687) | Supplementary Data 5（in vivo PM (CCLE) PM_high vs PM_low，上调 855 个；该文无单一离散分型签名，此为唯一公开 PM 相关基因列表） | Springer static-content MOESM7_ESM.xlsx | 855 |
| Yu2025_PM_malignant / Yu2025_LM_malignant | Yu et al., Cancer Lett 2025;629:217880 (PMID 40553883) | 补充 mmc5.xlsx：恶性细胞亚群 FindAllMarkers，PM 6 群（MMP7/LMO7/IGFBP7/TM4SF1/CA-like/REG4）并集 805 个，LM 6 群（RPL22L1/RPhigh/MTND1P23/CXCL14/MYC/GPRC5A）并集 1060 个 | Elsevier els-cdn 1-s2.0-S0304383525004483-mmc5.xlsx | 805 / 1060 |

注意：
- Gelli 61 基因中含多个 lncRNA/新转录本（LINC/AC 编号），与本研究 symbol 宇宙可匹配者有限（PM 侧 14/40）。
- Laoukili 旧式克隆号逗号小数（如 AC010969,1）已转为点号；符号匹配大小写不敏感。
- 超几何检验背景宇宙：主分析用本 meta 实际可检出基因（PM 宇宙 14,085；LM 宇宙 9,333），敏感性分析用 N=20,000 蛋白编码基因近似值。
