# Metastatic organotropism in colorectal cancer

Analysis code and display items for the manuscript:

> **Metastatic organotropism in colorectal cancer: a peritoneal stromal–adipocytic imprint, an APC-depleted genomic background, and the limits of primary-tumor prediction**
>
> Tao Wang, Hao Yi, Peng Qi

## Overview

This repository contains all analysis scripts, the manuscript, and all main and supplementary display items for a study combining:

- a patient-matched, multi-organ transcriptomic meta-analysis of five public GEO cohorts (GSE190609, GSE50760, GSE225182, GSE41258, GSE41568), and
- the MSK-MET clinical-genomic cohort (3,548 colorectal cancer samples, accessed via cBioPortal).

Key analyses: paired limma differential expression with patient blocking; random-effects meta-analysis (metafor, REML); CMS classification (CMSclassifier) with patient-aware GLMM; single-sample gene-set scoring with leave-one-cohort-out (LOCO) validation and 1,000 size- and detection-matched random gene-set permutation baselines; MSK-MET sample- and patient-level mutation associations, multivariable logistic regression, and Cox survival models; epithelial-fraction purity adjustment and organ-resident-cell admixture attribution; external-signature overlap (hypergeometric tests).

## Repository layout

| Path | Contents |
|---|---|
| `scripts/` | All analysis and figure scripts (numbered in pipeline order) |
| `manuscript/` | Main manuscript (`main_manuscript.md`), per-section files, blueprint, verified reference list, pre-review report |
| `main_figures/` | Figures 1–4 (600-dpi PNG + vector PDF) with collision/alignment audit records |
| `supplementary/figures/` | Supplementary Figures S1–S5 (600-dpi PNG + vector PDF) with audit records |
| `supplementary/tables/` | Supplementary Tables S1–S10 (CSV) with README |
| `results/` | Revision-stage analysis outputs (permutation baselines, LOCO, GLMM, contamination, overlaps) |
| `mskmet/` | MSK-MET revision outputs (deduplication, ovary-excluded sensitivity, logistic models) |

Raw expression data are not redistributed: all cohorts are publicly available from GEO, and MSK-MET from cBioPortal (see manuscript Data availability).

## Reproduction

Environment: R 4.2.0 (Bioconductor 3.15) with data.table, readxl, limma, metafor, survival, jsonlite, CMSclassifier, org.Hs.eg.db, AnnotationDbi, hgu133a.db, hgu133plus2.db, lme4 1.1.30, ggplot2, patchwork, scales, ggrepel; Python benchmark uses GSEApy 1.3.1. The random seed is fixed (`set.seed(123)`) in all stochastic steps.

Run scripts in numeric order (`scripts/10_*.R` … `scripts/35_*.R`); figure scripts are `20–23` (main figures) and `35` (supplementary figures); `33–34` build the supplementary tables.

## Citation

If you use this code, please cite the manuscript (citation to be updated upon publication).

## License

Code is released under the MIT License unless otherwise noted.
