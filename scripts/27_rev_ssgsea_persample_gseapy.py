# 修订 A4 附：GSEApy 第二基准（独立实现 ssGSEA：经验CDF秩次 + alpha 加权）
# 输入：ssgsea/revision_pool_matrix_for_gseapy.tsv.gz（合并 PT 矩阵，scripts/27 导出）
# 输出：ssgsea/revision_persample_benchmark_gseapy.tsv
import sys
import numpy as np
import pandas as pd
import gseapy as gp

np.random.seed(123)
proc = "data/processed/crc_organotropism"

pm = pd.read_csv(f"{proc}/meta/peritoneal_tropism_genes.tsv", sep="\t")["gene"].tolist()
lm = pd.read_csv(f"{proc}/meta/liver_tropism_genes.tsv", sep="\t")["gene"].tolist()
gsets = {"PM_core": pm, "LM_core": lm}

mat = pd.read_csv(f"{proc}/ssgsea/revision_pool_matrix_for_gseapy.tsv.gz", sep="\t", index_col=0)
labels = pd.read_csv(f"{proc}/ssgsea/revision_pool_labels.tsv", sep="\t")
print("matrix:", mat.shape, file=sys.stderr)

def run_ssgsea(df):
    ss = gp.ssgsea(data=df, gene_sets=gsets, outdir=None,
                   min_size=1, max_size=10000, permutation_num=0,
                   sample_norm_method="rank", threads=4, seed=123, verbose=False)
    res = ss.res2d
    # 每行一个 (term, sample)：列含 Term / ES / NES / Name(sample)
    es = res.pivot_table(index="Name", columns="Term", values="ES")
    return es

def auc(case_scores, ctrl_scores):
    case_scores = np.asarray(case_scores, float)
    ctrl_scores = np.asarray(ctrl_scores, float)
    n1, n2 = len(case_scores), len(ctrl_scores)
    allv = np.concatenate([case_scores, ctrl_scores])
    r = pd.Series(allv).rank().values[:n1]  # 与 R rank() 一致（ties 取平均）
    return (r.sum() - n1 * (n1 + 1) / 2) / (n1 * n2)

def wilcox_p(case_scores, ctrl_scores):
    from scipy.stats import mannwhitneyu
    return mannwhitneyu(case_scores, ctrl_scores, alternative="two-sided").pvalue

# ---------- 合并池打分 ----------
es_pool = run_ssgsea(mat)
lab_map = labels.set_index("sample")["label"]

# 与 R 两版分数的 Spearman 相关
r_ps = pd.read_csv(f"{proc}/ssgsea/revision_persample_scores.tsv.gz", sep="\t", index_col=0)
r_m3 = pd.read_csv(f"{proc}/ssgsea/pooled_PT_ssgsea.tsv.gz", sep="\t", index_col=0)
rows = []
for gs in ["PM_core", "LM_core"]:
    v_gp = es_pool[gs]
    v_ps = r_ps.loc[gs]
    v_m3 = r_m3.loc[gs]
    common = v_gp.index.intersection(v_ps.index)
    rho_ps = pd.Series(v_gp[common]).corr(pd.Series(v_ps[common]), method="spearman")
    rho_m3 = pd.Series(v_gp[common]).corr(pd.Series(v_m3[common]), method="spearman")
    rows.append(dict(analysis="score_correlation", tag=f"合并PT池: {gs} GSEApy vs R per-sample版",
                     n1=len(common), spearman_rho_scores=rho_ps))
    rows.append(dict(analysis="score_correlation", tag=f"合并PT池: {gs} GSEApy vs R M3归一化版",
                     n1=len(common), spearman_rho_scores=rho_m3))

cfgs = [
    ("PT_PM", "PT_LM", "PM_core", "主检验1: PM_core | 腹膜患者PT vs 肝转移患者PT"),
    ("PT_LM", "PT_PM", "LM_core", "主检验2: LM_core | 肝转移患者PT vs 腹膜患者PT"),
    ("PT_PM", "PT_unselected", "PM_core", "参考: PM_core | 腹膜患者PT vs 未选择PT"),
    ("PT_LM", "PT_unselected", "LM_core", "参考: LM_core | 肝转移患者PT vs 未选择PT"),
    ("PT_LM", "PT_unselected", "PM_core", "阴性对照: PM_core | 肝转移患者PT vs 未选择PT"),
]
for case, ctrl, gs, tag in cfgs:
    s_case = es_pool.loc[lab_map[lab_map == case].index, gs]
    s_ctrl = es_pool.loc[lab_map[lab_map == ctrl].index, gs]
    a = auc(s_case, s_ctrl)
    p = wilcox_p(s_case, s_ctrl)
    rows.append(dict(analysis="auc_gseapy", tag=tag, n1=len(s_case), n2=len(s_ctrl),
                     auc_gseapy=a, p_gseapy=p))
    print(f"{tag} | GSEApy AUC={a:.3f} (P={p:.3g})", file=sys.stderr)

# ---------- GSE190609 内部配置（子矩阵单独打分） ----------
m190 = pd.read_csv(f"{proc}/GSE190609_meta.tsv", sep="\t")
pt = m190[m190["organ"] == "primary tumor"].copy()
haslm = (m190[m190["organ"] != "primary tumor"].groupby("patient")["organ"]
         .apply(lambda s: (s == "liver metastasis").any()).rename("hasLM"))
pt = pt.merge(haslm, on="patient")
cols = [f"GSE190609|{t}" for t in pt["title"]]
es_190 = run_ssgsea(mat[cols])
s_case = es_190.loc[[f"GSE190609|{t}" for t in pt.loc[pt["hasLM"], "title"]], "LM_core"]
s_ctrl = es_190.loc[[f"GSE190609|{t}" for t in pt.loc[~pt["hasLM"], "title"]], "LM_core"]
a = auc(s_case, s_ctrl)
p = wilcox_p(s_case, s_ctrl)
rows.append(dict(analysis="auc_gseapy",
                 tag="附: GSE190609 PT 内部 LM_core PM+LM vs PM-only（探索性）",
                 n1=len(s_case), n2=len(s_ctrl), auc_gseapy=a, p_gseapy=p))
print(f"内部配置 | GSEApy AUC={a:.3f} (P={p:.3g})", file=sys.stderr)

out = pd.DataFrame(rows)
out.to_csv(f"{proc}/ssgsea/revision_persample_benchmark_gseapy.tsv", sep="\t", index=False)
print(out.to_string(), file=sys.stderr)
print("done", file=sys.stderr)
