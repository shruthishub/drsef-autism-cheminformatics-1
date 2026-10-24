"""22c_GSE298243_excluding_velmeshev_donors.py
Repeats the candidate-list enrichment of 22b after removing the 9 GSE298243 donors whose numeric IDs also appear in Velmeshev 2019
(a discovery cohort): 5144 5242 5294 5403 5408 5565 5577 5841 5978. (ID match only; confirm in the papers that these are the same people.)
Run: ~/miniforge3/envs/docking/bin/python scripts/22c_GSE298243_excluding_velmeshev_donors.py
"""
import pandas as pd, numpy as np
from scipy import stats
m = pd.read_csv("data_clean/GSE298243_pb_meta.csv"); v = pd.read_csv("data_clean/velmeshev_pb_meta.csv")
ov = sorted(set(m.donor.astype(str)) & set(v.individual.astype(str))); print("overlap donors:", ov, len(ov))
pb = pd.read_csv("data_clean/GSE298243_pb_counts.csv", index_col=0); cand = list(pd.read_csv("results/21_candidate_evidence_matrix.csv").gene)
def gene_t(Y, g, C):
    Q, _ = np.linalg.qr(C); Ry = Y - (Y @ Q) @ Q.T; rg = g - Q @ (Q.T @ g); gg = rg @ rg; num = Ry @ rg; rss = (Ry ** 2).sum(1) - num ** 2 / gg; df = Y.shape[1] - C.shape[1] - 1
    return (num / gg) / np.sqrt(rss / df / gg)
for fr in ("Pos", "Neg"):
    s = m[(m.fraction == fr) & (m.n_cells >= 20) & (~m.donor.astype(str).isin(ov))].dropna(subset=["diagnosis", "age"]).reset_index(drop=True)
    X = pb[s.group_id]; X = X[(X >= 5).sum(axis=1) >= 10]; L = np.log2(X / X.sum() * 1e6 + 1)
    g = (s.diagnosis == "ASD").astype(float).values; C = np.column_stack([np.ones(len(s)), s.age.astype(float), np.log(s.n_cells)])
    T = pd.Series(gene_t(L.values, g, C), index=L.index); c = T.reindex(cand).dropna(); o = T.drop(c.index)
    print(fr, f"excluding overlap: {int(g.sum())} v {int((1-g).sum())}; candidates median t {c.median():.2f} vs {o.median():.2f}, MWU p={stats.mannwhitneyu(c, o).pvalue:.3g}; frac>0 {np.mean(c>0):.2f} v {np.mean(o>0):.2f}")
