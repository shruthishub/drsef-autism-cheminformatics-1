"""22b_GSE298243_tests.py
Independent expression check in GSE298243 (BA9 snRNA, per-donor pseudobulk from 22a, NeuN+ ('Pos') and NeuN- ('Neg') fractions).
Gene t = autism coefficient from OLS: logCPM ~ autism + age + log(n_cells). Then, as for scripts 20/21:
 (1) are the 208 candidates enriched for t>0 vs all other genes?  (2) DDR1/MKNK2 and the 37 druggable candidates;  (3) myelin/oligodendrocyte sets (label permutation).
Caveats: all male; 9 donor IDs overlap Velmeshev (not removed); the Pos/Neg meaning is assumed; Neg pseudobulks are small (median ~50-110 nuclei).
Run: ~/miniforge3/envs/docking/bin/python scripts/22b_GSE298243_tests.py
"""
import numpy as np, pandas as pd
from scipy import stats
rng = np.random.default_rng(42); NPERM = 3000
pb = pd.read_csv("data_clean/GSE298243_pb_counts.csv", index_col=0); m = pd.read_csv("data_clean/GSE298243_pb_meta.csv")
M21 = pd.read_csv("results/21_candidate_evidence_matrix.csv").set_index("gene"); cand = list(M21.index)
GO = pd.read_csv("data_clean/go_bp_cc_sets.csv")
names = ["GOCC_MYELIN_SHEATH", "GOBP_OLIGODENDROCYTE_DIFFERENTIATION", "GOBP_OLIGODENDROCYTE_DEVELOPMENT", "GOBP_ENSHEATHMENT_OF_NEURONS", "GOBP_REGULATION_OF_MYELINATION",
         "GOBP_OXIDATIVE_PHOSPHORYLATION", "GOBP_INFLAMMATORY_RESPONSE", "GOBP_SYNAPTIC_SIGNALING"]
SETS = {n: set(GO.loc[GO.gs_name == n, "gene_symbol"]) for n in names}; SETS["MARKERS_oligo_curated"] = set("MBP PLP1 MOG MAG MOBP CNP OLIG1 OLIG2 SOX10 MYRF CLDN11 ST18 GJB1 ERMN TF".split())

def gene_t(Y, g, C):
    Q, _ = np.linalg.qr(C); g = g.reshape(len(g), -1); Ry = Y - (Y @ Q) @ Q.T; rg = g - Q @ (Q.T @ g)
    gg = (rg ** 2).sum(0); num = Ry @ rg; rss = (Ry ** 2).sum(1)[:, None] - num ** 2 / gg; df = Y.shape[1] - C.shape[1] - 1
    return (num / gg) / np.sqrt(np.maximum(rss, 1e-12) / df / gg)

allT, setrows, lines = {}, [], []
for fr in ("Pos", "Neg"):
    s = m[(m.fraction == fr) & (m.n_cells >= 20)].dropna(subset=["diagnosis", "age"]).reset_index(drop=True)
    X = pb[s.group_id]; X = X[(X >= 5).sum(axis=1) >= 10]; L = np.log2(X / X.sum() * 1e6 + 1)
    g = (s.diagnosis == "ASD").astype(float).values; C = np.column_stack([np.ones(len(s)), s.age.astype(float), np.log(s.n_cells)])
    A = L.values; t0 = gene_t(A, g, C)[:, 0]; T = pd.Series(t0, index=L.index); df = len(s) - C.shape[1] - 1; allT[fr] = T
    nA, nC = int(g.sum()), int((1 - g).sum())
    lines.append(f"\n===== {fr} fraction: {nA} ASD v {nC} CTL donors, {len(T)} genes tested; median nuclei per pseudobulk {s.n_cells.median():.0f}; global median t {T.median():.2f}, sd {T.std():.2f}")
    c = T.reindex(cand).dropna(); o = T.drop(c.index)
    lines.append(f"candidate list vs other genes: n={len(c)} median t {c.median():.2f} vs {o.median():.2f}; MWU p={stats.mannwhitneyu(c, o).pvalue:.3g}; frac t>0 {np.mean(c>0):.2f} vs {np.mean(o>0):.2f}")
    pas = c[(c > 0) & (2 * stats.t.sf(abs(c), df) < 0.05)]
    lines.append(f"candidates with t>0 and p<0.05: {len(pas)} (chance ~{0.025*len(c):.1f}): {', '.join(pas.index)}")
    dr = M21.index[M21.druggable]; cd = T.reindex(dr).dropna()
    lines.append(f"druggable candidates present: {len(cd)}; with t>0 & p<0.05: {', '.join(cd.index[(cd>0)&(2*stats.t.sf(abs(cd),df)<0.05)]) or 'none'}")
    for gname in ("DDR1", "MKNK2", "CSNK1E", "NR2F1", "PBX1", "NFIA"):
        if gname in T.index: lines.append(f"  {gname}: t={T[gname]:+.2f} p={2*stats.t.sf(abs(T[gname]),df):.3g} rank {int((T>=T[gname]).sum())}/{len(T)}")
        else: lines.append(f"  {gname}: not expressed enough")
    idx = {k: np.where(L.index.isin(v))[0] for k, v in SETS.items()}; perm = np.zeros((NPERM, len(idx)))
    for st in range(0, NPERM, 250):
        Gp = np.column_stack([rng.permutation(g) for _ in range(250)]); Tp = gene_t(A, Gp, C); perm[st:st + 250] = np.column_stack([Tp[i].mean(0) for i in idx.values()])
    for j, (k, i) in enumerate(idx.items()):
        if len(i) >= 5: obs = t0[i].mean(); setrows.append((fr, k, len(i), round(obs, 3), (1 + (np.abs(perm[:, j]) >= abs(obs)).sum()) / (1 + NPERM)))
S = pd.DataFrame(setrows, columns=["fraction", "gene_set", "n_genes", "mean_t(+ = up in ASD)", "p_label_permutation"]); S.to_csv("results/22_GSE298243_sets.csv", index=False)
pd.DataFrame(allT).to_csv("results/22_GSE298243_gene_t.csv")
print("\n".join(lines)); print("\n", S.to_string(index=False))
