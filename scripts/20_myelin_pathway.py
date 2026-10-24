"""20_myelin_pathway.py
Option 2: are myelin / oligodendrocyte genes changed in autism brain (cerebellum first, cortex for contrast)?
Per cohort and region: gene-wise t for 'autism' from OLS with that cohort's covariates (Frisch-Waugh, so exact), then for each gene set
  (a) SELF-CONTAINED test: statistic = mean t over set genes; null = permuting the autism/control labels (3000x), which keeps the correlation between genes;
  (b) COMPETITIVE check: where the set's mean t falls among 2000 random gene sets of the same size drawn from the genes tested
      (guards against a global shift, e.g. GSE38322 where many genes move together).
Positive t = higher in autism. Small n everywhere; contrast sets (astrocyte, inflammation, OXPHOS, synapse) show what a 'normal' result looks like.
Cohorts: cerebellum = GSE28521 (21), GSE38322 (26; shares 7 donors with GSE28521), GSE269105 (12); cortex = GSE28521 frontal / temporal, GSE102741 DLPFC,
GSE64018 STG (shares donors with GSE28521), Velmeshev snRNA oligodendrocyte-only pseudobulk (PFC, ACC; expression INSIDE oligodendrocytes, not composition).
Run: ~/miniforge3/envs/docking/bin/python scripts/20_myelin_pathway.py
"""
import gzip, io, numpy as np, pandas as pd
from scipy import stats
rng = np.random.default_rng(42); NPERM, NRAND = 3000, 2000
GO = pd.read_csv("data_clean/go_bp_cc_sets.csv")
sym_of_entrez = GO.drop_duplicates("ncbi_gene").set_index("ncbi_gene").gene_symbol
SETS = {n: set(GO.loc[GO.gs_name == n, "gene_symbol"]) for n in [
    "GOCC_MYELIN_SHEATH", "GOBP_OLIGODENDROCYTE_DIFFERENTIATION", "GOBP_OLIGODENDROCYTE_DEVELOPMENT", "GOBP_ENSHEATHMENT_OF_NEURONS",
    "GOBP_REGULATION_OF_MYELINATION", "GOBP_MYELIN_ASSEMBLY",
    "GOBP_ASTROCYTE_DIFFERENTIATION", "GOBP_INFLAMMATORY_RESPONSE", "GOBP_OXIDATIVE_PHOSPHORYLATION", "GOBP_SYNAPTIC_SIGNALING"]}
SETS["MARKERS_oligo_curated(15)"] = set("MBP PLP1 MOG MAG MOBP CNP OLIG1 OLIG2 SOX10 MYRF CLDN11 ST18 GJB1 ERMN TF".split())   # hand-picked, not data-derived
CONTRAST = {"GOBP_ASTROCYTE_DIFFERENTIATION", "GOBP_INFLAMMATORY_RESPONSE", "GOBP_OXIDATIVE_PHOSPHORYLATION", "GOBP_SYNAPTIC_SIGNALING"}

def gene_t(Y, g, C):
    """t for the group coefficient of every row of Y (genes x n) given group g (n,) or (n,B) and covariate matrix C (n x k, with intercept)."""
    Q, _ = np.linalg.qr(C); g = g.reshape(len(g), -1); Ry = Y - (Y @ Q) @ Q.T; rg = g - Q @ (Q.T @ g)
    gg = (rg ** 2).sum(0); num = Ry @ rg; rss = (Ry ** 2).sum(1)[:, None] - num ** 2 / gg; df = Y.shape[1] - C.shape[1] - 1
    return (num / gg) / np.sqrt(np.maximum(rss, 1e-12) / df / gg)

def analyse(name, Y, group, cov=None):
    """Y: DataFrame genes x samples (log scale); group: 0/1 array (1 = autism); cov: DataFrame of covariates (rows = samples) or None."""
    Y = Y.loc[Y.var(axis=1) > 1e-8]; g = np.asarray(group, float)
    C = np.ones((len(g), 1)) if cov is None else np.column_stack([np.ones(len(g)), cov.values.astype(float)])
    A = Y.values.astype(float); t0 = gene_t(A, g, C)[:, 0]; T = pd.Series(t0, index=Y.index)
    idx = {k: np.where(Y.index.isin(v))[0] for k, v in SETS.items()}
    perm = np.zeros((NPERM, len(idx)))
    for s in range(0, NPERM, 250):
        G = np.column_stack([rng.permutation(g) for _ in range(250)]); Tp = gene_t(A, G, C)
        perm[s:s + 250] = np.column_stack([Tp[i].mean(0) for i in idx.values()])
    out = []
    for j, (k, i) in enumerate(idx.items()):
        if len(i) < 5: continue
        obs = t0[i].mean(); pperm = (1 + (np.abs(perm[:, j]) >= abs(obs)).sum()) / (1 + NPERM)
        rnd = np.array([t0[rng.choice(len(t0), len(i), replace=False)].mean() for _ in range(NRAND)])
        pcomp = (1 + (np.abs(rnd - rnd.mean()) >= abs(obs - rnd.mean())).sum()) / (1 + NRAND)
        out.append((name, k, "contrast" if k in CONTRAST else "target", int((group == 1).sum()), int((group == 0).sum()), len(i), round(obs, 3), round(pperm, 4), round(pcomp, 4)))
    return out, T

def sym_collapse(df):                                            # several rows per symbol -> mean
    return df.groupby(level=0).mean()

res, tstats = [], {}
def run(name, Y, group, cov=None):
    o, T = analyse(name, sym_collapse(Y), group, cov); res.extend(o); tstats[name] = T; print("done", name, Y.shape)

# ---------- cerebellum ----------
m = pd.read_csv("data_clean/GSE28521_cerebellum_meta.csv"); x = pd.read_csv("data_clean/GSE28521_cerebellum_expr.csv", index_col=0)[m["sample"]]
run("cerebellum | GSE28521", x, (m.group == "autism").astype(int).values)
for reg, lab in (("frontal cortex", "frontal"), ("temporal cortex", "temporal")):
    mc = pd.read_csv("data_clean/GSE28521_cortex_meta.csv"); mc = mc[mc.region == reg]
    xc = pd.read_csv("data_clean/GSE28521_cortex_expr.csv", index_col=0)[mc["sample"]]
    run(f"{lab} cortex | GSE28521", xc, (mc.group == "autism").astype(int).values)

ann = gzip.open("data_raw/geo_extra/GPL10558.annot.gz", "rt").read().split("\n"); i0 = next(i for i, l in enumerate(ann) if l.startswith("ID\t"))
a = pd.read_csv(io.StringIO("\n".join(ann[i0:])), sep="\t", low_memory=False, usecols=["ID", "Gene symbol"]).dropna().set_index("ID")["Gene symbol"]
lines = gzip.open("data_raw/GSE38322/GSE38322_series_matrix.txt.gz", "rt").read().split("\n"); b0 = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_begin"))
mx = pd.read_csv(io.StringIO("\n".join(lines[b0 + 1:])), sep="\t", index_col=0).dropna(how="all"); mx = mx[~mx.index.str.startswith("!")]
tab = pd.read_csv("results/GSE38322_sample_table.csv").set_index("geo_accession")
mx = mx.loc[mx.index.isin(a.index)]; mx.index = a.loc[mx.index].str.split(";").str[0].values       # first symbol per probe
for reg in ("cerebellum", "occipital"):
    cols = [c for c in mx.columns if tab.loc[c, "region"] == reg]
    run(f"{reg} | GSE38322", mx[cols], (tab.loc[cols, "group"] == "autism").astype(int).values, tab.loc[cols, ["age (years):ch1"]].astype(float))

fp = pd.read_csv("data_raw/geo_extra/GSE269105_fpkm.txt.gz", sep="\t", index_col=0); fp.index = fp.index.astype(str)
tpm = fp[[c for c in fp.columns if c.startswith("tpm_")]]; tpm.columns = [c[4:] for c in tpm.columns]
st = dict(zip(["CBell_B4925","CBell_B5000","CBell_B5144","CBell_B5173","CBell_B5251","CBell_B5718","CBell_B5873","CBell_B6337","CBell_B6401","CBell_B6677","CBell_B6994","CBell_B7333","CBell_B7369"],
              ["ASD","ASD","ASD","ASD","unaffected","unaffected","unaffected","ASD","ASD","ASD","ADD","unaffected","unaffected"]))
tpm = tpm[[c for c in tpm.columns if st[c] != "ADD"]]; tpm = tpm[(tpm >= 1).sum(axis=1) >= 5]; L = np.log2(tpm + 1)
L.index = [sym_of_entrez.get(int(i), np.nan) if i.isdigit() else i for i in L.index]; L = L[~pd.isna(L.index)]      # rows are a mix of Entrez IDs and symbols
run("cerebellum | GSE269105", L, np.array([1 if st[c] == "ASD" else 0 for c in L.columns]))

# ---------- cortex (other cohorts) ----------
m = pd.read_csv("data_clean/GSE102741_meta.csv"); x = pd.read_csv("data_clean/GSE102741_expr.csv", index_col=0)[m["sample"]]
cov = pd.DataFrame({"age": m.age, "female": (m.sex == "Female").astype(float), "rin": m.rin, "black": (m.race == "AA").astype(float)})
run("DLPFC | GSE102741", x, (m.group == "autism").astype(int).values, cov)
m = pd.read_csv("data_clean/GSE64018_meta.csv"); c = pd.read_csv("data_clean/GSE64018_counts.csv", index_col=0)[m["sample"]]
c = c[(c >= 10).sum(axis=1) >= 6]; lc = np.log2(c / c.sum() * 1e6 + 1)
run("STG | GSE64018", lc, (m.group == "autism").astype(int).values, pd.DataFrame({"age": m.age, "male": (m.sex == "M").astype(float), "rin": m.rin}))

# ---------- Velmeshev oligodendrocyte-only pseudobulk ----------
vm = pd.read_csv("data_clean/velmeshev_oligo_meta.csv"); vx = pd.read_csv("data_clean/velmeshev_oligo_logcpm.csv", index_col=0)
for reg in ("PFC", "ACC"):
    s = vm[vm.region == reg]
    cov = pd.DataFrame({"age": s.age.values, "male": (s.sex == "M").astype(float).values, "log_cells": np.log(s.n_cells.values)})
    run(f"oligodendrocytes only, {reg} | Velmeshev", vx[s.group_id], (s.diagnosis == "ASD").astype(int).values, cov)

R = pd.DataFrame(res, columns=["cohort_region", "gene_set", "kind", "n_asd", "n_ctl", "n_genes_in_set", "mean_t(+ = up in autism)", "p_label_permutation", "p_vs_random_sets"])
R.to_csv("results/20_myelin_pathway.csv", index=False)
pd.DataFrame(tstats).to_csv("results/20_gene_t_stats_by_cohort.csv")
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500)
for k in SETS:
    print("\n==", k); print(R[R.gene_set == k].drop(columns=["gene_set", "kind"]).to_string(index=False))
