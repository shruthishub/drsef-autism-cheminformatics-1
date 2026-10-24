"""21_candidate_evidence_matrix.py
Option 1: run the DDR1 checks (independent expression, methylation, de novo mutations, literature count, druggability) on ALL 208 candidates in
results/12_candidates_full.csv, so DDR1 can be compared with every other candidate on the same footing.
Layers (flags declared BEFORE looking at results; every flag is 'supports a gene that is UP in autism being a real autism gene'):
  E1  GSE269105 cerebellum RNA-seq (7 v 5; independent donors, truly independent of discovery): t>0 and p<0.05
  E2  GSE38322 cerebellum (14 v 12; shares 7 donors with GSE28521; broad shift cohort so t is robust-z standardised across all genes): z>1.64
  M   promoter methylation (TSS1500/TSS200/5'UTR/1stExon probes averaged per gene): lower in autism (hypomethylation goes with higher expression)
      relative to each cohort's own genome-wide shift: GSE278285 cerebellum (30 v 45, adjusted age+sex) robust z<-1.64 AND lower in GSE53162 cerebellum (centred dBeta<0)
  D   de novo: >=1 autism coding LoF, or missense excess Poisson p<0.01 (denovo-db, length-only null)
Results are screening flags in small cohorts: expected number of chance passes per layer is printed so the observed counts can be judged.
Run: ~/miniforge3/envs/docking/bin/python scripts/21_candidate_evidence_matrix.py   (needs results/20_gene_t_stats_by_cohort.csv from script 20)
"""
import gzip, os, numpy as np, pandas as pd
from scipy import stats
G = "data_raw/geo_extra"
cand = pd.read_csv("results/12_candidates_full.csv"); genes = list(cand.gene)
T = pd.read_csv("results/20_gene_t_stats_by_cohort.csv", index_col=0)

# ---------- expression layers (gene t from script 20; positive = up in autism) ----------
def p_from_t(t, df): return 2 * stats.t.sf(abs(t), df)
e1 = T["cerebellum | GSE269105"].dropna(); e2 = T["cerebellum | GSE38322"].dropna(); e0 = T["cerebellum | GSE28521"].dropna()
rz = lambda s: (s - s.median()) / (1.4826 * (s - s.median()).abs().median())
M = pd.DataFrame(index=genes)
M["t_GSE269105"] = e1.reindex(genes); M["p_GSE269105"] = p_from_t(M.t_GSE269105, 10)
M["t_GSE38322"] = e2.reindex(genes); M["z_GSE38322_robust"] = rz(e2).reindex(genes)
M["t_GSE28521_cb"] = e0.reindex(genes)
M["E1"] = (M.t_GSE269105 > 0) & (M.p_GSE269105 < 0.05); M["E2"] = M.z_GSE38322_robust > 1.64
print(f"genes with expression data: GSE269105 {M.t_GSE269105.notna().sum()}/208, GSE38322 {M.t_GSE38322.notna().sum()}/208, GSE28521 cb {M.t_GSE28521_cb.notna().sum()}/208")
# Is the candidate LIST enriched for upregulation in independent data compared with all other genes?  (candidates were chosen as UP in cortex discovery cohorts)
enr = []
for nm, s in (("GSE269105 cerebellum", e1), ("GSE38322 cerebellum", rz(e2)), ("GSE28521 cerebellum (used for support in script 10)", e0)):
    c = s.reindex(genes).dropna(); o = s.drop(c.index, errors="ignore")
    enr.append((nm, len(c), round(c.median() - o.median(), 3), stats.mannwhitneyu(c, o).pvalue, round((c > 0).mean(), 3), round((o > 0).mean(), 3)))
enr = pd.DataFrame(enr, columns=["cohort", "n_candidates", "median_t_candidates_minus_other_genes", "MWU_p", "frac_candidates_t>0", "frac_other_t>0"])
print("\nCandidate list vs all other genes (cerebellum, independent of cortex discovery):\n", enr.to_string(index=False))

# ---------- methylation layer ----------
# Genome-wide promoter-mean methylation per gene is computed for ALL genes, so each cohort's candidate effects can be centred on that cohort's own
# genome-wide shift (GSE278285 shows a near-global +1% shift in autism, which is not gene-specific and must not be read as a DDR1 effect).
rows = []
with open(f"{G}/GPL13534.txt") as fh:
    for ln in fh:
        if ln.startswith("ID\tName"): cols = ln.rstrip("\n").split("\t"); break
    ix = {c: i for i, c in enumerate(cols)}
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        if len(f) < len(cols) or ln[0] in "!^": continue
        names = f[ix["UCSC_RefGene_Name"]].split(";"); grp = f[ix["UCSC_RefGene_Group"]].split(";")
        for k, nm in enumerate(names):
            if k < len(grp) and grp[k] in ("TSS1500", "TSS200", "5'UTR", "1stExon"): rows.append((f[0], nm))
pm = pd.DataFrame(rows, columns=["probe", "gene"]).drop_duplicates()
print(f"\npromoter probes mapped: {len(pm)} for {pm.gene.nunique()} genes; candidates covered: {len(set(genes) & set(pm.gene))}/208")

def gene_means(m):                                   # probes x samples -> genes x samples (mean of promoter probes)
    s = pm[pm.probe.isin(m.index)]
    return m.loc[s.probe].groupby(s.gene.values).mean()

def adj_ols(Y, g, cov):                              # all genes at once: beta ~ autism + covariates; returns dBeta, t
    X = np.column_stack([np.ones(len(g)), g] + ([cov] if cov is not None else []))
    Y = np.where(np.isnan(Y), np.nanmean(Y, axis=1, keepdims=True), Y)               # a few missing betas -> gene mean
    XtXi = np.linalg.inv(X.T @ X); B = Y @ X @ XtXi; R = Y - B @ X.T; df = len(g) - X.shape[1]
    se = np.sqrt((R ** 2).sum(1) / df * XtXi[1, 1]); return pd.DataFrame({"d": B[:, 1], "t": B[:, 1] / se, "p": 2 * stats.t.sf(abs(B[:, 1] / se), df)})

def read_series_all(path):
    meta, rows_ = {"char": []}, []; inb = False
    for ln in gzip.open(path, "rt"):
        if ln.startswith("!series_matrix_table_begin"): inb = True; continue
        if ln.startswith("!series_matrix_table_end"): break
        if inb: rows_.append([x.strip('"') for x in ln.rstrip("\n").split("\t")])
        elif ln.startswith("!Sample_geo_accession"): meta["gsm"] = [x.strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
        elif ln.startswith("!Sample_characteristics_ch1") or ln.startswith("!Sample_source_name_ch1"): meta["char"].append([x.strip('"') for x in ln.rstrip("\n").split("\t")[1:]])
    m = pd.DataFrame(rows_[1:], columns=rows_[0]).set_index("ID_REF").apply(pd.to_numeric, errors="coerce")
    txt = [" | ".join(t).lower() for t in zip(*meta["char"])]
    grp = ["autism" if "autism" in t and "control" not in t.replace("autism", "") else ("control" if ("control" in t or "normal" in t) else "?") for t in txt]
    reg = ["frontal" if ("ba10" in t or "prefrontal" in t or "frontal" in t) else "cerebellum" if "cerebell" in t else "other" for t in txt]
    return m, pd.DataFrame({"group": grp, "region": reg}, index=meta["gsm"])

meth = {}; ids = set(pm.probe)
hd = {}
for ln in gzip.open(f"{G}/GSE278285_series_matrix.txt.gz", "rt"):
    if ln.startswith("!Sample_title"): hd["title"] = [x.strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
    elif ln.startswith("!Sample_source_name_ch1"): hd["src"] = [x.strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
    elif ln.startswith("!Sample_characteristics_ch1"):
        v = [x.strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
        if v[0].lower().startswith("age"): hd["age"] = [float(x.split(":")[1]) for x in v]
        elif v[0].lower().startswith("sex"): hd["sex"] = [x.split(":")[1].strip() for x in v]
rws = []
for ln in gzip.open(f"{G}/GSE278285_beta.txt.gz", "rt"):
    k = ln[:ln.index("\t")]
    if k in ids or k == "ID_REF": rws.append(ln.rstrip("\n").split("\t"))
b = pd.DataFrame(rws[1:], columns=rws[0]).set_index("ID_REF").apply(pd.to_numeric, errors="coerce")[hd["title"]]
gm = gene_means(b); g_ = np.array([1.0 if "autism" in x.lower() else 0.0 for x in hd["src"]])
r = adj_ols(gm.values, g_, np.column_stack([hd["age"], [1.0 if s == "F" else 0.0 for s in hd["sex"]]])); r.index = gm.index
meth["GSE278285_cb"] = r; print(f"GSE278285 cerebellum {int(g_.sum())} v {int((1-g_).sum())}, adjusted age+sex, {len(gm)} genes")
for gse, regs in (("GSE53162", ("cerebellum", "frontal")), ("GSE53924", ("frontal",))):
    m, lab = read_series_all(f"{G}/{gse}_series_matrix.txt.gz"); gm = gene_means(m)
    for rg in regs:
        s = lab[(lab.region == rg) & lab.group.isin(["autism", "control"])]; s = s.loc[[i for i in s.index if i in gm.columns]]
        gg = (s.group == "autism").astype(float).values
        if gg.sum() < 3 or (1 - gg).sum() < 3: continue
        r = adj_ols(gm[s.index].values, gg, None); r.index = gm.index; meth[f"{gse}_{rg}"] = r; print(f"{gse} {rg}: {int(gg.sum())} v {int((1-gg).sum())}")
shift = []
for k, r in meth.items():
    zc = (r.t - r.t.median()) / (1.4826 * (r.t - r.t.median()).abs().median())
    shift.append((k, len(r), round(r.d.median(), 4), round((r.d > 0).mean(), 3), round((r.p < 0.05).mean(), 3)))
    M[f"dBeta_{k}"] = r.d.reindex(genes); M[f"dBeta_centred_{k}"] = (r.d - r.d.median()).reindex(genes); M[f"zrobust_{k}"] = zc.reindex(genes)
print("\nGenome-wide methylation shift per cohort (all genes' promoters; a median far from 0 / fraction >0 far from 0.5 / p<0.05 fraction far above 0.05 = global shift):")
print(pd.DataFrame(shift, columns=["cohort", "n_genes", "median_dBeta", "frac_dBeta>0", "frac_p<0.05"]).to_string(index=False))
# M flag: lower methylation than the cohort's own genome-wide centre, robust z<-1.64 in the largest cohort and the same direction (centred dBeta<0) in the second cerebellum cohort
M["M"] = (M["zrobust_GSE278285_cb"] < -1.64) & (M["dBeta_centred_GSE53162_cerebellum"] < 0)

# ---------- de novo layer ----------
dn = pd.read_csv("results/18_denovo_gene_counts.csv")
lof = dn[dn["class"] == "LoF"].set_index("gene"); mis = dn[dn["class"] == "missense"].set_index("gene")
M["denovo_LoF_autism"] = lof.autism_obs.reindex(genes); M["denovo_mis_autism"] = mis.autism_obs.reindex(genes); M["denovo_mis_p"] = mis.poisson_p_excess.reindex(genes)
M["D"] = (M.denovo_LoF_autism >= 1) | (M.denovo_mis_p < 0.01)

# ---------- attach candidate info ----------
keep = ["gene", "z_meta", "FDR_meta", "cereb_logFC", "cereb_P", "module", "known_sfari", "n_mech", "max_phase", "mech_drugs", "n_potent", "pubmed_asd", "pubmed_total"]
M = cand[keep].set_index("gene").join(M)
for c in ["E1", "E2", "M", "D"]: M[c] = M[c].fillna(False).astype(bool)
M["n_layers"] = M[["E1", "E2", "M", "D"]].sum(axis=1)
M["druggable"] = (M.n_mech.fillna(0) > 0) | (M.n_potent.fillna(0) >= 10)
M.reset_index().to_csv("results/21_candidate_evidence_matrix.csv", index=False)

print("\nDDR1 / MKNK2 methylation, raw vs centred:"); print(M.loc[["DDR1","MKNK2"],[c for c in M.columns if c.startswith(("dBeta","zrobust"))]].round(4).T.to_string())
print("\nLayer pass counts among 208 candidates (expected by chance in brackets; E1 ~ 208x0.025, E2 ~ 208x0.05, M ~ 208x0.05x0.5, D from script 18 ~ few):")
for c, exp in (("E1", 208 * 0.025), ("E2", 208 * 0.05), ("M", 208 * 0.05 * 0.5)): print(f"  {c}: {int(M[c].sum())}  [{exp:.1f}]")
print(f"  D: {int(M.D.sum())}   (11 genes with an autism LoF and 41 with a missense were counted in script 18)")
print("  number of layers passed:", M.n_layers.value_counts().sort_index().to_dict(), "| druggable candidates:", int(M.druggable.sum()))
pd.set_option("display.width", 260); pd.set_option("display.max_colwidth", 40)
show = ["z_meta", "t_GSE269105", "p_GSE269105", "z_GSE38322_robust", "dBeta_centred_GSE278285_cb", "zrobust_GSE278285_cb", "dBeta_centred_GSE53162_cerebellum", "denovo_LoF_autism", "denovo_mis_autism", "E1", "E2", "M", "D", "n_layers", "n_mech", "n_potent", "pubmed_asd"]
print("\nDRUGGABLE candidates (n_mech>0 or >=10 potent compounds), sorted by layers passed then discovery z:")
print(M[M.druggable].sort_values(["n_layers", "z_meta"], ascending=False)[show].round(3).to_string())
print("\nTop non-druggable candidates with >=2 layers:"); print(M[(~M.druggable) & (M.n_layers >= 2)].sort_values(["n_layers", "z_meta"], ascending=False)[show].round(3).to_string())
