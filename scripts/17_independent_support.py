"""17_independent_support.py
Does DDR1 (and MKNK2 for comparison) show autism-vs-control differences in data NOT used to nominate it?
 (1) expression: GSE269105 (cerebellum RNA-seq, 7 ASD / 5 control; "ADD" sample B6994 dropped) and
                 GSE38322 (cerebellum 14/12, occipital 4/6; lumi-normalised Illumina arrays)
 (2) methylation (450K): GSE53162 (PFC / temporal / cerebellum), GSE53924 (BA10 frontal / BA24 cingulate), GSE278285 (cerebellum)
Caveats are printed with the results: GSE38322 shares 7 donors with GSE28521; every test is small-n; probe p-values are not corrected
across genome (only within gene), so these are descriptive checks, not discoveries.
Run:  ~/miniforge3/envs/docking/bin/python scripts/17_independent_support.py SCRATCH_DIR
  where SCRATCH_DIR contains GSE269105_fpkm.txt.gz, GPL10558.annot.gz, GPL13534.txt (downloaded earlier).
"""
import sys, os, gzip, io, subprocess, numpy as np, pandas as pd
from scipy import stats
S = sys.argv[1]; GENES = {"DDR1": "780", "MKNK2": "2872"}      # Entrez gene IDs (GSE269105 rows are keyed by Entrez ID)
os.makedirs("results", exist_ok=True)

def tt(a, b):                                   # autism minus control; Welch t + Mann-Whitney (missing values dropped)
    a, b = np.asarray(a, float), np.asarray(b, float); a, b = a[~np.isnan(a)], b[~np.isnan(b)]
    t = stats.ttest_ind(a, b, equal_var=False)
    return a.mean() - b.mean(), t.statistic, t.pvalue, stats.mannwhitneyu(a, b).pvalue

out = []
# ---------- (1a) GSE269105 ----------
hdr = pd.read_csv(f"{S}/GSE269105_fpkm.txt.gz", sep="\t", index_col=0); hdr.index = hdr.index.astype(str)
tpm = hdr[[c for c in hdr.columns if c.startswith("tpm_")]]; tpm.columns = [c[4:] for c in tpm.columns]
status = dict(zip(["CBell_B4925","CBell_B5000","CBell_B5144","CBell_B5173","CBell_B5251","CBell_B5718","CBell_B5873","CBell_B6337","CBell_B6401","CBell_B6677","CBell_B6994","CBell_B7333","CBell_B7369"],
                  ["ASD","ASD","ASD","ASD","unaffected","unaffected","unaffected","ASD","ASD","ASD","ADD","unaffected","unaffected"]))   # from the GEO sample table (checked 2026-10-19)
asd = [c for c in tpm.columns if status[c] == "ASD"]; ctl = [c for c in tpm.columns if status[c] == "unaffected"]
L = np.log2(tpm + 1); keep = (tpm >= 1).sum(axis=1) >= 5
allt = pd.Series({g: stats.ttest_ind(L.loc[g, asd], L.loc[g, ctl], equal_var=False).statistic for g in L.index[keep]}).dropna()
for g, eid in GENES.items():
    d, t, p, pm = tt(L.loc[eid, asd], L.loc[eid, ctl])
    out.append(("expression", "GSE269105", "cerebellum", g, f"{len(asd)} v {len(ctl)}", round(d, 3), round(p, 4), round(pm, 4),
                f"t={t:.2f}; rank of t among {len(allt)} expressed genes: {int((allt >= t).sum())+1} (top {100*((allt >= t).mean()):.0f}% = up)"))

# ---------- (1b) GSE38322 ----------
gz = gzip.open(f"{S}/GPL10558.annot.gz", "rt").read().split("\n"); i0 = next(i for i, l in enumerate(gz) if l.startswith("ID\t"))
a = pd.read_csv(io.StringIO("\n".join(gz[i0:])), sep="\t", low_memory=False, usecols=["ID", "Gene symbol"])
probes = {g: set(a.loc[a["Gene symbol"] == g, "ID"]) for g in GENES}
lines = gzip.open("data_raw/GSE38322/GSE38322_series_matrix.txt.gz", "rt").read().split("\n")
b0 = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_begin"))
mx = pd.read_csv(io.StringIO("\n".join(lines[b0 + 1:])), sep="\t", index_col=0).dropna(how="all"); mx = mx[~mx.index.str.startswith("!")]
tab = pd.read_csv("results/GSE38322_sample_table.csv").set_index("geo_accession")
for reg in ["cerebellum", "occipital"]:
    cols = [c for c in mx.columns if tab.loc[c, "region"] == reg]; A = [c for c in cols if tab.loc[c, "group"] == "autism"]; C = [c for c in cols if tab.loc[c, "group"] == "control"]
    ts = mx[cols].apply(lambda r: stats.ttest_ind(r[A], r[C], equal_var=False).statistic, axis=1).dropna()
    for g in GENES:
        pr = [p for p in probes[g] if p in mx.index]
        if not pr: out.append(("expression", "GSE38322", reg, g, f"{len(A)} v {len(C)}", np.nan, np.nan, np.nan, "no probe passed detection filter")); continue
        v = mx.loc[pr, cols].mean(axis=0); d, t, p, pm = tt(v[A], v[C])
        out.append(("expression", "GSE38322", reg, g, f"{len(A)} v {len(C)}", round(d, 3), round(p, 4), round(pm, 4), f"t={t:.2f}; {len(pr)} probe(s); rank of probe-t among {len(ts)}: {int((ts >= ts.loc[pr].mean()).sum())+1}"))

# ---------- (2) methylation ----------
pl = []
with open(f"{S}/GPL13534.txt") as fh:
    for ln in fh:
        if ln.startswith("ID\tName"): cols = ln.rstrip("\n").split("\t"); break
    ix = {c: i for i, c in enumerate(cols)}
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        if len(f) < len(cols) or ln.startswith("!") or ln.startswith("^"): continue
        names = f[ix["UCSC_RefGene_Name"]].split(";")
        for g in GENES:
            if g in names:
                grp = f[ix["UCSC_RefGene_Group"]].split(";"); k = names.index(g)
                pl.append((f[0], g, f[ix["CHR"]], f[ix["MAPINFO"]], grp[k] if k < len(grp) else "", f[ix["Relation_to_UCSC_CpG_Island"]]))
pm_ = pd.DataFrame(pl, columns=["probe", "gene", "chr", "pos", "region_in_gene", "cpg_island_rel"]).drop_duplicates("probe").set_index("probe")
pm_["promoter"] = pm_.region_in_gene.isin(["TSS1500", "TSS200", "5'UTR", "1stExon"])
print("450K probes per gene:", pm_.groupby("gene").size().to_dict(), "| promoter-class probes:", pm_.groupby("gene").promoter.sum().to_dict())
pm_.to_csv("results/17_methylation_probe_map.csv")

def fetch_matrix(gse):
    f = f"{S}/{gse}_series_matrix.txt.gz"
    if not os.path.exists(f):
        d = gse[:-3] + "nnn"; subprocess.run(["curl", "-s", "-L", "-o", f, f"https://ftp.ncbi.nlm.nih.gov/geo/series/{d}/{gse}/matrix/{gse}_series_matrix.txt.gz"], check=True)
    return f
def read_series(gse):
    meta, rows, ids = {}, [], set(pm_.index); inb = False
    for ln in gzip.open(fetch_matrix(gse), "rt"):
        if ln.startswith("!series_matrix_table_begin"): inb = True; continue
        if ln.startswith("!series_matrix_table_end"): break
        if inb:
            k = ln[:ln.index("\t")].strip('"')
            if k in ids or k == "ID_REF": rows.append(ln.rstrip("\n").split("\t"))
        elif ln.startswith("!Sample_geo_accession"): meta["gsm"] = [x.strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
        elif ln.startswith("!Sample_characteristics_ch1") or ln.startswith("!Sample_source_name_ch1"):
            meta.setdefault("char", []).append([x.strip('"') for x in ln.rstrip("\n").split("\t")[1:]])
    m = pd.DataFrame([[x.strip('"') for x in r] for r in rows[1:]], columns=[c.strip('"') for c in rows[0]]).set_index("ID_REF")
    m = m.apply(pd.to_numeric, errors="coerce")                      # blanks / "null" -> NaN
    return meta, m
def sample_labels(meta, gse):
    g = meta["gsm"]; txt = [" | ".join(t) for t in zip(*meta["char"])]
    grp = ["autism" if "autism" in t.lower() and "control" not in t.lower().replace("autism", "") else ("control" if ("control" in t.lower() or "normal" in t.lower()) else "?") for t in txt]
    reg = []
    for t in txt:
        tl = t.lower()
        reg.append("PFC/frontal" if ("ba10" in tl or "prefrontal" in tl or "frontal" in tl) else "cingulate" if ("ba24" in tl or "cingulate" in tl) else "temporal" if "temporal" in tl else "cerebellum" if "cerebell" in tl else "?")
    return pd.DataFrame({"group": grp, "region": reg}, index=g)

def read_278285():
    """GSE278285: betas are in a supplementary file whose columns are the sample TITLES; group from source_name, age/sex from characteristics."""
    hd = {}
    for ln in gzip.open(f"{S}/GSE278285_series_matrix.txt.gz", "rt"):
        if ln.startswith("!Sample_title"): hd["title"] = [x.strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
        elif ln.startswith("!Sample_source_name_ch1"): hd["src"] = [x.strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
        elif ln.startswith("!Sample_characteristics_ch1"):
            v = [x.strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
            if v[0].lower().startswith("age"): hd["age"] = [float(x.split(":")[1]) for x in v]
            elif v[0].lower().startswith("sex"): hd["sex"] = [x.split(":")[1].strip() for x in v]
    lab = pd.DataFrame({"group": ["autism" if "autism" in x.lower() else "control" for x in hd["src"]], "region": "cerebellum",
                        "age": hd["age"], "female": [1.0 if s_ == "F" else 0.0 for s_ in hd["sex"]]}, index=hd["title"])
    rows, ids = [], set(pm_.index)
    for ln in gzip.open(f"{S}/GSE278285_beta.txt.gz", "rt"):
        k = ln[:ln.index("\t")]
        if k in ids or k == "ID_REF": rows.append(ln.rstrip("\n").split("\t"))
    m = pd.DataFrame(rows[1:], columns=rows[0]).set_index("ID_REF").apply(pd.to_numeric, errors="coerce")
    return m, lab

def adj_test(y, lab, A, C):
    """beta ~ autism + age + female (OLS); returns group coefficient (dBeta) and its p-value."""
    idx = [i for i in A.tolist() + C.tolist() if not np.isnan(y[i])]
    X = np.column_stack([np.ones(len(idx)), [1.0 if i in set(A) else 0.0 for i in idx], lab.loc[idx, "age"].values, lab.loc[idx, "female"].values])
    yy = y[idx].values.astype(float); b, *_ = np.linalg.lstsq(X, yy, rcond=None); r = yy - X @ b
    s2 = r @ r / (len(yy) - X.shape[1]); se = np.sqrt(s2 * np.linalg.inv(X.T @ X)[1, 1]); t = b[1] / se
    return b[1], 2 * stats.t.sf(abs(t), len(yy) - X.shape[1])

for gse in ["GSE53162", "GSE53924", "GSE278285"]:
    if gse == "GSE278285": m, lab = read_278285()
    else: meta, m = read_series(gse); lab = sample_labels(meta, gse); lab = lab.loc[[c for c in m.columns if c in lab.index]]
    print(f"\n{gse}: {m.shape[0]} probes x {m.shape[1]} samples | groups {lab.group.value_counts().to_dict()} | regions {lab.region.value_counts().to_dict()}")
    for reg in sorted(lab.region.unique()):
        s = lab[lab.region == reg]; A = s.index[s.group == "autism"]; C = s.index[s.group == "control"]
        if len(A) < 3 or len(C) < 3: continue
        for g in GENES:
            pr = [p for p in pm_.index[pm_.gene == g] if p in m.index]
            res = []
            for p in pr:
                if gse == "GSE278285": d, pv = adj_test(m.loc[p], lab, A, C)          # adjusted for age + sex
                else: d, t, pv, pmw = tt(m.loc[p, A], m.loc[p, C])
                res.append((p, d, pv, pm_.loc[p, "promoter"]))
            r = pd.DataFrame(res, columns=["probe", "dBeta", "p", "promoter"])
            if r.empty: continue
            r["fdr"] = stats.false_discovery_control(r.p) if len(r) > 1 else r.p
            prm = r[r.promoter]; bod = r[~r.promoter]
            for lab_, sub in (("promoter", prm), ("body/other", bod)):
                if len(sub): out.append(("methylation_split", gse, reg, g, f"{len(A)} v {len(C)}", round(sub.dBeta.mean(), 4), f"{sub.p.min():.2e}", f"{(sub.p<0.05).sum()}/{len(sub)} probes p<0.05", lab_))
            out.append(("methylation", gse, reg, g, f"{len(A)} v {len(C)}", round(r.dBeta.mean(), 4), round(r.p.min(), 4), round(r.fdr.min(), 3),
                        f"{len(r)} probes; promoter probes {len(prm)}; mean promoter dBeta {prm.dBeta.mean():.4f}; probes with raw p<0.05: {int((r.p<0.05).sum())}; max |dBeta| {r.dBeta.abs().max():.3f}"))

R = pd.DataFrame(out, columns=["layer", "dataset", "region", "gene", "n_asd_v_ctl", "effect(log2 diff or mean dBeta)", "min_or_only_p", "MWU_p_or_min_FDR", "notes"])
R.to_csv("results/17_independent_support.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 130)
print("\n" + R.to_string(index=False))
