"""22a_pseudobulk_GSE298243.py
GSE298243 (snmCT-seq, BA9, 26 ASD / 23 CTL donors, all male): sums single-nucleus RNA counts per donor x sorted fraction
(cell-name suffix 'Pos' = NeuN+ (assumed neurons), 'Neg' = NeuN- (assumed glia); the README does not define this, so it is an assumption)
and writes data_clean/GSE298243_pb_counts.csv + GSE298243_pb_meta.csv.  Cell-type labels are NOT in the GEO count matrix.
Run: ~/miniforge3/envs/docking/bin/python scripts/22a_pseudobulk_GSE298243.py
"""
import gzip, re, numpy as np, pandas as pd
ser = {}
for ln in gzip.open("data_raw/candidates/GSE298243_series_matrix.txt.gz", "rt"):
    if ln.startswith("!Sample_title"): ser["title"] = [x.strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
    elif ln.startswith("!Sample_characteristics_ch1"):
        v = [x.strip('"') for x in ln.rstrip("\n").split("\t")[1:]]; k = v[0].split(":")[0].strip()
        if k in ("donor", "diagnosis", "age", "Sex") or k.lower() in ("rin", "pmi"): ser[k] = [x.split(":", 1)[1].strip() if ":" in x else "" for x in v]
meta = pd.DataFrame({k: v for k, v in ser.items() if k != "title"}); print(meta.head(3).to_string()); print(meta.diagnosis.value_counts().to_dict())
donors = set(meta.donor)
f = "data_raw/GSE298243/GSE298243_single_cell_countMatrix.tsv.gz"
cols = gzip.open(f, "rt").readline().rstrip("\n").split("\t")[1:]
def cell_group(c):
    parts = c.split("_"); d = next((p for p in parts if p in donors), None)
    fr = "Pos" if "Pos" in c or "Pos2" in c else ("Neg" if "Neg" in c else None)
    return None if d is None or fr is None else f"{d}|{fr}"
grp = pd.Series([cell_group(c) for c in cols]); print("cells:", len(cols), "| unmatched:", int(grp.isna().sum()), "| groups:", grp.nunique())
ncell = grp.value_counts()
codes, uniq = pd.factorize(grp); ok = codes >= 0
import scipy.sparse as sp
S = sp.csr_matrix((np.ones(ok.sum()), (np.where(ok)[0], codes[ok])), shape=(len(cols), len(uniq)))
out, names = [], []
for ch in pd.read_csv(f, sep="\t", index_col=0, chunksize=3000):
    out.append(pd.DataFrame(np.asarray(sp.csr_matrix(ch.values) @ S.toarray() if False else ch.values @ S.toarray()), index=ch.index, columns=uniq)); print(ch.index[0], end=" ", flush=True)
pb = pd.concat(out); pb = pb[~pb.index.duplicated()]
pb.to_csv("data_clean/GSE298243_pb_counts.csv")
m = pd.DataFrame({"group_id": uniq, "n_cells": ncell.reindex(uniq).values}); m["donor"] = m.group_id.str.split("|").str[0]; m["fraction"] = m.group_id.str.split("|").str[1]
m = m.merge(meta, on="donor", how="left"); m.to_csv("data_clean/GSE298243_pb_meta.csv", index=False)
print("\n", m.groupby(["fraction", "diagnosis"]).agg(n_pb=("donor", "size"), median_cells=("n_cells", "median")).to_string(), "\ndonors:", m.donor.nunique())
