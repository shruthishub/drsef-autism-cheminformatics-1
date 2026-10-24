"""23_genetics_first_druggability.py
Genetics-first view: start from genes with MUTATION support instead of expression, then ask which are druggable.
Universe = SFARI score-1 genes (245)  U  genes with genome-wide de novo LoF excess in denovo-db (Poisson p<1e-3, script 18)  U  the 13 expression candidates with a de novo flag (script 21, flag D).
For each gene (public APIs, no login; query date stored): UniProt accession/protein/EC/keywords/PDB count, ChEMBL single-protein target, drugs with a recorded mechanism
(action type, max phase), number of potent compounds (pChEMBL>=6).
IMPORTANT interpretation limit stored in the output: most of these genes are hit by LoF variants, i.e. LOSS of function. For them an INHIBITOR is the wrong direction
(an activator/agonist/positive modulator or gene therapy would be needed); inhibitor docking only makes sense if gain-of-function / over-activity is plausible.
Run: ~/miniforge3/envs/docking/bin/python scripts/23_genetics_first_druggability.py
"""
import requests, datetime, pandas as pd, numpy as np
from concurrent.futures import ThreadPoolExecutor
s = pd.read_csv("data_raw/reference/sfari_genes.csv").rename(columns={"gene-symbol": "gene", "gene-score": "sfari_score"})
dn = pd.read_csv("results/18_denovo_gene_counts.csv"); lof = dn[dn["class"] == "LoF"].set_index("gene"); mis = dn[dn["class"] == "missense"].set_index("gene")
M21 = pd.read_csv("results/21_candidate_evidence_matrix.csv")
src = {}
for g in s[s.sfari_score == 1].gene: src.setdefault(g, []).append("SFARI1")
for g in lof.index[lof.poisson_p_excess < 1e-3]: src.setdefault(g, []).append("denovoLoF_p<1e-3")
for g in M21[M21.D].gene: src.setdefault(g, []).append("expr_candidate+denovo")
genes = sorted(src); print("universe:", len(genes))
ses = requests.Session()
def J(url):
    for _ in range(3):
        try:
            r = ses.get(url, timeout=40)
            if r.ok: return r.json()
        except Exception: pass
    return None
def one(g):
    r = dict(gene=g, source=";".join(src[g]))
    u = J(f"https://rest.uniprot.org/uniprotkb/search?query=gene_exact:{g}+AND+organism_id:9606+AND+reviewed:true&fields=accession,protein_name,ec,keyword,xref_pdb&format=json&size=1")
    if not u or not u.get("results"): return r
    x = u["results"][0]; r["uniprot"] = x["primaryAccession"]
    r["protein"] = x.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value")
    r["keywords"] = ";".join(k["name"] for k in x.get("keywords", []) if k["name"] in ("Kinase", "Receptor", "Ion channel", "Transcription regulation", "Hydrolase", "Transferase", "Transport", "Chromatin regulator", "Ubl conjugation", "Protease", "Oxidoreductase", "G-protein coupled receptor", "Ligand-gated ion channel", "Voltage-gated channel", "Methyltransferase", "Helicase", "Phosphatase"))
    r["n_pdb"] = sum(1 for c in x.get("uniProtKBCrossReferences", []) if c["database"] == "PDB")
    t = J(f"https://www.ebi.ac.uk/chembl/api/data/target.json?target_components__accession={r['uniprot']}&target_type=SINGLE%20PROTEIN&limit=1")
    if t and t.get("targets"):
        tid = t["targets"][0]["target_chembl_id"]; r["chembl_target"] = tid
        me = J(f"https://www.ebi.ac.uk/chembl/api/data/mechanism.json?target_chembl_id={tid}&limit=100")
        if me:
            ms = me["mechanisms"]; r["n_mech"] = len(ms); ph = [float(m_["max_phase"]) for m_ in ms if m_.get("max_phase") is not None]; r["max_phase"] = max(ph) if ph else np.nan
            r["action_types"] = ";".join(sorted({m_["action_type"] for m_ in ms if m_.get("action_type")})); r["mech_molecules"] = ";".join(sorted({m_["molecule_chembl_id"] for m_ in ms}))[:300]
        a = J(f"https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id={tid}&pchembl_value__gte=6&limit=1")
        r["n_potent"] = a["page_meta"]["total_count"] if a else np.nan
    return r
with ThreadPoolExecutor(8) as ex: rows = list(ex.map(one, genes))
D = pd.DataFrame(rows)
D["sfari_score"] = D.gene.map(s.drop_duplicates("gene").set_index("gene").sfari_score); D["syndromic"] = D.gene.map(s.drop_duplicates("gene").set_index("gene").syndromic)
D["denovo_LoF_autism"] = D.gene.map(lof.autism_obs); D["denovo_LoF_p"] = D.gene.map(lof.poisson_p_excess); D["denovo_mis_autism"] = D.gene.map(mis.autism_obs); D["denovo_mis_p"] = D.gene.map(mis.poisson_p_excess)
D = D.merge(M21[["gene", "z_meta", "n_layers"]].rename(columns={"z_meta": "discovery_z_if_candidate"}), on="gene", how="left")
D["druggable"] = (D.n_mech.fillna(0) > 0) | (D.n_potent.fillna(0) >= 10)
D["lof_driven"] = D.denovo_LoF_autism.fillna(0) >= D.denovo_mis_autism.fillna(0)
D["query_date"] = str(datetime.date.today())
D.to_csv("results/23_genetics_first_druggability.csv", index=False)
print(f"queried {len(D)} genes | UniProt found {D.uniprot.notna().sum()} | with PDB {int((D.n_pdb>0).sum())} | ChEMBL target {D.chembl_target.notna().sum()} | druggable {int(D.druggable.sum())} | with approved/clinical mechanism drug {int((D.n_mech>0).sum())}")
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 45)
print("\nDruggable genes with genetic support:"); print(D[D.druggable].sort_values(["n_potent"], ascending=False)[["gene", "source", "sfari_score", "denovo_LoF_autism", "denovo_mis_autism", "keywords", "n_pdb", "n_mech", "max_phase", "action_types", "n_potent", "lof_driven"]].to_string(index=False))
