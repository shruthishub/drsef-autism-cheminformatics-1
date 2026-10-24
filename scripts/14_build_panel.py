"""14_build_panel.py
Build the ligand panel for docking against MKNK2 (ChEMBL target CHEMBL4204), from ChEMBL (queried live).
 group 1  known MNK inhibitors   : potent ChEMBL actives (pChEMBL >= 7.5), max 3 per Murcko scaffold, up to 50 (POSITIVE CONTROLS)
 group 2  brain-disorder drugs   : schizophrenia / bipolar / autism-irritability / Down-syndrome-trial drugs, by name
 group 3  approved kinase inhibitors: from results/12_target_drug_pairs.csv (small molecules) + named extras
 decoys   : phase-4 small-molecule drugs, NOT kinase inhibitors, property-matched (MW, logP, HBD, HBA, rotatable bonds)
            to the group-1 actives, 4 per active, no reuse. APPROXIMATION of DUD-E style decoys (those are not available offline).
Run:  ~/miniforge3/envs/docking/bin/python scripts/14_build_panel.py [CHEMBL_TARGET_ID NAME N_ACTIVES DECOYS_PER_ACTIVE]
      default = CHEMBL4204 MKNK2 50 4
Output: docking/panel.csv (MKNK2) or docking/panel_<NAME>.csv   (columns: group, name, chembl_id, smiles, pchembl, MW, logP, HBD, HBA, RotB, heavy)
"""
import re, sys, time, json, datetime
import numpy as np, pandas as pd, requests
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
RDLogger.DisableLog("rdApp.*")
B = "https://www.ebi.ac.uk/chembl/api/data"
TID = sys.argv[1] if len(sys.argv) > 1 else "CHEMBL4204"
TNAME = sys.argv[2] if len(sys.argv) > 2 else "MKNK2"
NACT = int(sys.argv[3]) if len(sys.argv) > 3 else 50
NDEC = int(sys.argv[4]) if len(sys.argv) > 4 else 4
from rdkit.Chem import inchi
def pkey(m): return inchi.MolToInchiKey(m)[:14]            # connectivity block = parent compound (salts/stereo ignored)

def get(path, **params):
    for i in range(4):
        try:
            time.sleep(0.25 * (i + 1))
            r = requests.get(f"{B}/{path}.json", params=params, timeout=90)
            if r.ok: return r.json()
        except Exception: pass
    raise RuntimeError(f"ChEMBL failed: {path} {params}")

def clean(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None: return None
    frags = Chem.GetMolFrags(m, asMols=True)
    m = max(frags, key=lambda x: x.GetNumHeavyAtoms())          # strip salts / counter-ions
    return m

def props(m):
    return dict(MW=Descriptors.MolWt(m), logP=Crippen.MolLogP(m), HBD=Lipinski.NumHDonors(m),
                HBA=Lipinski.NumHAcceptors(m), RotB=rdMolDescriptors.CalcNumRotatableBonds(m), heavy=m.GetNumHeavyAtoms())

rows = []
# ---------- group 1 ----------
acts = get("activity", target_chembl_id=TID, pchembl_value__gte=7.5, limit=1000)["activities"]
g1 = {}
for a in acts:
    s, mid = a.get("canonical_smiles"), a["molecule_chembl_id"]
    if not s: continue
    p = float(a["pchembl_value"]); 
    if mid not in g1 or p > g1[mid][1]: g1[mid] = (s, p, a.get("molecule_pref_name") or mid)
scaf_n, picked = {}, []
for mid, (s, p, nm) in sorted(g1.items(), key=lambda kv: -kv[1][1]):
    m = clean(s)
    if m is None or not (150 < Descriptors.MolWt(m) < 650): continue
    sc = MurckoScaffold.MurckoScaffoldSmiles(mol=m)
    if scaf_n.get(sc, 0) >= 3: continue
    scaf_n[sc] = scaf_n.get(sc, 0) + 1
    picked.append(("1_MNK_inhibitor", nm, mid, Chem.MolToSmiles(m), p))
    if len(picked) == NACT: break
rows += picked
print("group 1 actives:", len(picked), "(from", len(g1), "potent compounds)")

# ---------- group 2 & named extras ----------
def by_name(name):
    j = get("molecule/search", q=name, limit=10)["molecules"]
    for mol in j:
        if (mol.get("pref_name") or "").lower() == name.lower() and mol.get("molecule_structures"):
            return mol
    return None
brain = ["risperidone","aripiprazole","clozapine","olanzapine","quetiapine","haloperidol","ziprasidone","lurasidone","paliperidone",
         "valproic acid","lamotrigine","carbamazepine","memantine","donepezil","galantamine","rivastigmine","lithium carbonate",
         "methylphenidate","guanfacine","atomoxetine","fluoxetine","sertraline","buspirone","epigallocatechin gallate","harmine","leucettine l41"]
for nm in brain:
    mol = by_name(nm)
    if mol is None: print("  not found in ChEMBL by exact name:", nm); continue
    m = clean(mol["molecule_structures"]["canonical_smiles"])
    rows.append(("2_brain_drug", nm, mol["molecule_chembl_id"], Chem.MolToSmiles(m), np.nan))

# ---------- group 3 ----------
pairs = pd.read_csv("results/12_target_drug_pairs.csv")
ids = [i for i in pairs["molecule"].unique()]
seen = {r[2] for r in rows}
for mid in ids:
    j = get(f"molecule/{mid}")
    if j.get("molecule_type") != "Small molecule" or not j.get("molecule_structures"): continue   # skip antibodies/proteins
    m = clean(j["molecule_structures"]["canonical_smiles"]);  nm = j.get("pref_name") or mid
    if mid not in seen: rows.append(("3_kinase_or_target_drug", nm, mid, Chem.MolToSmiles(m), np.nan))
for nm in ["tinodasertib", "imatinib", "dasatinib", "staurosporine", "nilotinib", "ponatinib"]:
    mol = by_name(nm)
    if mol and mol["molecule_chembl_id"] not in {r[2] for r in rows}:
        rows.append(("3_kinase_or_target_drug", nm, mol["molecule_chembl_id"], Chem.MolToSmiles(clean(mol["molecule_structures"]["canonical_smiles"])), np.nan))

df = pd.DataFrame(rows, columns=["group","name","chembl_id","smiles","pchembl"])
P = pd.DataFrame([props(Chem.MolFromSmiles(s)) for s in df.smiles]); df = pd.concat([df, P], axis=1)
print(df.group.value_counts().to_string())

# ---------- decoys ----------
kin = re.compile(r"(tinib|nib$|ciclib|sertib|rafenib|degib|parib|lisib|limus|stat$|mab$|cept$)", re.I)
pool, off = [], 0
while True:
    j = get("molecule", max_phase=4, molecule_type="Small molecule", limit=1000, offset=off)
    pool += j["molecules"]
    if not j["page_meta"]["next"]: break
    off += 1000
print("phase-4 small molecules fetched:", len(pool))
used = {r.chembl_id for r in df.itertuples()}
used_keys = {pkey(Chem.MolFromSmiles(x)) for x in df.smiles}
cand = []
for mol in pool:
    nm = mol.get("pref_name") or ""
    if mol["molecule_chembl_id"] in used or not mol.get("molecule_structures") or kin.search(nm): continue
    m = clean(mol["molecule_structures"]["canonical_smiles"])
    if m is None or not (150 < Descriptors.MolWt(m) < 650): continue
    k = pkey(m)
    if k in used_keys: continue                      # same parent compound already in the panel or already a decoy
    used_keys.add(k)
    cand.append((nm or mol["molecule_chembl_id"], mol["molecule_chembl_id"], Chem.MolToSmiles(m), props(m)))
keys = ["MW","logP","HBD","HBA","RotB"]
C = np.array([[c[3][k] for k in keys] for c in cand]); sd = C.std(axis=0)
taken, dec = set(), []
for a in df[df.group == "1_MNK_inhibitor"].itertuples():
    v = np.array([getattr(a, k) for k in keys]); d = (np.abs(C - v) / sd).sum(axis=1)
    n = 0
    for i in np.argsort(d):
        if i in taken: continue
        taken.add(i); dec.append(("decoy", *cand[i][:3], np.nan, *[cand[i][3][k] for k in keys + ["heavy"]])); n += 1
        if n == NDEC: break
D = pd.DataFrame(dec, columns=["group","name","chembl_id","smiles","pchembl"] + keys + ["heavy"])
out = pd.concat([df, D], ignore_index=True)
outf = "docking/panel.csv" if TNAME == "MKNK2" else f"docking/panel_{TNAME}.csv"
out.to_csv(outf, index=False)
print("panel written:", outf, len(out), "ligands | queried", datetime.date.today())
print(out.groupby("group")[keys + ["heavy"]].mean().round(1).to_string())
