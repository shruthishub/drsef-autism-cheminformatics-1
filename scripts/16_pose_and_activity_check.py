"""16_pose_and_activity_check.py
Follow-up to the MKNK2 panel docking (scripts 14-15). For selected drugs:
 (a) look up MEASURED activity against MNK2 (CHEMBL4204) and MNK1 (CHEMBL4718) in ChEMBL
 (b) re-dock into 8P9B keeping 9 poses, then ask whether any pose makes the kinase hinge hydrogen bond
     (ligand N/O within 3.4 A of the hinge backbone atoms that the CRYSTAL ligand contacts) and sits in the ATP site
Controls: tinodasertib + staurosporine (real kinase inhibitors), two top-scoring DECOYS (eltrombopag, delamanid).
Run:  ~/miniforge3/envs/docking/bin/python scripts/16_pose_and_activity_check.py
"""
import os, time, numpy as np, pandas as pd, requests
from multiprocessing import Pool
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy, PDBQTMolecule, RDKitMolCreate
from vina import Vina
RDLogger.DisableLog("rdApp.*")
# Settings come from environment variables so the same script works for any target (defaults = MKNK2 run):
#   DDIR=docking/4BKJ_STI PANEL=docking/panel_DDR1.csv BOX=30 TAG=DDR1 TARGETS=CHEMBL5319:DDR1 NAMES="aripiprazole,risperidone,..."
D = os.environ.get("DDIR", "docking/8P9B_X8K"); SEED, EXH, NP = 42, 16, 9
BOX = float(os.environ.get("BOX", 24)); TAG = os.environ.get("TAG", "")
TARGETS = [t.split(":") for t in os.environ.get("TARGETS", "CHEMBL4204:MNK2,CHEMBL4718:MNK1").split(",")]
NAMES = os.environ.get("NAMES", "lurasidone,risperidone,paliperidone,buspirone,ziprasidone,donepezil,tinodasertib,staurosporine,ELTROMBOPAG,DELAMANID").split(",")
panel = pd.read_csv(os.environ.get("PANEL", "docking/panel.csv"))
sel = panel[panel.name.str.lower().isin([n.lower() for n in NAMES])].drop_duplicates("name").reset_index(drop=True)

# ---- protein + crystal ligand geometry ----
def atoms(path, rec):
    out = []
    for ln in open(path):
        if ln.startswith(rec):
            out.append((ln[12:16].strip(), ln[17:20].strip(), int(ln[22:26]), ln[21], np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]), ln[76:78].strip()))
    return out
prot = atoms(f"{D}/receptor.pdb", "ATOM"); xl = atoms(f"{D}/crystal_ligand.pdb", "HETATM")
bb = [a for a in prot if a[0] in ("N", "O")]                       # backbone H-bond atoms
xl_no = [a for a in xl if a[5] in ("N", "O")]
hinge = {}
for l in xl_no:
    for b in bb:
        if np.linalg.norm(l[4] - b[4]) <= 3.4: hinge[(b[1], b[2], b[0])] = b[4]
print("Crystal ligand H-bond partners on the protein BACKBONE (taken as hinge):", sorted(hinge))
xc = np.mean([a[4] for a in xl], axis=0)

# ---- (a) ChEMBL measured activity ----
def chembl(path, **p):
    for i in range(4):
        time.sleep(0.3 * (i + 1))
        try:
            r = requests.get(f"https://www.ebi.ac.uk/chembl/api/data/{path}.json", params=p, timeout=60)
            if r.ok: return r.json()
        except Exception: pass
    return None
rows = []
for r in sel.itertuples():
    for tid, tn in TARGETS:
        j = chembl("activity", molecule_chembl_id=r.chembl_id, target_chembl_id=tid, limit=50)
        n = None if j is None else j["page_meta"]["total_count"]
        vals = [] if not j else [f'{a["standard_type"]} {a["standard_relation"] or ""}{a["standard_value"]} {a["standard_units"]}' for a in j["activities"] if a.get("standard_value")]
        rows.append((r.name, tn, n, "; ".join(vals[:3])))
act = pd.DataFrame(rows, columns=["drug", "target", "n_records", "examples"])
print("\n(a) ChEMBL measured activity records (None = query failed):"); print(act.to_string(index=False))

# ---- (b) pose check ----
def dock(row):
    name, smi = row
    m = Chem.AddHs(Chem.MolFromSmiles(smi)); AllChem.EmbedMolecule(m, randomSeed=SEED, useRandomCoords=True); AllChem.MMFFOptimizeMolecule(m, maxIters=500)
    pq, ok, _ = PDBQTWriterLegacy.write_string(MoleculePreparation().prepare(m)[0])
    v = Vina(sf_name="vina", cpu=2, seed=SEED, verbosity=0)
    v.set_receptor(f"{D}/receptor.pdbqt"); v.set_ligand_from_string(pq); v.compute_vina_maps(center=xc.tolist(), box_size=[BOX, BOX, BOX])
    v.dock(exhaustiveness=EXH, n_poses=NP); f = f"{D}/poses_{name.replace(' ', '_')}.pdbqt"; v.write_poses(f, n_poses=NP, overwrite=True)
    return name, f, v.energies(n_poses=NP)[:, 0].tolist()

if __name__ == "__main__":
    with Pool(4) as pool: res = pool.map(dock, list(zip(sel.name, sel.smiles)), chunksize=1)
    out = []
    for name, f, en in res:
        mol = RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule.from_file(f, skip_typing=True))[0]
        for k in range(mol.GetNumConformers()):
            pos = mol.GetConformer(k).GetPositions(); el = [a.GetSymbol() for a in mol.GetAtoms()]
            heavy = [i for i, e in enumerate(el) if e in ("N", "O")]
            hb = sorted({key[0] + str(key[1]) + ":" + key[2] for key, bpos in hinge.items() for i in heavy if np.linalg.norm(pos[i] - bpos) <= 3.4})
            cen = pos[[i for i, e in enumerate(el) if e != "H"]].mean(axis=0)
            out.append((name, k + 1, round(en[k], 2), len(hb), ",".join(hb), round(float(np.linalg.norm(cen - xc)), 1)))
    P = pd.DataFrame(out, columns=["drug", "pose", "score", "n_hinge_hbonds", "hinge_partners", "centroid_dist_to_xtal_A"])
    P.to_csv(f"results/16_pose_check_all_poses{TAG}.csv", index=False)
    print("\n(b) For each drug: top-scored pose, and the best-scoring pose that makes >=1 hinge H-bond inside the ATP site (centroid <5 A of crystal ligand):")
    for name in NAMES:
        q = P[P.drug.str.upper() == name.upper()]
        if q.empty: continue
        top = q.iloc[0]; hb = q[(q.n_hinge_hbonds >= 1) & (q.centroid_dist_to_xtal_A < 5)]
        s = f"  {name:<14} top pose {top.score:>6} (hinge H-bonds: {top.n_hinge_hbonds}, {top.centroid_dist_to_xtal_A} A from site)"
        s += f" | best hinge pose: rank {int(hb.iloc[0].pose)} score {hb.iloc[0].score} [{hb.iloc[0].hinge_partners}]" if len(hb) else " | NO pose with a hinge H-bond in the ATP site"
        print(s)
    act.to_csv(f"results/16_chembl_activity_check{TAG}.csv", index=False)
