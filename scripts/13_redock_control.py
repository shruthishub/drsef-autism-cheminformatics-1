"""13_redock_control.py
Re-docking control: take a crystal structure that has a ligand, remove the ligand, dock it back
from a RANDOM starting shape with AutoDock Vina, and compare to the crystal pose (heavy-atom RMSD,
symmetry aware). A target only counts as dockable if the crystal pose is recovered (RMSD < 2 A).

Run with the docking environment:
  ~/miniforge3/envs/docking/bin/python scripts/13_redock_control.py PDBID LIGAND_RESNAME [CHAIN] [BOX_SIZE_A]
Outputs go to docking/<PDBID>_<LIG>/ and one summary line is appended to results/13_redock_controls.csv
"""
import sys, os, json, subprocess, urllib.request, csv, ssl
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from meeko import MoleculePreparation, PDBQTWriterLegacy, PDBQTMolecule, RDKitMolCreate
from vina import Vina

pdb_id, lig = sys.argv[1].upper(), sys.argv[2].upper()
chain = sys.argv[3] if len(sys.argv) > 3 else "A"
BOX = float(sys.argv[4]) if len(sys.argv) > 4 else 24.0   # cube edge in Angstrom; use ~30 for long type-II inhibitors
OB = os.path.expanduser("~/miniforge3/envs/docking/bin/obabel")
out = f"docking/{pdb_id}_{lig}"; os.makedirs(out, exist_ok=True)
SEED = int(os.environ.get("SEED", 42))          # override: SEED=7 EXH=64 python ...
EXH = int(os.environ.get("EXH", 32))

def fetch(url, path):
    if not os.path.exists(path):
        r = subprocess.run(["curl", "-s", "-L", "-o", path, url], check=True)
    return path

pdb = fetch(f"https://files.rcsb.org/download/{pdb_id}.pdb", f"data_raw/structures/{pdb_id}.pdb")
ccd = json.load(open(fetch(f"https://data.rcsb.org/rest/v1/core/chemcomp/{lig}", f"data_raw/structures/{lig}_ccd.json")))
desc = ccd["rcsb_chem_comp_descriptor"]
smiles = desc.get("SMILES_stereo") or desc.get("SMILES")
print("ligand", lig, "SMILES:", smiles)

# ---- split protein / ligand from the crystal file ----
prot, ligl = [], []
for ln in open(pdb):
    if ln.startswith("ATOM") and ln[21] == chain:
        prot.append(ln)
    elif ln.startswith("HETATM") and ln[17:20].strip() == lig and ln[21] == chain:
        ligl.append(ln)
assert prot and ligl, "protein or ligand not found"
open(f"{out}/receptor.pdb", "w").write("".join(prot) + "END\n")
open(f"{out}/crystal_ligand.pdb", "w").write("".join(ligl) + "END\n")

# ---- crystal ligand with correct bond orders (template = CCD SMILES) ----
raw = Chem.MolFromPDBFile(f"{out}/crystal_ligand.pdb", removeHs=True)
tmpl = Chem.MolFromSmiles(smiles)
ref = AllChem.AssignBondOrdersFromTemplate(tmpl, raw)
ref = Chem.RemoveHs(ref)
centre = ref.GetConformer().GetPositions().mean(axis=0)
print("heavy atoms:", ref.GetNumHeavyAtoms(), "| pocket centre:", np.round(centre, 2))

# ---- receptor -> pdbqt (Open Babel: adds H, charges; rigid) ----
subprocess.run([OB, f"{out}/receptor.pdb", "-O", f"{out}/receptor.pdbqt", "-xr", "-h"], check=True, stderr=subprocess.DEVNULL)

# ---- ligand: random 3D conformer from SMILES (NOT the crystal pose) ----
lm = Chem.AddHs(Chem.MolFromSmiles(smiles))
assert AllChem.EmbedMolecule(lm, randomSeed=SEED) == 0
AllChem.MMFFOptimizeMolecule(lm)
setup = MoleculePreparation().prepare(lm)[0]
pdbqt, ok, err = PDBQTWriterLegacy.write_string(setup)
assert ok, err

v = Vina(sf_name="vina", seed=SEED, verbosity=0)
v.set_receptor(f"{out}/receptor.pdbqt")
v.set_ligand_from_string(pdbqt)
v.compute_vina_maps(center=centre.tolist(), box_size=[BOX, BOX, BOX])
v.dock(exhaustiveness=EXH, n_poses=9)
v.write_poses(f"{out}/poses.pdbqt", n_poses=9, overwrite=True)
energies = v.energies(n_poses=9)[:, 0]

pm = PDBQTMolecule.from_file(f"{out}/poses.pdbqt", skip_typing=True)
mols = RDKitMolCreate.from_pdbqt_mol(pm)[0]
rms = []
for i in range(mols.GetNumConformers()):
    conf = Chem.Mol(mols); conf.RemoveAllConformers(); conf.AddConformer(Chem.Conformer(mols.GetConformer(i)), assignId=True)
    conf = Chem.RemoveHs(conf)
    rms.append(rdMolAlign.CalcRMS(conf, ref))
rms = np.array(rms)
print("\npose  vina_kcal/mol  RMSD_to_crystal_A")
for i, (e, r) in enumerate(zip(energies, rms), 1):
    print(f"{i:>4}  {e:>12.2f}  {r:>10.2f}")
top_rmsd, best_rmsd = rms[0], rms.min()
verdict = "PASS" if top_rmsd < 2 else ("PARTIAL (a lower-ranked pose is < 2 A)" if best_rmsd < 2 else "FAIL")
print(f"\nTop-scored pose RMSD {top_rmsd:.2f} A | best of 9 {best_rmsd:.2f} A | verdict: {verdict}")

os.makedirs("results", exist_ok=True)
f = "results/13_redock_controls.csv"; new = not os.path.exists(f)
with open(f, "a", newline="") as fh:
    w = csv.writer(fh)
    if new: w.writerow(["pdb", "ligand", "chain", "heavy_atoms", "top_score", "top_rmsd", "best_rmsd", "verdict", "seed", "exhaustiveness"])
    w.writerow([pdb_id, lig, chain, ref.GetNumHeavyAtoms(), round(float(energies[0]), 2), round(float(top_rmsd), 2), round(float(best_rmsd), 2), verdict, SEED, EXH])
