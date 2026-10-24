"""15_dock_panel.py
Dock every ligand in docking/panel.csv into one prepared receptor with AutoDock Vina and score the
enrichment of known actives over decoys.
Run:  ~/miniforge3/envs/docking/bin/python scripts/15_dock_panel.py TARGET_NAME DOCK_DIR
  e.g. scripts/15_dock_panel.py MKNK2 docking/8P9B_X8K [PANEL_CSV] [EXHAUSTIVENESS] [BOX_A]
Needs DOCK_DIR/receptor.pdbqt and DOCK_DIR/crystal_ligand.pdb (made by script 13).
Output: results/15_<TARGET>_panel_scores.csv  and a printed summary.
Vina scores grow with molecule size, so ligand efficiency (score / heavy atoms) is reported too.
"""
import sys, os, numpy as np, pandas as pd
from multiprocessing import Pool
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from vina import Vina
RDLogger.DisableLog("rdApp.*")
target, ddir = sys.argv[1], sys.argv[2].rstrip("/")
EXH = int(sys.argv[4]) if len(sys.argv) > 4 else 16
BOX = float(sys.argv[5]) if len(sys.argv) > 5 else 24.0
PANEL = sys.argv[3] if len(sys.argv) > 3 else "docking/panel.csv"
NPOSE, SEED = 5, 42
centre = Chem.MolFromPDBFile(f"{ddir}/crystal_ligand.pdb").GetConformer().GetPositions().mean(axis=0)

def dock(row):
    idx, smi = row
    try:
        m = Chem.AddHs(Chem.MolFromSmiles(smi))
        if AllChem.EmbedMolecule(m, randomSeed=SEED) != 0 and AllChem.EmbedMolecule(m, randomSeed=SEED, useRandomCoords=True) != 0:
            return idx, np.nan, "embed_failed"
        AllChem.MMFFOptimizeMolecule(m, maxIters=500)
        pdbqt, ok, err = PDBQTWriterLegacy.write_string(MoleculePreparation().prepare(m)[0])
        if not ok: return idx, np.nan, "meeko_failed"
        v = Vina(sf_name="vina", cpu=2, seed=SEED, verbosity=0)
        v.set_receptor(f"{ddir}/receptor.pdbqt"); v.set_ligand_from_string(pdbqt)
        v.compute_vina_maps(center=centre.tolist(), box_size=[BOX, BOX, BOX])
        v.dock(exhaustiveness=EXH, n_poses=NPOSE)
        return idx, float(v.energies(n_poses=1)[0][0]), "ok"
    except Exception as e:
        return idx, np.nan, f"error:{type(e).__name__}"

if __name__ == "__main__":
    p = pd.read_csv(PANEL)
    # checkpoint: every finished ligand is appended to a partial file, so a sleep/crash/kill loses almost nothing
    part = f"results/15_{target}_partial.csv"
    done = pd.read_csv(part, header=None, names=["i", "score", "status"]) if os.path.exists(part) else pd.DataFrame(columns=["i", "score", "status"])
    todo = [(i, s) for i, s in zip(p.index, p.smiles) if i not in set(done.i.astype(int))]
    print(f"resuming: {len(done)} done, {len(todo)} to go", flush=True)
    with Pool(4) as pool, open(part, "a") as fh:
        for n, (i, sc, st) in enumerate(pool.imap_unordered(dock, todo), 1):
            fh.write(f"{i},{sc},{st}\n"); fh.flush()
            if n % 10 == 0: print(f"  {len(done) + n}/{len(p)} ligands docked", flush=True)
    r = pd.read_csv(part, header=None, names=["i", "score", "status"]).drop_duplicates("i", keep="last").set_index("i")
    p["score"], p["status"] = r.score, r.status
    p["LE"] = p.score / p.heavy                      # kcal/mol per heavy atom (more negative = more efficient)
    out = f"results/15_{target}_panel_scores.csv"; p.to_csv(out, index=False)
    ok = p[p.status == "ok"]
    print("docked:", len(ok), "of", len(p), "| failures:", p.status.value_counts().drop("ok", errors="ignore").to_dict())
    a, d = ok[ok.group == "1_MNK_inhibitor"], ok[ok.group == "decoy"]
    from scipy.stats import mannwhitneyu
    def auc(x, y):                                      # P(score of active better than decoy); lower score = better
        return float(np.mean([(xi < yj) + 0.5 * (xi == yj) for xi in x for yj in y]))
    for lab, col in [("raw Vina score", "score"), ("ligand efficiency", "LE")]:
        print(f"{lab:>18}: AUC = {auc(a[col], d[col]):.2f} | actives mean {a[col].mean():.2f} | decoys mean {d[col].mean():.2f} | MWU p = {mannwhitneyu(a[col], d[col]).pvalue:.1e}")
    # early enrichment: actives among top 10% by score
    top = ok.sort_values("score").head(int(0.1 * len(ok)))
    print(f"top 10% by score ({len(top)} ligands): {int((top.group=='1_MNK_inhibitor').sum())} actives, {int((top.group=='decoy').sum())} decoys")
    print("\nGroup means (score / LE):"); print(ok.groupby("group")[["score", "LE"]].mean().round(2).to_string())
    print("\nBrain drugs + kinase/target drugs, ranked by score:")
    print(ok[ok.group.isin(["2_brain_drug", "3_kinase_or_target_drug"])].sort_values("score")[["group", "name", "heavy", "score", "LE"]].round(2).to_string(index=False))
