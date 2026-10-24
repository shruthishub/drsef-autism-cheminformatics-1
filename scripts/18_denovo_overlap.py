"""18_denovo_overlap.py
Mutation layer (option a): are the candidate genes hit by de novo coding variants in autism probands (denovo-db, open access)?
Test: observed de novo variants in autism probands vs the number expected from the gene's coding length (Poisson), separately for
loss-of-function (frameshift / stop-gained / splice) and missense. Controls shown for comparison. Length alone is a crude null (ignores
mutability, GC, expression), so this is a screen, not a rate test. Studies overlap between SSC and non-SSC files; variants are de-duplicated
on (sample, chr, position, variant).
Run:  ~/miniforge3/envs/docking/bin/python scripts/18_denovo_overlap.py SCRATCH_DIR     (SCRATCH_DIR has dn_ssc.tsv.gz, dn_non-ssc.tsv.gz)
"""
import sys, numpy as np, pandas as pd
from scipy import stats
S = sys.argv[1]
d = pd.concat([pd.read_csv(f"{S}/dn_{k}.tsv.gz", sep="\t", comment="#", header=None, low_memory=False,
               names="SampleID StudyName PubmedID NumProbands NumControls SequenceType PrimaryPhenotype Validation Chr Position Variant rsID DbsnpBuild AncestralAllele 1000GenomeCount ExacFreq EspAaFreq EspEaFreq Transcript codingDnaSize Gene FunctionClass cDnaVariant ProteinVariant ExonIntron PolyPhenHDiv PolyPhenHVar SiftScore CaddScore LofScore LrtScore".split())
               for k in ["ssc", "non-ssc"]]).drop_duplicates(["SampleID", "Chr", "Position", "Variant"])
print("phenotypes:", d.PrimaryPhenotype.value_counts().head(8).to_dict()); print("function classes:", d.FunctionClass.value_counts().head(14).to_dict())
cds = d[d.codingDnaSize > 0].groupby("Gene").codingDnaSize.max()          # coding length from ALL variants (any class) so zero-hit genes are included
d["cls"] = d.FunctionClass.map(lambda x: "LoF" if x in ("frameshift", "stop-gained", "splice-donor", "splice-acceptor", "start-lost", "stop-lost") else ("missense" if x == "missense" else None))
d = d[d.cls.notna() & d.PrimaryPhenotype.isin(["autism", "control"])]
nprob = d[d.PrimaryPhenotype == "autism"].drop_duplicates("SampleID").shape[0]; nctl = d[d.PrimaryPhenotype == "control"].drop_duplicates("SampleID").shape[0]
print(f"coding de novo variants kept: autism {int((d.PrimaryPhenotype=='autism').sum())} (from {nprob} probands with >=1 variant), control {int((d.PrimaryPhenotype=='control').sum())} ({nctl} controls with >=1 variant)")
rows = []
for cls in ["LoF", "missense"]:
    a = d[(d.PrimaryPhenotype == "autism") & (d.cls == cls)].groupby("Gene").size(); c = d[(d.PrimaryPhenotype == "control") & (d.cls == cls)].groupby("Gene").size()
    tot, L = a.sum(), cds.sum()
    for g in cds.index:
        ob = int(a.get(g, 0)); ex = tot * cds[g] / L
        rows.append((g, cls, ob, int(c.get(g, 0)), round(ex, 3), stats.poisson.sf(ob - 1, ex) if ob > 0 else 1.0))
R = pd.DataFrame(rows, columns=["gene", "class", "autism_obs", "control_obs", "autism_expected_by_length", "poisson_p_excess"])
R.to_csv("results/18_denovo_gene_counts.csv", index=False)
cand = pd.read_csv("results/12_candidates_full.csv"); keep = set(cand.gene)
q = R[R.gene.isin(keep)]
print("\nSanity - strongest LoF excess genome-wide (known autism genes should lead):"); print(R[R["class"]=="LoF"].sort_values("poisson_p_excess").head(8).to_string(index=False))
print("\nDDR1 / MKNK2 (coding length DDR1 %s, MKNK2 %s bp in table):" % (cds.get("DDR1"), cds.get("MKNK2"))); print(R[R.gene.isin(["DDR1", "MKNK2"])].to_string(index=False))
print(f"\nOf {len(keep)} candidates ({q.gene.nunique()} present in denovo-db), with >=1 autism de novo variant:")
for cls in ["LoF", "missense"]:
    z = q[(q["class"] == cls) & (q.autism_obs > 0)]; print(f"  {cls}: {z.gene.nunique()} genes hit; with Poisson p<0.01: {int((z.poisson_p_excess<0.01).sum())}")
print("\nCandidates with >=2 autism de novo LoF or p<0.01 (any class):"); w = q[(q.autism_obs >= 2) & (q["class"] == "LoF") | (q.poisson_p_excess < 0.01)].sort_values("poisson_p_excess"); print(w.head(15).to_string(index=False))
# whole-genome sanity: how many genes genome-wide have p<0.01 for LoF at this null?
lo = R[R["class"] == "LoF"]; print(f"\ngenome-wide sanity (LoF): {int((lo.poisson_p_excess<0.01).sum())} of {len(lo)} genes have p<0.01 (chance ~{0.01*len(lo):.0f} if null held)")
