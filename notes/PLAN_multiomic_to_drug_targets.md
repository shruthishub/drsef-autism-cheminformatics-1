# Plan: brain multi-omics -> target nomination -> cheminformatics (drafted 2026-10-19)

Status: PLAN ONLY. Nothing below has been run except the metadata checks in week0_findings.md.
Interferon/immune is no longer a hypothesis; it is read from the data like any other pathway.

## Stage 1. Expression, per region (frontal, cerebellum)
Cohorts: GSE28521 (frontal/temporal/cerebellum), GSE102741, GSE38322, GSE64018 (sensitivity), GSE269105 (cerebellum, 7-8 v 5), Velmeshev pseudobulk.
Output per gene per region: log fold change, SE, class = UP / DOWN / UNCHANGED.
UNCHANGED must mean "confidently flat" (small |LFC| with a tight CI, equivalence test), not just "not significant".
Existing scripts 03-08 are reusable. Donor overlap (GSE28521/GSE38322/GSE64018) stays flagged.

## Stage 2. Methylation, per region
Cohorts: GSE53162 (PFC/temporal/cerebellum), GSE53924 (BA10/BA24, repeated donors), GSE278285 (cerebellum), GSE109875 (BA9 WGBS), GSE38608 (cerebellum/occipital, same donors as GSE38322), GSE298243 (BA9 snmC, same donors as the snRNA counts).
450K: minfi QC, cell-composition/age/sex covariates, probe -> region (DMR). Output class = HYPER / HYPO / UNCHANGED.
GSE107012 is 5hmC, a different mark: keep separate.

## Stage 3. Integration (activated vs deactivated)
Per region, per gene: promoter methylation direction x expression direction.
- hyper + down = candidate silenced; hypo + up = candidate activated
- discordant = reported as its own group
Same-donor pairs (GSE38608 x GSE38322; GSE298243) allow person-level tests. All other cohort pairs are group-level only.

## Stage 4. "Mutated genes"
(a) Overlap: enrich UP/DOWN/DMR genes for SFARI (data_raw/reference/sfari_genes.csv) and published de novo / rare-variant gene lists (to be sourced and cited, not typed from memory).
(b) Variant calling from raw RNA-seq (SRA). Needs: sra-tools, aligner (STAR/HISAT2), samtools/bcftools or GATK. NONE are installed. Only expressed genes, noisy, no matched controls for germline vs somatic.
    Feasible scope on this Mac (39 GB free, no conda/brew): GSE269105 only (13 samples), if tools are installed.
    GSE102741 (52) + GSE64018 (24) are too big for local disk; need a cluster/cloud or a smaller subset.
    Any variant call is a CANDIDATE, not a validated mutation.

## Stage 5. Target nomination
Rank genes that are (1) strongly UP in at least one region, (2) direction-consistent across cohorts where possible, (3) have supporting methylation or SFARI evidence, (4) are druggable (protein with a known inhibitable class, not a transcription factor with no pocket).
Inhibitors only make sense for UP genes. DOWN genes need activators/replacement, so are out of scope for this stage.
Weakness to state up front: replication across cohorts has been poor (see week0_findings.md), so the target list will be hypotheses, not findings.

## Stage 6. Cheminformatics
1. Known drugs: query ChEMBL / DrugBank / DGIdb / Open Targets for approved or clinical compounds against each target.
2. Docking: protein structure from PDB or AlphaFold; Vina-type docking with a redocking control (known ligand must be recovered, RMSD < 2 A) and a decoy set before trusting any score.
3. Output = ranked repurposing candidates with binding-score caveats.
4. "New therapeutic target and dock them": novel-compound design is a separate, later step. Needs medicinal-chemistry input; docking scores alone do not show a compound inhibits a target or helps autism.
Tools missing: RDKit, AutoDock Vina, Open Babel, ChemmineR (nothing installed).

## Decisions needed from user
1. Variant calling (b): scope to GSE269105 only, or get compute for the bigger RNA-seq sets?
2. OK to install tools (sra-tools, STAR, samtools, bcftools; RDKit, Vina, Open Babel)? Needs a package manager (conda/brew not present).
3. Source for de novo / rare-variant gene lists for (a).
