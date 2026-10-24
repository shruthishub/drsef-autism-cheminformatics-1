# Week 0 findings (2026-10-13), from script 01_download_metadata.R

Counts below come from the GEO series-matrix files, not from the paper. Labels were guessed by regex and then checked by eye on the printed values.

## MISMATCH 1 (important): GSE59288 has no human controls
- GSE59288 = 34 autism PFC samples (GPL11154) + 39 CHIMPANZEE samples (GPL16809). Total 73, which matches "GEO lists 73".
- The "40 controls" in the plan are probably the 40 RHESUS samples in GSE51264 (GPL14954), not human controls.
- Human samples are in the companion GSE51264: 38 "human_XX" PFC samples (GPL11154, same platform as the autism samples).
- Likely real design: 34 autism (GSE59288) vs 38 human controls (GSE51264). This matches the "34/38" table in the plan.
- TO VERIFY: that all 38 human GSE51264 samples are healthy controls (check the paper or GEO text).
- Risk: cases and controls are in different GEO series, so there may be a batch effect. Check this in QC.

## GSE28521 (microarray GPL6883): 79 samples
- 39 autism / 40 control. By region: cerebellum 21, frontal 32, temporal 26.
- Plan said 58 cortex (29/29) after QC. The raw GEO table has more, probably technical replicates. Needs checking.

## GSE28475 (Chow 2012): two platforms
- GPL13388 (20 samples): only a technical series (reference RNA, heat-treatment tests). No diagnosis. Do not use.
- GPL6883 (123 samples): autism 52, control 61, 10 reference RNA. Includes technical replicates ("rep1", Frozen vs Fixed tissue, DASL vs IVT assay). No brain region field in the metadata.
- Not independent samples. Needs deduplication, and a check for overlap with GSE28521.

## GSE102741 (RNA-seq): 13 autism / 39 control, all DLPFC. Matches the plan.
- Has age, sex, RIN, pH, race and brainid columns. No postmortem interval column.

## Other
- WGCNA is not installed in R yet (needed for weeks 3-4).
- Python pandas is broken on this Mac (missing dateutil). Fix before the plotting steps.

## Update: discovery set build and QC (script 02)
- Downloaded: GSE59288_exp_mRNA.txt.gz (34 autism + 39 chimp; chimp dropped), GSE51264_rpkm_human.tsv.gz (38 human), GSE102741_log2RPKMcounts.xlsx (not yet loaded).
- Human only, 12,541 shared Ensembl genes (GSE51264 IDs have version suffixes, which were stripped). Final: 34 autism / 38 control. Both PFC.
- Ages: autism 2-60 y, control 0.005-61.5 y. Controls skew younger (median 11 vs 18 y). Sex: autism 10F/24M, control 12F/26M.
- BATCH PROBLEM: PC1 explains 58% of variance and correlates r = -0.99 with diagnosis (age r = -0.11). Cases and controls come from different GEO series, processed differently (the autism table has only 12,557 genes, the control table 58k). A DE test on this matrix would only find the processing difference.
- Conclusion: the GSE59288 + GSE51264 pair CANNOT be used for DE as it stands. It would have to be re-quantified from raw reads (SRA) with one pipeline.
- Outlier flagged: aut16 (mean sample correlation 0.50, others >= 0.62).

## GEO search for other human brain ASD data (E-utilities, 2026-10-13)
- Nearly all hits are organoids, iPSC models, mouse, or zebrafish. NOT usable.
- Human postmortem but not bulk-expression: GSE298243 (single-cell methylation/expression ASD PFC, 49 samples), GSE242427 (methylation array, 12), GSE278285 (cerebellum methylation, 75).
- Not found by this search but known: GSE64018 (human, STG, 12 vs 12, selection bias, kept only as a sensitivity check).
- PsychENCODE (Parikshak 2016, Gandal 2022) needs a Synapse login. Not downloaded.

## Update: more human brain cohorts (script 01b, 2026-10-13)
Downloaded to data_raw/. All human postmortem. Counts are from GEO / UCSC metadata tables.

| Cohort | Data | Region | Autism / Control | Notes |
|---|---|---|---|---|
| Velmeshev 2019 (UCSC cell browser "autism") | snRNA-seq, raw counts (rawMatrix.zip 755 MB) + meta.tsv | PFC + ACC | 15 / 16 INDIVIDUALS (22 / 19 samples; 52k / 53k cells) | Ages only 4-22 y. Pseudobulk per person per cell type. Public. |
| GSE298243 | snRNA counts from snmCT-seq (also methylation) | BA9 (PFC) | 26 / 23 donors, all male | Ages 2-60. A SEPARATE cohort, not a re-upload of Velmeshev. Controlled-access raw reads, GEO has processed counts. |
| GSE64018 | RNA-seq counts | Superior temporal gyrus | 12 / 12 | Samples picked for largest expression difference (selection bias). Sensitivity check only. |
| GSE38322 | Illumina HT12 array (raw report files) | Cerebellum 14/12, occipital 4/6 | 18 / 18 | NO frontal or temporal cortex. Secondary only. |

Not usable for this question: GSE113834 (poly(A) tail fractions), GSE236761 (translatome, BA19), GSE211154 (Purkinje cells only), GSE36315 (8 samples), GSE30573 (3 vs 3), GSE80017 (methylation).

## DONOR OVERLAP (same person in more than one dataset)
- GSE28521 and GSE64018: 11 donors in common. GSE28521 and GSE38322: 7. GSE38322 and GSE64018: 2.
- GSE298243 and Velmeshev: 9 donor IDs match (numeric IDs, so confirm in the papers that these are the same people). GSE298243 and GSE64018: 1.
- GSE102741 shares no donors with any other cohort I checked.
- GSE28475 has no donor ID in its metadata, so overlap with GSE28521 is still unchecked (52 autism and 61 control arrays, 26 of them Fixed tissue).
- Consequence: GSE28521, GSE38322 and GSE64018 cannot all be counted as independent replications of each other.

## Aim 1 + Aim 2 results (scripts 03-08) - GATE 1 REACHED
Design: each cohort analysed on its own, then combined. PRIMARY independent cohorts = GSE28521 cortex (58 samples, 29/29, 32 donors) and GSE102741 (13/39, DLPFC). Sensitivity = GSE64018 (shares 11 donors with GSE28521, selected for large differences). Third check = Velmeshev snRNA pseudobulk (15 ASD / 16 control people). GSE59288+GSE51264 pair excluded (diagnosis = batch). GSE28475 excluded (technical frozen/fixed comparison, no region or donor IDs, independence from GSE28521 unverifiable). data_clean/discovery_* files come from the excluded pair; do not use.

Aim 1 (DE):
- GSE28521 cortex: 1,241 genes FDR<0.05 (432 up, 809 down), inflation lambda 3.6. No age/sex/RIN/PMI in GEO, so adjusted with 5 surrogate variables + region + donor correlation (0.56).
- GSE102741: 1 gene FDR<0.05 (lambda 0.32). Adjusted for age, sex, RIN, race.
- Velmeshev total pseudobulk: 0 genes FDR<0.05 (lambda 1.5).
- Replication GSE28521 -> GSE102741: 59.8% of 911 discovery hits point the same way (binomial p=3e-9), but 0 are p<0.05 in the replication cohort (1.1 expected by chance). Spearman r of z across genes = 0.18. Weak directional agreement, no gene-level replication.
- Velmeshev vs GSE28521: r=0.05, 52.5% same direction (p=0.12). Velmeshev vs GSE102741: r=0.055.
- Primary meta-analysis (weighted Stouffer): 3 genes FDR<0.05, 0 fully replicated. Adding GSE64018 gives 1,078, but it is not independent.
- Pathways (ranked meta z): oxidative phosphorylation / electron transport DOWN (NES -2.6, padj 4e-14); TNFA/NFkB, complement, IL6-JAK-STAT3, interferon gamma, inflammatory response UP. Interferon alpha NES 1.54, padj 0.063. 9 of 68 synapse GO sets down. Driven mostly by GSE28521.

Aim 2 (WGCNA on GSE28521 cortex, 6000 variable genes, power 9, 21 modules):
- 10 modules differ by diagnosis (FDR<0.05): turquoise (1457 genes) down; green, black, pink, brown, magenta, greenyellow, red up; yellow, lightcyan down.
- Enrichment (permutation, FDR<0.05): magenta = interferon gamma (15 vs 2.6 expected), inflammatory, interferon alpha. royalblue = interferon alpha (21 of 53 genes, ~60x), IFN gamma. purple = SFARI (28 vs 11 expected for score 1-2). So H2 is supported inside GSE28521.
- BUT: the SFARI module (purple) and the interferon-alpha module (royalblue) do NOT differ by diagnosis on average (p=0.87, p=0.58). Only magenta (IFN gamma/inflammatory) is up in autism.
- Replication: 0 of 10 diagnosis-associated modules replicate. GSE102741: all module P > 0.4. Velmeshev: only 1 of 10 at p<0.05 (and in the module-level correlation r=0.07). Preservation Zsummary: turquoise 9.1, a few 2-5, magenta -1.2, green -0.9 (not preserved).
- H3 hint to test next: royalblue (IFN alpha) is not different on average, which is exactly what a subgroup effect would look like.

Cell composition (Velmeshev): protoplasmic astrocytes are ~2x more frequent in ASD samples (10.8% vs 5.2% of nuclei, FDR 0.010). No other type FDR<0.8. Compositional, region-sampling dependent; needs care.

## PIVOT (2026-10-19): drop interferon/immune framing, ask what is up/down in autism BLOOD
Brain work (scripts 03-08) is kept as-is but is no longer the main line. Interferon enrichment is not a hypothesis any more; blood DE is run unbiased, and pathways are read from the result.

GEO E-utilities search for human blood/leukocyte/PBMC ASD expression data (metadata only, nothing downloaded or checked yet; counts are the GEO "n_samples" field, which may include non-ASD groups):
| Series | N | Material | Notes |
|---|---|---|---|
| GSE25507 | 146 | PBMC, array | Alter 2011 |
| GSE42133 | 147 | leukocytes, array | Pramparo 2015, toddlers |
| GSE111175 | 141 | leukocytes, array | Pramparo/Lombardo group |
| GSE111176 | 288 | array | same group; likely overlaps GSE42133 / GSE111175 |
| GSE26415 | 84 | leukocytes, array | |
| GSE6575 | 56 | whole blood, array | Gregg 2008 |
| GSE140702 | 130 | monocytes, RNA-seq | immune-cell specific |
| GSE87847 | 93 | blood, array | ASD with ileocolitis |
| GSE77103 | 8 | PBMC, array | infants, too small |
| GSE217850 | 14 | PBMC, scRNA-seq | small |
| GSE29919 / 29918 / 29691 | 29/14/15 | LCL + blood, array | developmental delay with autistic features |
Not expression, ignore for now: GSE27044, GSE109905, GSE113967, GSE164563 (methylation); GSE89596, GSE200197, GSE222046 (miRNA/circRNA/exosome).
TO CHECK before use: autism vs control labels, platform, age/sex, and DONOR OVERLAP between the Pramparo/Lombardo series (same issue as the brain cohorts).

## Brain-only search for more frontal + cerebellum ASD data (2026-10-19)
Re-ran GEO E-utilities (autism x cerebellum/cortex/postmortem). About 130 hits, nearly all iPSC/organoid/mouse/methylation. Metadata only; nothing downloaded.
NEW human postmortem expression candidates:
- GSE269105 (13 samples, GPL28038): cerebellum, LINE-1/Alu chimeric transcripts in ASD. Check if it is plain RNA-seq with case/control, and which donors.
- GSE62098 (12, RNA-seq): corpus callosum, 6 ASD / 6 control. White matter, not frontal or cerebellum. Weak fit.
- GSE38609 (72): SuperSeries of GSE38322 (expression, 36) + GSE38608 (methylation, 36). Already have GSE38322, so no new expression samples.
- GSE222259 (30): co-expression of ASD genes, unclear whether case/control brain. Check before use.
ALREADY HAVE / EXCLUDED: GSE28521 (frontal, temporal, cerebellum), GSE102741 (DLPFC), GSE64018 (STG), GSE38322 (cerebellum, occipital), GSE211154 (Purkinje only), GSE36315 (8 samples), GSE30573 (3 v 3).
NOT postmortem brain: GSE63524 (nasal olfactory stem cells), GSE76852 (controls only), GSE113524 (Alzheimer's brains).
Conclusion: GEO has very little left. Bigger frontal + cerebellum cohorts are outside GEO (see list in chat): PsychENCODE/Synapse, and others that need verifying.

## Checked the two small GEO leads (2026-10-19, series-matrix metadata only)
- GSE222259: NOT brain. Bulk RNA-seq of iPSC-derived induced neurons, SCN2A Het vs WT (one cell line, GM8330). Excluded.
- GSE269105: postmortem CEREBELLUM, RNA-seq (GPL28038), Harvard Brain Tissue Resource Center, 13 samples (IDs like B4925). Labels in the sample table: ASD 7, unaffected 5, "ADD" 1 (B6994). The series text says ASD n=8, so the "ADD" sample is probably counted as ASD in the paper but is labelled differently in GEO. TO RESOLVE: what ADD means, and drop it or check the paper. Processed counts + FPKM in GSE269105_Processed_data_gene_expression_read_and_fpkm.txt.gz (not downloaded yet). No age/sex/RIN in the sample table. Donor overlap with GSE28521 / GSE38322 not checked (those use Autism Tissue Program IDs, this uses Harvard B-numbers).

## Multi-omic plan (2026-10-19): expression + methylation + "mutated" genes, brain only
Brain methylation series checked from series-matrix metadata (labels, region, N). Nothing downloaded yet.
| Series | Region | Autism / Control | Platform | Notes |
|---|---|---|---|---|
| GSE53162 | PFC 6/5, temporal 6/10, cerebellum 7/6 | 19 / 21 | 450K | Best multi-region set; has age, sex |
| GSE53924 | BA10 frontal + BA24 cingulate | 22 / 23 (45 samples) | 450K | Several samples per donor likely (25 donors) |
| GSE109875 | BA9 DLPFC | 10 / 6 | WGBS | Small, has age/sex |
| GSE107012 | Cerebellum | 17 / 19 | 5hmC (hydroxymethylation, not 5mC) | Different mark |
| GSE278285 | Cerebellum | 30 / 45 | 450K | Largest cerebellum methylation set |
| GSE38608 | Cerebellum 11/9, occipital 8/8 | 19 / 17 | 27K | SAME donors as expression set GSE38322 |
| GSE298243 | BA9 | 26 / 23 donors | single-nucleus methyl | Same donors as the snRNA counts |
Excluded: GSE242427 (dorsal raphe), GSE131706 (subventricular zone), GSE140730 (cord blood).
MUTATIONS: none of the brain sets above has genotypes. Options: (a) overlap DE/DMR genes with SFARI / de novo gene lists (already have SFARI), (b) call variants from RNA-seq raw reads (SRA), weak and only for expressed genes.

## Step 1 of target plan: cross-cohort UP list (script 09, 2026-10-19)
Independent cohorts: GSE28521 cortex, GSE102741, Velmeshev pseudobulk (5,884 genes measured in all 3). GSE64018 shown, not counted.
- Equal-weight Stouffer z across 3 cohorts; "up in all 3" = 1,142 genes (19.4%, vs 12.5% if pure noise). 58 genes up in all 3 with FDR<0.05 (8 in SFARI); 208 with P_meta<0.01 (25 SFARI). Top: CPNE3, ELK1, RAB23, PLOD1, NR2F1, KLF15, ZFP36L1, IGFBP5.
- CAVEAT: GSE28521 has inflation lambda 3.6 (Velmeshev 1.5). After deflating each cohort's z by sqrt(lambda) (genomic control), 0 genes reach FDR<0.05 and 0 up genes reach FDR<0.10. Top-25 genes get ~60% of their meta z from GSE28521 alone. Median mean-logFC of all-3-up genes is only 0.085.
- Conclusion: the list is a RANKED set of hypotheses, not a validated set. Files: results/09_cross_cohort_all_genes.csv, results/09_up_candidates.csv.
- "Known vs novel" is SFARI only so far; the literature cross-check (PubMed) is still to do.

## Step 2: supporting evidence + druggability for the 208 UP candidates (scripts 10, 11; 2026-10-19)
Candidate set = up in all 3 independent cohorts AND meta P<0.01 (208 genes; results/09_up_candidates.csv).
Support (script 10, results/10_candidates_with_support.csv):
- Cerebellum (GSE28521 cerebellum, 10 autism / 11 control, 2 SVs, lambda 1.77, 0 genes FDR<0.05): 130 of 192 candidates also up (96 expected by chance); 23 with P<0.05. SAME COHORT as the cortex data, so not independent replication.
- WGCNA: 166 of 208 sit in a module that is up in autism (FDR<0.05). CIRCULAR: the candidates are driven by GSE28521, and those modules were built and tested in GSE28521. Do not count as separate evidence.
Druggability (script 11, UniProt + ChEMBL queried 2026-10-19, results/11_candidates_druggability.csv):
- 208/208 have a reviewed UniProt entry and an AlphaFold model; 132 have >=1 PDB structure; 91 have a ChEMBL single-protein target; 14 have a drug/clinical mechanism; 21 have >=100 potent (pChEMBL>=6) compounds.
- Max clinical phase 4 (approved drug with recorded mechanism): RAF1, CSNK1E, LYN, F3, XPO1, ACVR2A, ABCA1. Phase 1-3: NQO1, MERTK, CXCL12, GPNMB, NCAM1, FYN, TNFRSF1A, MKNK2(?) see csv. Other tractable kinases/enzymes with many potent compounds but no mechanism entry: MKNK1/2, DDR1, DYRK2, SPR, AKR1C3, MAT2A, BRPF1, HSP90B1.
- Transcription factors among top hits (ELK1, NR2F1, CREB1, CXXC4) have no ChEMBL compounds: poor docking targets.
- CAVEAT: at least one ChEMBL call timed out (curl status 28), so some "no ChEMBL target" (117) may be failures, not true absences (e.g. KLF15, ZFP36L1, PHF21A show NA). Re-query NA rows before concluding.
- "novel" is still vs SFARI only; PubMed check pending.

## Step 3: ChEMBL re-query, drug names, PubMed novelty (script 12, 2026-10-22)
Files: results/12_target_drug_pairs.csv (25 records, 24 molecules, 14 genes), results/12_candidates_full.csv (208 genes).
- Re-query: the 117 genes with no ChEMBL target are TRUE absences (0 found, 0 failures after retries). So 91 targets, 14 with a recorded drug mechanism. 27 genes are "druggable" (mechanism or >=100 potent compounds).
- PubMed (title/abstract, gene AND autism terms, queried 2026-10-22): median 0 papers; 108 of 208 genes have 0, 67 have 1-4, 33 have >=5. Search is by gene SYMBOL, so short or ambiguous symbols (F3, SPR, VHL, ABCA1) give noisy counts; 0 papers does not prove no autism link (could be in full text only).
- Drug-target pairs of note: RAF1 (sorafenib, regorafenib), LYN (bosutinib), CSNK1E (umbralisib), XPO1 (selinexor), ABCA1 (probucol), NQO1 (vatiquinone, MeSH indications include Tourette and pervasive child development disorders; VERIFY what this trial was), FYN (JNJ-26483327). Most are cancer drugs: repurposing into brain would need blood-brain-barrier and safety checks. My "brain_indication" flag is a loose regex and gives false positives (e.g. "neuro" in neoplasm names): do not trust it without reading indications.
- DOCKING SHORTLIST (druggable, up in all 3 cohorts, few autism papers): MKNK2 (cereb up, 1 paper, 15 PDB), DDR1 (cereb P=0.014, 1 paper), DYRK2 (0 papers, 26 PDB), BRPF1 (0 papers, cereb P=0.046), LYN (cereb P=0.012, approved inhibitor bosutinib), CSNK1E (SFARI, approved inhibitor umbralisib: use as KNOWN-GENE control), HSP90B1, AKR1C3, MAT2A. RAF1 and FYN are heavily studied but have many ligands (good docking controls).
- Honest status: nothing here is validated. Candidates come from weak, partly inflated expression evidence.

## Docking tools installed (2026-10-22)
Miniforge (conda 26.7.2) at ~/miniforge3, Apple M1 (arm64). Environment "docking" (python 3.11): RDKit 2026.03.1, AutoDock Vina 1.2.7, Meeko 0.8.0, Open Babel 3.2.1, Biopython 1.88, pandas, numpy, scipy, requests.
Use without activating:  ~/miniforge3/envs/docking/bin/python  (also vina, obabel, mk_prepare_ligand.py, mk_prepare_receptor.py in the same bin folder).
Disk after install: 37 GB free (env + miniforge about 1.4 GB). Nothing was installed system-wide; R and the system Python are untouched.
Next: pick a PDB structure with a bound ligand for each shortlisted target and run the re-docking control before any scoring.

## Docking step 1: MKNK2 re-docking control (script 13, 2026-10-22)
Method: crystal ligand removed, re-docked from a random RDKit conformer with Vina 1.2.7 (exhaustiveness 32, seed 42, 24 A box on the crystal ligand centre, rigid receptor from Open Babel, chain A, zinc/waters removed). Pass = top-scored pose < 2 A heavy-atom RMSD (symmetry aware) to the crystal pose.
| PDB | Ligand | Res. | Protein | Top score (kcal/mol) | Top-pose RMSD | Best of 9 | Verdict |
|---|---|---|---|---|---|---|---|
| 8P9B | X8K (31 heavy atoms) | 2.59 A | MNK2 D228G | -11.23 | 0.70 A | 0.70 A | PASS |
| 6JLR | BV9 (21 heavy atoms) | 2.90 A | MNK2 wild type | -9.42 | 9.30 A | 1.60 A | PARTIAL (top pose wrong; 6th pose right) |
- Use 8P9B as the MKNK2 docking receptor. 6JLR (small ligand, many near-equal scores) shows Vina cannot reliably rank poses of small flat ligands here, so single-pose conclusions for small compounds are weak.
- 8P9B is the D228G mutant (a point change in the activation loop used for crystallisation), not the wild-type sequence. Note this when reporting.
- Output: results/13_redock_controls.csv, docking/<PDB>_<LIG>/ (poses, prepared receptor). Structures cached in data_raw/structures/.
- A pass shows the setup recovers a known pose; it does NOT show Vina scores rank different drugs correctly. That needs the decoy test (known inhibitors vs property-matched decoys).

## Docking step 2: panel vs MKNK2 (scripts 14, 15; 2026-10-23)
Receptor 8P9B (MNK2 D228G, validated by re-docking, 0.70 A). Panel = 291 ligands from ChEMBL (queried 2026-10-22): 50 potent known MNK inhibitors (pChEMBL>=7.5, max 3 per scaffold), 200 property-matched decoys (approved non-kinase drugs, 4 per active; approximation of DUD-E decoys), 23 brain-disorder drugs, 18 kinase/target drugs. Lithium dropped (metal ion, cannot be docked); EGCG and leucettine L41 not found by exact name in ChEMBL. Vina 1.2.7, exhaustiveness 16, seed 42, rigid receptor, ligands from RDKit default protonation (no pH assignment). All 291 docked, 0 failures.
- DECOY TEST: actives separate from decoys. AUC 0.86 (raw Vina score; mean -9.53 vs -7.76 kcal/mol, p=2.8e-15); AUC 0.79 on ligand efficiency (score per heavy atom), so it is not just a size effect (actives 33.2 vs decoys 32.8 heavy atoms on average). Top 10% by score (29 ligands): 19 actives, 5 decoys (base rate of actives 17%).
- Positive control tinodasertib (eFT508) scores best of everything: -11.42.
- Brain-disorder drugs, raw score: lurasidone -10.05, risperidone -9.92, paliperidone -8.90, buspirone -8.70, ziprasidone -8.57, donepezil -8.55 ... aripiprazole -8.20, quetiapine -6.42, memantine -4.75, valproic acid -5.10.
- Known kinase drugs: umbralisib -10.00, regorafenib -9.77, sorafenib -9.32, imatinib -9.24, dasatinib -8.07, bosutinib -7.51 (an approved kinase inhibitor scoring below the decoy mean: a warning about scoring accuracy). Selinexor -9.87 binds XPO1 covalently: its MKNK2 score is not meaningful.
- WHAT THIS DOES NOT SHOW: a docking score is not binding. No pose inspection (hinge hydrogen bond) done. Rigid receptor from ONE crystal conformation. Lurasidone/risperidone are D2/5-HT antagonists, not known MNK inhibitors; need a ChEMBL activity check and pose check before any claim. MKNK2 is also a weak autism candidate (P_meta 1.6e-4 inflated by GSE28521; cerebellum P=0.77).
- Output: results/15_MKNK2_panel_scores.csv; log notes/15_mknk2_dock.log.

## Docking step 3: do the top brain drugs look like real MNK2 binders? (script 16, 2026-10-23)
(a) ChEMBL MEASURED activity on MNK2 (CHEMBL4204) / MNK1 (CHEMBL4718): lurasidone, risperidone, paliperidone, ziprasidone, buspirone, donepezil = 0 records on both. 0 records means never tested/curated, NOT proven inactive. Controls have data: tinodasertib IC50 86 nM (MNK2), 10.6-64 nM (MNK1); staurosporine Kd 17-22 nM (MNK2).
(b) Poses (8P9B, 9 poses/ligand, seed 42). Hinge taken from the crystal ligand = Met162 backbone N (hinge segment Val158-Ser166 checked in a second pass). Test = ligand N/O within 3.4 A of hinge backbone N/O, pose centroid <5 A of crystal ligand.
- Tinodasertib: top pose makes the hinge bonds (Glu160 O, Met162 N). Only ligand whose TOP pose does.
- All six brain drugs: top pose makes NO hinge H-bond. Hinge-contacting poses exist only at lower ranks: risperidone rank 5 (Glu160 O + Met162 N, score -9.3 vs top -9.44), buspirone rank 5, donepezil rank 2; lurasidone, paliperidone, ziprasidone only touch other segments (Ser166 N / Glu160 O).
- Top decoy eltrombopag scores -10.12 (similar to lurasidone) with no proper hinge binding: high Vina scores here can come from generic hydrophobic filling of the pocket.
- Limit of the test: staurosporine (real potent binder) shows no hinge H-bond in any in-site pose, so Vina/this rigid receptor misplaces at least one true binder. The hinge test can miss real binders; it cannot prove binding.
CONCLUSION: no support for lurasidone/risperidone/other brain drugs as MNK2 inhibitors (no measured data, top poses lack the hinge bond, a decoy scores as high). Do NOT present them as candidates. Risperidone's rank-5 hinge pose is the only weak lead and is not worth following without wet-lab data.
Also found: decoy set contains duplicates of one drug in salt forms (eltrombopag free base + choline + olamine salts); strip duplicates by parent compound in a rerun.
Files: results/16_pose_check_all_poses.csv, results/16_chembl_activity_check.csv, docking/8P9B_X8K/poses_*.pdbqt.

## Docking controls for LYN and DDR1 (script 13, 2026-10-23)
- DDR1, 4BKJ chain A, ligand STI (imatinib), 1.7 A, box 30 A, exh 32, seed 42: top pose 0.80 A -> PASS. Receptor for DDR1 panel = docking/4BKJ_STI.
- LYN, 5XY1 chain A, ligand 8H0 (30 heavy atoms), 2.7 A (only LYN kinase-domain structure with a ligand; 8WFF is the SH3 domain, 6NMW/1W1F/1WA7 SH3): seed 42 exh 32: top 7.86 A. Retried at exh 64 with seeds 1/2/3: top-pose RMSD 2.43 / 7.88 / 2.38 A; best-of-9 1.52 / 1.94 / 1.57 A. Never below 2 A at rank 1 -> NOT validated. Vina finds the right pose but cannot rank it first. LYN panel scores would not be trustworthy; LYN is on hold (could be revisited with a different structure/flexible residues, or a different scoring function).
- script 13 now takes SEED and EXH as environment variables and an optional box size (4th argument).
- DDR1 panel (docking/panel_DDR1.csv, 189 ligands: 30 DDR1 actives pChEMBL>=7.5 [label still reads 1_MNK_inhibitor in the file], 120 decoys deduplicated by parent compound, 23 brain drugs, 16 kinase/target drugs) docking at exhaustiveness 8, box 30; checkpointed to results/15_DDR1_partial.csv.

## Docking step 4: DDR1 panel + checks (scripts 14-16, 2026-10-23)
Receptor 4BKJ chain A (DDR1 kinase domain + imatinib, DFG-out, 1.7 A; re-docking top pose 0.80 A). Panel docking/panel_DDR1.csv: 30 potent DDR1 actives (pChEMBL>=7.5), 120 decoys (deduplicated by parent), 23 brain drugs, 16 kinase/target drugs; Vina exhaustiveness 8, box 30 A; 188/189 docked (buprenorphine failed 3D embedding).
- Enrichment: AUC 0.89 raw score (actives -11.53 vs decoys -9.19, p=2.7e-11), 0.76 on ligand efficiency; top 10% (18 ligands) = 12 actives, 0 decoys.
- Top of ranking = known DDR1 binders: ponatinib -13.79, nilotinib -13.73, regorafenib -13.58, BAY-826 -13.20, sorafenib -13.14, imatinib -12.51. ChEMBL measured DDR1 activity: ponatinib IC50 9.4 nM, imatinib IC50 43 nM / Kd 0.7 nM, dasatinib Kd 0.69 nM.
- Receptor bias: 4BKJ is DFG-out, so TYPE II inhibitors score high and TYPE I are penalised: dasatinib (real Kd 0.69 nM) scored only -8.89, staurosporine -9.64, bosutinib -8.90. Scores reflect fit to ONE rigid conformation.
- Decoy contamination: 3 of the top-5 "decoys" are real kinase inhibitors that slipped past the name filter (netarsudil and belumosudil = ROCK; tivozanib = VEGFR) so true enrichment is slightly under-estimated. Next rebuild: filter decoys by ChEMBL kinase activity, not by name suffix.
- Brain drugs (raw score, rank of 188): aripiprazole -11.01 (38), risperidone -10.50 (50), lurasidone -10.34 (54), paliperidone -10.25 (56), ziprasidone -9.71 (79). Mean of brain drugs -8.53 < decoy mean -9.19 (they are smaller). ChEMBL measured DDR1 activity: 0 records for all of them.
- Poses (script 16, tagged _DDR1): the only backbone contacts found from the crystal ligand were Asp784 N (the DFG Asp, NOT the Met704 hinge), Met704 N, Val763 O, His764 O; hits were all Asp784 N. Aripiprazole's top pose sits in the site and touches Asp784 N, as do ponatinib, pranlukast (decoy) and netarsudil (decoy). Risperidone, lurasidone, paliperidone, ziprasidone top poses lie 6-8 A from the ATP site (outside the pocket).
CONCLUSION: the docking setup discriminates real kinase inhibitors from decoys (MKNK2 AUC 0.86, DDR1 AUC 0.89) and re-finds known inhibitors, but gives NO credible support for repurposing the brain-disorder drugs onto MKNK2 or DDR1 (no measured activity, mostly off-site or non-specific poses, decoys of similar score). Aripiprazole/DDR1 is the only in-site pose and still has no data behind it.
What the docking CAN support: if DDR1 is a real autism target, approved/clinical type-II inhibitors with measured DDR1 potency (ponatinib, nilotinib, imatinib, regorafenib, sorafenib) are the repurposing leads. All are oncology drugs; brain penetration and safety are unaddressed.

## Is DDR1 a credible autism target? Three checks (scripts 17, 18; PubMed; 2026-10-23)
### 1. Independent support (script 17; results/17_independent_support.csv, 17_methylation_probe_map.csv; log notes/17_independent_support.log)
Expression, autism minus control, log2 scale:
- GSE269105 cerebellum RNA-seq (7 ASD v 5 control, "ADD" sample dropped): DDR1 -0.08, p=0.83 (NULL). MKNK2 +0.21, p=0.67.
- GSE38322 cerebellum (14 v 12): DDR1 +0.28, p=4e-4 (t=4.22); MKNK2 +0.43, p=0.0014. BUT (a) this cohort shares 7 donors with GSE28521, (b) DDR1's probe sits at rank ~4,780 of 47,323 probes, so large t values are common in this dataset (broad shift), not specific to DDR1. Occipital (4 v 6): DDR1 +0.23, p=0.30 (not significant).
Methylation (450K; DDR1 has 89 probes, 58 promoter-class; effects are mean delta-beta, autism minus control):
- Cerebellum promoter: GSE278285 (30 v 45, adjusted for age+sex) +0.013 (33/58 probes p<0.05, min p 2.5e-7); GSE53162 (7 v 6) +0.0075. Same sign in both, but ~1% methylation difference.
- Frontal promoter: GSE53162 -0.007 (6 v 5) vs GSE53924 BA10 +0.023 (12 v 12): OPPOSITE signs. Cingulate (GSE53924) -0.025. Temporal -0.003.
- Probes within a gene are correlated, so "many probes p<0.05" is not many independent hits. No cell-type composition correction. Direction note: promoter hypermethylation would normally go with LOWER expression, the opposite of the upregulation we saw.
### 2. Literature (PubMed title/abstract, queried 2026-10-23)
- DDR1 total 846 papers; with autism 1; brain/neuro 63; myelin/oligodendrocyte 20; schizophrenia/bipolar 13; cerebellum 2.
- The one autism paper (PMID 35215271, 70 Arab ASD patients, blood RNA-seq + exome genotyping, 2022): a regulatory variant (rs1940475 etc.) linked to DOWNregulation of DDR1. Opposite direction to our brain finding and a different tissue.
- DDR1 is NOT novel in neuropsychiatry: DDR1 hypermethylation in psychosis (leukocytes and DLPFC, PMID 31920096), DDR1 methylation in bipolar disorder, DDR1 coexpression with oligodendrocyte/schizophrenia-risk genes, DDR1 in oligodendrocyte differentiation and myelin repair, a 2025 review "DDR1 as a therapeutic target for neurological and psychiatric disorders" (PMID 41045082), and a 2025 nilotinib (a DDR1 inhibitor) trial in dementia with Lewy bodies (verify). Novelty for autism specifically is plausible (1 paper), novelty for psychiatry in general is not. A full literature search (full text, preprints, other aliases) has not been done.
### 3. Mutation layer (script 18; denovo-db open files, SSC + non-SSC; results/18_denovo_gene_counts.csv)
- 5,099 autism and 1,898 control coding de novo variants (LoF + missense) after de-duplication. Sanity check passes: CHD8, ARID1B, SCN2A, ADNP, SYNGAP1, ASH1L, DYRK1A, SHANK3 lead the LoF excess.
- DDR1: 0 LoF and 0 missense de novo variants in autism probands (expected 0.06 and 0.32 by coding length) and 0 in controls. MKNK2: 0 LoF, 0 autism missense (1 control missense).
- Across the 208 candidates (199 present in denovo-db): 11 genes with an autism LoF, 41 with an autism missense; only CSNK1E and NR2F1 (2 missense each, Poisson p~0.01, uncorrected over ~400 tests so expect ~4 by chance) reach nominal excess. No candidate has a convincing mutation signal. Length-only null; ignores mutability.
VERDICT: DDR1 is a hypothesis, not a validated autism target. Evidence: up in 3 discovery cohorts (inflated), one partly overlapping replication (GSE38322), a null in the one independent cerebellum RNA-seq (GSE269105), small and inconsistent methylation shifts, opposite direction in the single published autism DDR1 result, no coding mutations. Its real interest is the link to the psychiatric DDR1 literature (myelination, hypermethylation) and existing DDR1-inhibitor drugs.

## Option 1 + 2 follow-up (scripts 19-21; 2026-10-24)
### Option 2: myelin / oligodendrocyte gene sets in autism brain (scripts 19, 20; results/20_myelin_pathway.csv; log notes/20_myelin_pathway.log)
Method: per cohort/region, gene-wise t (autism vs control, cohort covariates), set statistic = mean t, tested by permuting sample labels (3000x; keeps gene-gene correlation). A second "vs random size-matched sets" p is printed but is anti-conservative (genes in a GO set are correlated; random sets are not), so ONLY the label-permutation p is used for conclusions. Cohorts: cerebellum GSE28521 (10 v 11), GSE38322 (14 v 12), GSE269105 (7 v 5); cortex GSE28521 frontal/temporal, GSE102741 DLPFC, GSE64018 STG; Velmeshev oligodendrocyte-only pseudobulk (PFC 13 v 10, ACC 9 v 9).
- Myelin sheath (GOCC_MYELIN_SHEATH): cerebellum GSE28521 t -0.02 (p 0.97); GSE269105 -0.16 (p 0.70); GSE38322 +1.32 (p 0.001). Cortex: all |t|<0.6, no p<0.2 (DLPFC -0.27, Velmeshev oligodendrocytes PFC -0.28 / ACC +0.05). Same pattern for oligodendrocyte differentiation/development, ensheathment of neurons, regulation of myelination, 15 hand-picked oligodendrocyte markers.
- The only "hit" is GSE38322 cerebellum, and it is NOT specific: the same cohort moves oxidative phosphorylation by t -2.5 (p 3e-4), astrocyte differentiation +1.0 (p 3e-4), inflammatory response +0.33 (p 0.007), i.e. a broad shift (as already seen for DDR1's probe rank). It shares 7 donors with GSE28521, where the myelin signal is absent, and the one fully independent cerebellum RNA-seq (GSE269105) points the other way (small, non-significant).
- Inside oligodendrocytes (Velmeshev snRNA), myelin genes are not changed (permutation p 0.27 PFC, 0.80 ACC). Oligodendrocyte fraction is lower in ASD (10.8% v 17.4%) but not significant (FDR 0.81, script 08), so composition does not carry it either.
- Context that DOES agree with earlier work: oxidative phosphorylation is lower in autism in 7 of 8 bulk cohort-regions (permutation p<0.05 in GSE28521 frontal, GSE38322 cerebellum/occipital, GSE64018; GSE102741 is flat at t=+0.01). Immune/inflammatory sets are mildly up in most.
CONCLUSION (option 2): no support for a myelination/oligodendrocyte change in autism cerebellum or cortex in these open datasets. The DDR1 -> myelin story does not carry over to expression data. Caveat: n is 12-26 per cohort, microarray + postmortem confounds, so this rules out only a large effect.

### Option 1: DDR1 checks on all 208 candidates (script 21; results/21_candidate_evidence_matrix.csv; log notes/21_candidate_evidence_matrix.log)
Flags declared before looking: E1 GSE269105 cerebellum t>0 and p<0.05; E2 GSE38322 cerebellum robust-z>1.64 (t standardised against all genes to remove the broad shift); M promoter methylation lower than the cohort's own genome-wide centre (robust z<-1.64 in GSE278285 and same direction in GSE53162 cerebellum); D de novo autism LoF, or missense Poisson p<0.01. 37 of 208 are "druggable" (mechanism drug in ChEMBL, or >=10 potent compounds).
- Candidate list vs all other genes in independent GSE269105 cerebellum: NOT enriched (median t +0.05, p=0.19; 42.6% v 40.7% positive). In GSE38322 (donor overlap, circular) and GSE28521 cerebellum the list is strongly up (83% and 73% positive v ~50% and 45%). So the discovery signal does not reappear in the one independent cohort.
- Pass counts (expected by chance in brackets): E1 0 [~5], E2 53 [~10, but donor overlap], M 0 [~5], D 13. Layers passed: 0 -> 148 genes, 1 -> 54, 2 -> 6, none with 3. E1 and M, the two truly independent layers, are empty.
- CORRECTION to the DDR1 methylation entry above: GSE278285 has a near-global shift (median promoter dBeta +0.0086; 70% of all genes higher in autism; 64% of ALL genes p<0.05) and GSE53924 frontal +0.020 (99% of genes up). The DDR1 "hypermethylation" (+1.3%) is that shift. Centred on the cohort's genome-wide median DDR1 is +0.004 (robust z -0.11), i.e. nothing. DDR1 methylation is therefore not evidence for or against it. (Likely cell-type composition or batch; not resolved.)
- The 2-layer genes are E2 + de novo LoF: PHF21A, PTPN13, NFIA, PBX1, TNRC6B, NLGN4X. 5 of 6 are SFARI genes (PTPN13 is not). None is druggable. 25 of the 208 candidates are SFARI genes; 10 of the 13 D genes are SFARI genes, so D mostly rediscovers known autism genes (which is a sanity check, not a new finding).
- Druggable candidates: 9 pass one layer, none two. 8 of them pass only E2 (MKNK2, ABCA1, BRPF1, INPPL1, DDR1, RAF1, YAP1, BCL6; the layer with donor overlap). Only CSNK1E passes D (2 missense, p=0.0095; SFARI-listed), but in GSE269105 CSNK1E is DOWN (t -2.36, p=0.04). ChEMBL lists umbralisib for CK1e (its main target is thought to be PI3K-delta; verify before using).
- DDR1 has no special standing: same evidence profile as MKNK2, RAF1, YAP1 (E2 only). 
CONCLUSION (option 1): none of the 37 druggable candidates has more than one layer of weak support, and the two independent layers are empty for all 208. No druggable candidate is better supported than DDR1.

### Where this leaves the project
- Defensible claims: (1) pipeline works end to end (DE -> WGCNA -> candidate list -> docking validated by decoy test and re-docking); (2) the up-regulated candidate list does not replicate in the independent cerebellum RNA-seq or in methylation after correcting for the global shift; (3) mutation data recover known autism genes but none of the druggable candidates; (4) myelin/oligodendrocyte pathway not supported.
- Not defensible: any named drug-target pair as an autism candidate.
- Remaining legitimate routes: report as methods + negative result; or flip the candidate source to genetically supported genes (SFARI/de novo, e.g. CSNK1E, NR2F1, PBX1, NFIA) and look for druggability there; GSE298243 snRNA (not yet used) is one more independent expression check.
