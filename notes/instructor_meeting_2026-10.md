# Instructor meeting: what was said and how it applies here (notes written 2026-10-19)

Source: pasted transcript of the call. Gut microbiome / metagenomics part deliberately left out (user asked to ignore it).
These are the instructor's views as spoken, not verified facts. Marked [CLAIM] where it is a biological or technical statement I have not checked.

## Why the meeting happened
Machine-learning attempts and the expression analysis gave null results. Instructor: do NOT go back to ML (chasing a hypothesis that may not exist). Instead extract what is in the data with standard workflows, then find novelty by removing what is already known in the literature.

## The directions the instructor proposed
1. EXPRESSION -> biomarkers. Differential expression gives highly up/down genes. Cross-check each with literature. Genes not known in autism = potential novel therapeutic target. "The simplest" route.
2. GENOMICS -> mutations. Needs whole-exome or whole-genome sequencing data (not RNA-seq). Goal: driver genes = genes carrying the most mutations in autism. Instructor notes autism is genetic, so reviewers may say expression targets only suppress symptoms, not cure.
3. EPIGENOMICS -> DMRs. Hyper/hypo methylated regions replace up/down. A DMR is a chromosome region that could be a target.
4. CHEMINFORMATICS, layered on direction 1.

## Cheminformatics, as described
- Input is ONLY the NAME of the highly upregulated gene(s) from differential expression.
- Get that gene's protein structure from a free database (no institutional ID needed).
- Take known drugs, especially ones already used for brain disorders (schizophrenia, bipolar, similar). Ask: what do these drugs target, how do they work, and how is that target related to our differentially expressed genes.
- Molecular docking: bind the protein and the drug; a significant binding score suggests the drug could act as an INHIBITOR of that target = candidate drug.
- Output = a (new) therapeutic target + a candidate drug. Call them "candidate", because nothing is tested in a wet lab.
- Example given: paracetamol / COX-2 (a known drug with a known target).
- Instructor caveat [CLAIM]: for mutation-driven disease you cannot fix the gene with a binding drug; drugs would suppress symptoms. CRISPR was discussed and judged not a good idea here. Infectious-disease logic (block a hijacked host protein) does not carry over.
- The biomarker half has a second use: a highly changed gene can be a DIAGNOSTIC marker (BRCA1/2 example), or a cluster of co-expressed genes can be a marker panel.

## Workflow the instructor wants shown (methods layers)
differential expression -> pathway analysis -> network analysis -> co-expression (WGCNA) -> (cheminformatics docking). More layers = more complexity and more to present.
Novelty claim: after cross-checking the literature for each approach, whatever is left is the novel part.

## Practical constraints stated
- Timeline: about two months.
- Data access: no institutional ID, so controlled-access resources are OUT. The user's data is public GEO (and similar open sources like SRA, TCGA-type open data). "Don't waste time" on institutional-login sites (this includes the Synapse/PsychENCODE route I suggested earlier).
- Compute: raw whole-genome / shotgun data needs heavy compute (cloud or server). Check data size before committing.
- Instructor advice: use BRAIN data, not blood. Narrow to the brain regions most related to autism (cortex etc.).
- Decide first: (A) find a treatment/target or (B) explain disease mechanism (systems biology). Data and tools follow from that choice.
- Homework: read autism omics papers from the last year to see which omics and questions are current. Instructor to send a written step list per idea.

## How this maps onto what we have (drsefProjectv3)

| Instructor step | Status here |
|---|---|
| Differential expression | Done (scripts 04-05), but results are mostly null/weak; only GSE28521 has many hits and it is inflated (lambda 3.6) |
| Pathway analysis | Done (fgsea, 05) |
| WGCNA / co-expression | Done on GSE28521 (07); 10 modules differ by diagnosis, 0 replicate |
| Network analysis | NOT done (e.g. PPI network, hub genes of modules) |
| DMR / methylation | Datasets identified (notes/week0_findings.md), not downloaded |
| Mutation data | Instructor says WES/WGS. None of our brain data has it. See below |
| Literature cross-check for novelty | NOT done; needs a curated "known autism genes" list beyond SFARI |
| Cheminformatics | Not started; no tools installed |

Points where the instructor's guidance changes my earlier suggestions:
1. PsychENCODE/Synapse: drop it. It needs the institutional login. Earlier I recommended it; that does not fit.
2. Mutations: calling variants from RNA-seq (my option b) is a weak substitute for the WES/WGS the instructor means. For direction 2, search GEO/SRA/dbGaP-open-tier for public autism exome/genome VCFs or summary tables. Most autism WES cohorts are controlled access, so this may be a dead end. The realistic open route is published gene lists from exome papers (supplementary tables), used as an overlap check.
3. Cheminformatics only needs gene names + free structure databases (UniProt, PDB, AlphaFold DB) + free drug databases (ChEMBL, DrugBank open data, DGIdb). No login problem there. Local tools still must be installed.
4. Blood: instructor says brain, which matches the move back to brain-only.

## Fit with our data: the main risk
Instructor's step "target the highly upregulated genes" assumes a trustworthy list of them. Our current data gives 3 meta-analysis genes with FDR<0.05 and none replicated. Docking against genes that do not replicate would produce a candidate list built on weak evidence. Ways to protect the result:
- choose targets by agreement across cohorts and regions, not by fold change in one cohort
- use co-expression modules / pathways (as the instructor suggested) so the target is a group, not one noisy gene
- state clearly that candidates are in silico only
- novelty check: remove genes already linked to autism in the literature

## Open questions for the instructor / user
1. Direction: treatment/target (A) or mechanism (B)?
2. Is mutation direction wanted, given likely no open WES/WGS? Overlap-with-published-gene-lists is the open alternative.
3. Which last-year autism omics papers to read first (the user was asked to do this).
