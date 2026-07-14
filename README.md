# Strict-cis proteogenomic atlas across cardiometabolic and neuropsychiatric outcomes

This repository contains the minimal code and derived inputs needed to inspect
the principal analysis definitions and redraw the four main figures for the
manuscript **"A strict-cis proteogenomic atlas for prioritizing circulating
proteins across cardiometabolic and neuropsychiatric outcomes."**

The repository is intentionally limited. It does not redistribute source GWAS,
QTL, LD-reference, participant-level, or institutional data, and it is not a
complete archive of every supplementary or audit analysis.

## Repository contents

- `R/01_reproduce_strict_cis_core.R`: reusable definitions for coding-gene
  +/-500-kb cis regions, instrument filtering and clumping parameters,
  exact-variant MR, pairwise colocalization, and the 16-point ranking rubric.
- `R/02_rebuild_main_figures.R`: redraws Figures 1-4 from the released derived
  input table.
- `data/main_figure_inputs.tsv`: aggregate, non-participant-level values used
  by the four main figures, including the post hoc protein-altering-variant
  sensitivity summary shown in Figure 2.
- `environment.yml`: principal R and command-line dependencies.
- `LICENSE`: MIT License for code.
- `LICENSE-DATA`: CC BY 4.0 terms for the derived figure input table.

## Data availability

The derived values required to redraw the main figures are openly available in
`data/main_figure_inputs.tsv` under CC BY 4.0. The file contains aggregate MR,
regional posterior, evidence-state, and ranking results only. It contains no
participant-level records, identifiers, credentials, local paths, or source
summary-statistic caches.

The source summary statistics remain available from their original providers
under provider-specific access and redistribution terms. The primary protein
resource and seven canonical outcome datasets were accessed through MRC IEU
OpenGWAS. The canonical outcome identifiers were:

| Outcome | OpenGWAS identifier |
|---|---|
| Major depressive disorder | `ieu-b-102` |
| Alzheimer disease | `ebi-a-GCST90027158` |
| Anxiety disorders | `finn-b-KRA_PSY_ANXIETY` |
| All-cause dementia | `finn-b-F5_DEMENTIA` |
| Cognitive performance | `ebi-a-GCST006572` |
| Heart failure | `ebi-a-GCST009541` |
| Type 2 diabetes | `ebi-a-GCST90018926` |

Additional resources included eQTLGen, GTEx, UKB-PPP, European LD-reference
data, SMR BESD resources, external disease GWAS, and tissue, CSF, or brain-cell
QTL datasets. These source files and their local caches are not included here.
Users must obtain them from the cited providers and comply with the applicable
terms. An OpenGWAS JWT, when required, must be supplied through the
`OPENGWAS_JWT` environment variable; credentials must never be committed.

## Installation

Create the base environment:

```bash
conda env create -f environment.yml
conda activate strict-cis-proteogenomic-atlas
```

The analysis used R 4.5.1, `data.table` 1.18.2.1, `TwoSampleMR` 0.7.0,
`ieugwasr` 1.1.0.9000, `coloc` 5.2.3, `hyprcoloc` 0.0.2, and PLINK
v1.90b7.2. Installation notes for packages not consistently available through
conda are included in `environment.yml`.

## Method self-check

The core script can be checked without downloading source data:

```bash
Rscript R/01_reproduce_strict_cis_core.R --self-check --print-outcomes
```

For local LD clumping, set `PLINK_BIN` to the PLINK executable and `LD_BFILE`
to the prefix of a permitted PLINK-format reference panel. The repository does
not provide an LD panel.

## Redraw the main figures

```bash
Rscript R/02_rebuild_main_figures.R
```

The script creates `Figure1` through `Figure4` as vector PDF and 600-dpi PNG
files in `output/`. A different destination can be supplied with the
`OUTPUT_DIR` environment variable.

## Interpretation caveats

- This is a summary-data target-prioritization analysis, not evidence of drug
  efficacy, safety, neuroprotection, or a molecular mechanism.
- The cardiometabolic and neuropsychiatric labels are outcome-defined strata;
  the analyses did not support biological separation between them.
- Cardiometabolic-outcome proteins should not be interpreted as
  neuroprotective targets.
- The post hoc protein-altering-variant sensitivity did not update discovery
  classes, scores, or tiers. A retained association does not exclude other
  assay-binding or pleiotropic effects.
- The evidence score is a heuristic ranking rubric with partly related
  components. It is not a probability of target validity or treatment success.
- Most source datasets were predominantly of European ancestry. The bounded
  non-European analyses were instrument-transportability checks, not
  cross-ancestry MR or colocalization replication.
- The released plotting table supports verification and redrawing of the main
  figures. It does not replace the complete numerical Supplementary Tables.

## License

Code is released under the MIT License. The derived figure input table is
released under CC BY 4.0. Third-party source data remain governed by their
original terms and are not relicensed by this repository.

## Citation

Citation details will be added after journal publication.

> Study 10 authors. A strict-cis proteogenomic atlas for prioritizing
> circulating proteins across cardiometabolic and neuropsychiatric outcomes.
> Manuscript under review.
