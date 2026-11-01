# Pre-registration for an out-of-sample test of the RFO result

Written **before any expression data from the candidate studies were quantified or seen**.
Only SRA/ENA record metadata (accessions, library names, BioSample structure, release
dates, cultivar and genotype fields) had been inspected at the time of writing. The point
of fixing this now is that the discovery cohort's panel, background and test were all
settled after results had been seen; nothing here may be revised once counts exist.

## Candidate studies

Both were present in the original 62-BioProject screen and were excluded by the automated
labeller, for reasons that are recoverable rather than disqualifying.

| dataset | why it was not in the discovery cohort | what is actually there |
|---|---|---|
| PRJNA1380757 | only 6 experiments were public when screened | 12 libraries, 12 BioSamples, released 2025-12-17 |
| PRJNA931112 | all 12 libraries scored "unlabeled": the labeller does not read `Ct`/`Dr` | 12 libraries, 2 BioSamples, released 2025-01-02 |
| Zenodo 15040335 | never in the screen: it is not in SRA, GEO or ArrayExpress | 48 samples, 4 genotypes x 2 treatments x leaf and root x 3 replicates |

None shares a BioProject, BioSample or experiment accession with the 13 discovery strata.

### The third candidate is admitted on a different footing

Zenodo 15040335 (Juan-Cabot et al., 2025) deposits **count data only** — no reads — so it
cannot be re-quantified with this project's salmon/ITAG4.1 pipeline, and its matrix is keyed
on Entrez GeneIDs produced by the authors' own annotation. That breaks the uniform
re-quantification rule the discovery cohort follows, and it is admitted anyway for one
specific reason: the test applied here is **competitive and within-study**. The RFO set and
its background are drawn from the same matrix, quantified the same way, so a difference in
reference or pipeline shifts target and background together and largely cancels. The same
foreign matrix would **not** be admissible in the pooled effect-size meta-analysis, and it
is not used there.

Only the 21 panel genes need mapping from Solyc to Entrez; the background is every other row
of their matrix, so no genome-wide identifier mapping is involved. The mapping will be
recorded gene by gene, and panel members that cannot be mapped unambiguously are dropped and
reported rather than guessed.

Two further things must be said in advance about this dataset. Its treatment is **deficit
irrigation at 50% of the well-watered volume**, not water withholding — and this project's
own subgroup analysis found the RFO response to be a water-withholding response that the
three deficit-irrigation strata did not show. A null result here is therefore the outcome
our own moderator analysis predicts, and it will be reported as a test of that moderator
claim rather than as a failure of the pathway result. Second, it shares the landrace
'Lucariello' with PRJNA931112; the two differ in year, in the other genotypes and in tissue
coverage, and no accession links them, but neither can be shown to be free of shared plant
material, so they are **not** treated as two independent observations of that landrace.

Leaf and root are analysed separately. Leaf is the pre-registered comparison; root is
exploratory, matching how root is treated in the discovery cohort.

## Hypothesis

The RFO biosynthetic arm is displaced **upward** relative to an expression-matched
background in drought-versus-control leaf comparisons. One-sided; the direction is fixed by
the discovery result and may not be reversed after the fact.

## Gene set, background and test — all inherited unchanged

- **Target set**: `rfo_biosynthetic` exactly as declared in `config/geneset_sets.csv`
  (9 genes). Members not estimable in a study are dropped and reported, never imputed as
  zero, and the set is not re-specified to improve the result.
- **Background**: every gene estimable in that study and absent from the 21-gene panel in
  `config/rfo_panel_full.csv`.
- **Test**: `18_competitive_rotation.R` — CAMERA with the inter-gene correlation estimated
  from the libraries, run per study on voom-transformed counts, with the same
  `filterByExpr` and blocking rules the discovery cohort used. `fry` is reported alongside
  as the self-contained companion but is not the criterion.
- **Negative control**: the `gt8_control` set, under the same test. A positive control
  result accompanied by a positive GT8 result is a failure, not a success.

## Eligibility, decided now

- Leaf tissue, drought against a concurrently sampled control.
- At least 3 libraries per arm after any blocking restriction.
- PRJNA1380757 qualifies as 6 against 6 with genotype (`WT`, `OE15`, `OE18`) blocked. The
  wild-type arm alone is 2 against 2 and is therefore **not** analysed separately.
- PRJNA931112 qualifies as 6 against 6 with cultivar (`Corbarino`, `Lucariello`) blocked,
  **conditional** on the check below.

## Condition on PRJNA931112

Its 12 libraries carry only 2 BioSamples — one per cultivar — so the records do not
distinguish the control from the drought libraries, and the treatment assignment rests
entirely on the experiment alias (`Ct1-3`, `Dr1-3`). Before it is used, the assignment must
be corroborated from the counts by the same procedure applied to PRJNA1164923: the alias
groups must separate on a within-study expression axis. If they do not separate, the study
is reported as unusable and excluded; it is not silently relabelled.

## Success criterion

Declared in advance, and all three parts must hold:

1. The biosynthetic arm's one-sided CAMERA p-value is below 0.05 in the combined
   out-of-sample analysis (Stouffer over the qualifying studies, as in the discovery
   cohort).
2. Its direction is **Up** in a majority of qualifying studies.
3. The GT8 negative control does not reach p < 0.05 in the same analysis.

The criterion is evaluated on the **water-withholding** datasets (PRJNA1380757, and
PRJNA931112 if its condition is met). Zenodo 15040335 is reported alongside but scored
separately, because its deficit-irrigation treatment is the condition our own subgroup
analysis says does not carry the response; it cannot count towards success and a null there
cannot count as failure.

Anything else is a failure to replicate and will be reported as one. With so few studies the
test is underpowered, so a null result bounds the effect rather than refuting it — but a
null result may not be re-described as support.

## What this cannot establish

Two studies, one of them dominated by transgenic overexpression lines, do not make a
validated biomarker. A positive outcome raises confidence that the discovery result is not
an artefact of the discovery corpus; it says nothing about raffinose accumulation,
biosynthetic flux, or drought tolerance, none of which are measured here.
