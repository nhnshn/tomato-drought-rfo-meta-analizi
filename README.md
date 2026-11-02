# Tomato drought RNA-seq meta-analysis — RFO / galactinol pathway

Analysis code for a reproducible meta-analysis of drought response in tomato
(*Solanum lycopersicum*), focused on the raffinose family oligosaccharide (RFO)
biosynthetic arm. 63 candidate BioProjects were screened; 27 were reanalysed from raw
reads (30 strata, 435 libraries), and two were held out as an out-of-sample replication
cohort.

All reported results come from **one sealed end-to-end run** of commit `b6319cb`, tagged
`sealed-run-2026-11-01` (2026-11-01, 20:01:59–21:56:36 UTC; 44 tracked files, 61 inputs,
174 outputs, exit 0).

**This repository holds code only.** The count matrices, the frozen result tables, the
figures and the run logs are deposited separately at Zenodo under ⟨DOI:10.5281/zenodo.23873311⟩. The manuscript
text is not distributed in either place.

## Layout

| path | contents |
|---|---|
| `Makefile` | the whole pipeline; `make complete` runs every stage and seals the run |
| `scripts/` | quantification, matrix building, differential expression, meta-analysis, set tests, sensitivity and calibration, figures, prose checker |
| `config/` | sample sheet, gene-set definitions, genotype arm definitions, prose-claim registry |
| `docs/RUN_MANIFEST.json` | manifest of the sealed run: per-stage commands, start and end times, exit codes, SHA-256 of every input and output |

## Reproducing

The pipeline starts from gene-level count matrices, not FASTQ. Download
`tomato-drought-count-matrices.tar.gz` from the Zenodo deposit, unpack it so the matrices
sit under `matrices/`, then:

```
make complete
```

`make complete` chains `all → sensitivity → robustness → out_of_sample →
background_sensitivity → calibration → seal`. The seal step refuses to finish if any
tracked code or input file changed while the run was in flight, and writes
`logs/RUN_MANIFEST.json`. Comparing that manifest against `docs/RUN_MANIFEST.json` is the
reproduction check.

`make prose` checks every quantitative sentence of the manuscript against the frozen
tables using the claim registry in `config/prose_claims.csv` (67 claims). It needs the
manuscript, which is not distributed here, so it cannot be run from a clone; the registry
is kept because it documents what was checked and against which table and column.

Environment: Python 3.11 with pandas 2.3.3, numpy 2.4.6, scipy 1.17.1, statsmodels 0.15.0,
matplotlib 3.11.2, pymare 0.0.13; R 4.5.3 with limma 3.68.5, edgeR 4.10.7, statmod 1.5.2.
Quantification used salmon 2.8.0 (k = 31) against the ITAG4.1 cDNA reference, decoy-free.

## What the analysis finds

The RFO biosynthetic arm (7 testable genes) has a mean pooled log2 fold change of
**+1.34** across the 13 concurrent-control leaf/shoot strata, against an
expression-matched background null of −0.32 to +0.37; the variance-inflation-corrected
empirical p is **0.0028**, and the correlation-aware rotation test combined across strata
gives **p = 0.00047**. A related glycosyltransferase control set (GT8) sits inside its own
null (−0.05, p = 0.94).

## Limits the tables do not remove

- Of the 19 testable panel genes, only one galactinol synthase (Solyc02g084980) clears
  genome-wide FDR; under the HKSJ variance estimator no gene does.
- Panel, cohort and exclusion decisions were all fixed after the results were seen, so no
  p-value computed inside this corpus is a calibrated error rate.
- Sample labels, control harvest times and biological-replicate independence rest on
  archive metadata and have not been verified against the originating laboratories.
- The calibration simulation tests specific synthetic conditions; it does not calibrate
  very small p-values.

## Post-seal note

`config/prose_claims.csv` is the only tracked file that changed after the seal: the claim
registry can only be written once the prose exists, which is after the run. It enters no
result. The manifest records its hash at run start (`5ecc71e4499e1aec…`); the version here
is `88d2f6bf899c9866…`. The other 43 tracked files are byte-identical to the sealed run.
