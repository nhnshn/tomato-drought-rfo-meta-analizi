# Tomato drought analysis: execution and provenance

This README accompanies the correction overlay for the fourth analysis package. The existing FROZEN tables and old manifest remain historical results; installing changed code does **not** regenerate or validate them. The RNA-seq analyses, eligibility and study-level data are unchanged by this overlay.

## Dependencies and inputs

The conda environment files describe the Python and quantification dependencies; the R interpreter must separately provide `limma` and `edgeR`. Record and check actual versions. Restore the count archive into `matrices`, stripping its first component. For external testing, put the two external-study count matrices into `rep/matrices` and the external design file at `rep/stratum_design.csv`. These paths are the inputs expected by the Makefile, not new biological data. Preserve the original archive and work in an empty directory.

## Running

- `make all RSCRIPT=<Rscript>` runs the main workflow and seals its manifest; it does **not** run the calibration, external-cohort, or background-sensitivity targets.
- `make complete RSCRIPT=<Rscript>` runs `all`, then external-cohort testing, background sensitivity, calibration and a final seal, in that order. The calibration can be time-consuming. An absent R package, external input or other command failure aborts the relevant stage. Run it in a clean working directory; `complete` does not remove stale files automatically.
- `make currency EXPORT_DIR=<exported-package>` compares exported FROZEN tables and PNG files to files in the working tree. Ad-hoc outputs require a separate check; absence of flagged mismatches does not prove inferential or scientific validity.

Before execution place all inputs in their final locations. The manifest records hashes of code, primary matrices and any external matrices/design present **at start**. The seal checks that these inputs and tracked code still match and records hashes of result files and logs, together with the finish time. Rule logs include a command, timestamps and exit status. The manifest does not establish that a pre-existing file was generated in this run. Preserve the complete logs, environment and original user ZIP alongside the new run; compare regenerated values and PNGs with the manuscript before publication.

## Scientific qualification

The 13-stratum primary cohort is a working classification, not a fully verified collection of time-matched and independent biological units. PRJNA1164923 has a conflicting sample label, PRJNA1307411 and PRJNA1381306 have unresolved sampling-time evidence, and PRJNA472702 and PRJNA644865 share BioSamples. Treat exclusions and alternative labels as exploratory sensitivity analyses, preserve the eight-stratum gene coverage rule or transparently mark analyses as unestimable, and retain study-level outputs. Do not infer biological replicate independence from PCA separation. Discovery-panel selection and gene-set testing were not prospectively registered against all prior results.

The gene-level DL/normal-z, PM/normal-z and HKSJ summaries report different inferential choices and should remain together. A rescaled gene-resampling null based on across-study correlations is not CAMERA; its p-value is exploratory. The external study GT8 control is not testable in either external study, so it does not constitute clean replication.
