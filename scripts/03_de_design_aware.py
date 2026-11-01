#!/usr/bin/env python3
"""03b_de_design_aware.py - per-stratum differential expression with a design chosen
from the curated metadata instead of a single two-group model everywhere.

Rules, applied in this order and recorded per stratum in the summary:

1. Shared-level restriction. If some control samples sit at timepoints where no
   drought sample exists (a partially baseline-referenced design), the stratum is
   restricted to the timepoints both arms share, provided >=3v3 survive.
2. Baseline-referenced designs. If control and drought share NO timepoint, the
   contrast is not a drought-vs-concurrent-control effect but a change from a
   baseline sample. The stratum is kept, tagged estimand="baseline_referenced",
   and excluded from the primary meta-analysis.
3. Blocking. Any of timepoint / genotype / subtissue that still varies after the
   restriction enters the design as a blocking term, so the condition effect is
   estimated within blocks. Terms that make the design rank-deficient are dropped
   one at a time, least informative first, and the drop is recorded.
4. The contrast is always ["treatment_class", "drought", "control"].

Usage:
  python 03b_de_design_aware.py --matrices matrices_full --samples samplesheet.csv \
      --out de_design --min-reps 3
"""
import argparse, glob, os, re
import numpy as np, pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference

PUBLIC_PREFIX = "PUBLIC_"
BLOCK_CANDIDATES = ["subtissue", "timepoint", "genotype"]


def full_rank(meta, terms):
    """True if an intercept + dummy coding of terms + treatment_class is full rank."""
    X = pd.get_dummies(meta[terms + ["treatment_class"]].astype(str), drop_first=True).astype(float)
    X.insert(0, "intercept", 1.0)
    return np.linalg.matrix_rank(X.values) == X.shape[1] and X.shape[0] > X.shape[1]


def curated_estimand(meta, notes, why):
    """Fall back to the curated control_timing field when the timepoints cannot decide.

    Used both when no sample carries a timepoint and when the timed samples alone do
    not establish a shared level. A blank timepoint is an unknown time, never evidence
    that two arms were sampled together.
    """
    curated = (meta["control_timing"].dropna().unique().tolist()
               if "control_timing" in meta.columns else [])
    ev = (meta["control_timing_evidence"].dropna().unique().tolist()
          if "control_timing_evidence" in meta.columns else [])
    src = ev[0] if ev else "see samplesheet"
    if curated == ["concurrent"]:
        notes.append(f"{why}; concurrency taken from curation: {src}")
        return "concurrent_control"
    if curated == ["baseline_referenced"]:
        notes.append(f"{why}; baseline design taken from curation: {src}")
        return "baseline_referenced"
    if curated == ["unresolved"]:
        # A deliberate curation decision, not an absence of one. Recording it separately
        # keeps "we read the paper and it does not settle the timing" distinguishable
        # from "nobody looked", which matters when the stratum is later re-examined.
        notes.append(f"{why}; left unresolved by curation: {src}")
        return "unresolved_timing"
    notes.append(f"{why} and no curated control timing: cannot establish whether the "
                 "control is concurrent or a baseline")
    return "unresolved_timing"


def dropped_composition(meta, kept_mask):
    """Describe the libraries a shared-level restriction removes, by arm."""
    gone = meta[~kept_mask]
    n_c = int((gone.treatment_class == "control").sum())
    n_d = int((gone.treatment_class == "drought").sum())
    if n_d == 0:
        return f"{len(gone)} librar{'y' if len(gone)==1 else 'ies'}, all control-only"
    if n_c == 0:
        return f"{len(gone)} librar{'y' if len(gone)==1 else 'ies'}, all drought-only"
    return f"{len(gone)} libraries, {n_c} control and {n_d} drought"


def choose_design(meta, min_reps):
    """Return (meta_subset, blocking_terms, estimand, notes)."""
    notes = []
    tp = meta["timepoint"]
    # Absence of a timepoint in both arms is absence of evidence, not evidence that the
    # arms were sampled together: filling the blanks with one label would make every
    # untimed study look concurrently controlled. Such strata fall back to the curated
    # control_timing field, and are left unresolved when that is empty too.
    if tp.isna().all():
        est = curated_estimand(meta, notes, "no timepoint recorded in either arm")
        terms = [c for c in BLOCK_CANDIDATES
                 if c in meta.columns and meta[c].astype(str).nunique() > 1]
        while terms and not full_rank(meta, terms):
            drop = min(terms, key=lambda c: (meta[c].astype(str).nunique(), BLOCK_CANDIDATES.index(c)))
            terms.remove(drop)
            notes.append(f"blocking term '{drop}' dropped: design not full rank with it")
        if terms:
            notes.append("blocked on " + " + ".join(terms))
        return meta, terms, est, notes
    # Some samples carry a timepoint and some do not. A blank is an unknown time, so it
    # is excluded from the shared-level computation rather than collapsed to one "NA"
    # label -- collapsing would let two arms "share" the blank and make an unresolved
    # design look concurrently controlled.
    timed = tp.notna()
    n_untimed = int((~timed).sum())
    if n_untimed:
        notes.append(f"{n_untimed} of {len(meta)} libraries carry no timepoint; they cannot "
                     "establish a shared level and are excluded from that decision")
    ctl = set(tp[timed & (meta.treatment_class == "control")])
    dro = set(tp[timed & (meta.treatment_class == "drought")])
    shared = ctl & dro
    if not shared:
        if n_untimed:
            # the timed samples alone decide nothing; defer to curation, do not assume
            est = curated_estimand(
                meta, notes, "the timed libraries share no timepoint and the rest are untimed")
            terms = [c for c in BLOCK_CANDIDATES
                     if c in meta.columns and meta[c].astype(str).nunique() > 1]
            while terms and not full_rank(meta, terms):
                drop = min(terms, key=lambda c: (meta[c].astype(str).nunique(),
                                                 BLOCK_CANDIDATES.index(c)))
                terms.remove(drop)
                notes.append(f"blocking term '{drop}' dropped: design not full rank with it")
            if terms:
                notes.append("blocked on " + " + ".join(terms))
            return meta, terms, est, notes
        return meta, [], "baseline_referenced", [
            "control and drought share no timepoint: the contrast measures change from a "
            "baseline sample, not a difference from a concurrent control"]
    keep = timed & tp.isin(shared)
    if not keep.all():
        sub = meta[keep]
        n_c = (sub.treatment_class == "control").sum()
        n_d = (sub.treatment_class == "drought").sum()
        if n_c >= min_reps and n_d >= min_reps:
            dropped = sorted(set(tp[timed]) - shared)
            notes.append(f"restricted to timepoints present in both arms; dropped {dropped} "
                         f"({dropped_composition(meta, keep)})")
            meta = sub
        else:
            notes.append(f"timepoints {sorted(set(tp[timed]) - shared)} are not shared by both "
                         f"arms but restricting would leave {n_c}v{n_d}; kept unrestricted "
                         "and flagged")
            return meta, [], "time_confounded", notes
    terms = [c for c in BLOCK_CANDIDATES if c in meta.columns and meta[c].astype(str).nunique() > 1]
    while terms and not full_rank(meta, terms):
        # Drop the least informative term first: the one with the fewest levels. In a
        # genotype panel the sampling stage is nested inside genotype, so dropping the
        # coarser term keeps the pairing that actually blocks the comparison.
        dropped = min(terms, key=lambda c: (meta[c].astype(str).nunique(), BLOCK_CANDIDATES.index(c)))
        terms.remove(dropped)
        notes.append(f"blocking term '{dropped}' ({meta[dropped].astype(str).nunique()} levels) dropped: "
                     f"design not full rank with it")
    if terms:
        notes.append("blocked on " + " + ".join(terms))
    return meta, terms, "concurrent_control", notes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrices", default="matrices_full")
    ap.add_argument("--samples", required=True)
    ap.add_argument("--out", default="de_design")
    ap.add_argument("--min-total-count", type=int, default=10)
    ap.add_argument("--min-reps", type=int, default=3)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    inf = DefaultInference(n_cpus=1)
    sheet = pd.read_csv(args.samples)
    summary, used = [], []

    for cfile in sorted(glob.glob(os.path.join(args.matrices, "*_counts.csv.gz"))):
        study = os.path.basename(cfile).replace("_counts.csv.gz", "")
        counts = pd.read_csv(cfile, index_col=0)
        meta_all = sheet[(sheet.study == study) & sheet.sample_id.isin(counts.columns)
                         & sheet.treatment_class.isin(["control", "drought"])]
        for tissue, m0 in meta_all.groupby("tissue_h"):
            tag = f"{study}__{re.sub(r'[^0-9A-Za-z]+', '', str(tissue))}"
            n_c0 = (m0.treatment_class == "control").sum()
            n_d0 = (m0.treatment_class == "drought").sum()
            if n_c0 < args.min_reps or n_d0 < args.min_reps:
                summary.append(dict(stratum=tag, n_control=n_c0, n_drought=n_d0,
                                    status="skipped: <min reps")); continue
            m, terms, estimand, notes = choose_design(m0, args.min_reps)
            n_c = (m.treatment_class == "control").sum()
            n_d = (m.treatment_class == "drought").sum()
            X = counts[m.sample_id].T.round().astype(int)
            keep = X.sum(axis=0) >= args.min_total_count
            if int(keep.sum()) == 0:
                summary.append(dict(stratum=tag, n_control=n_c, n_drought=n_d, genes_tested=0,
                                    estimand=estimand, design="-",
                                    status="excluded: no gene reaches the count floor")); continue
            md = m.set_index("sample_id").copy()
            for t in terms:
                md[t] = md[t].astype(str)
            design = "~" + " + ".join(terms + ["treatment_class"])
            dds = DeseqDataSet(counts=X.loc[:, keep], metadata=md, design=design,
                               refit_cooks=True, inference=inf, quiet=True)
            dds.deseq2()
            st = DeseqStats(dds, contrast=["treatment_class", "drought", "control"],
                            inference=inf, quiet=True)
            st.summary()
            res = st.results_df.copy()
            res["max_count"] = X.loc[:, keep].max(axis=0)
            res["n_control"], res["n_drought"] = n_c, n_d
            res["stratum"], res["study"], res["estimand"] = tag, study, estimand
            res.index.name = "gene"
            res.to_csv(os.path.join(args.out, f"{PUBLIC_PREFIX}{tag}_droughtVScontrol.csv"))
            # Record exactly which libraries entered this fit and under which terms. Any
            # later analysis that works from the counts rather than from the effect table
            # has to reproduce this restriction, and reconstructing it from the samplesheet
            # would silently drift from the design actually fitted here.
            used.extend(dict(stratum=tag, study=study, sample_id=sid, estimand=estimand,
                             design=design, treatment_class=row.treatment_class,
                             **{t: row[t] for t in terms})
                        for sid, row in m.set_index("sample_id").iterrows())
            sig = int((res.padj < 0.05).sum())
            summary.append(dict(stratum=tag, n_control=n_c, n_drought=n_d,
                                n_control_before=n_c0, n_drought_before=n_d0,
                                genes_tested=int(keep.sum()), n_padj_lt_05=sig,
                                estimand=estimand, design=design,
                                design_notes="; ".join(notes), status="ok"))
            print(f"{tag}: {n_c}v{n_d} {design} [{estimand}] -> {sig} at FDR<0.05", flush=True)
    pd.DataFrame(summary).to_csv(os.path.join(args.out, f"{PUBLIC_PREFIX}de_summary.csv"), index=False)
    pd.DataFrame(used).to_csv(
        os.path.join(args.out, f"{PUBLIC_PREFIX}stratum_design.csv"), index=False)


if __name__ == "__main__":
    main()
