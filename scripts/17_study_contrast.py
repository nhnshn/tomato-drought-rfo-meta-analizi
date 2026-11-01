#!/usr/bin/env python3
"""Study-level contrast of the RFO biosynthetic arm against an expression-matched background.

Why this test exists
--------------------
`11_geneset_test.py` compares the arm's pooled effect with a null built by resampling
individual genes. That null treats genes as exchangeable draws, so co-expression inside the
arm makes the null too narrow; the variance-inflation correction applied there widens it
only approximately, using a single average correlation.

This script avoids the problem by construction rather than by correction. The unit of
analysis is the *stratum*, not the gene. Within each stratum the arm's mean log2FC and a
matched background's weighted mean log2FC are computed from the same study, and their
difference is one observation. Whatever shared response couples the arm's genes inside a
stratum is absorbed into that single number; the spread *across* strata then carries the
uncertainty. No inter-gene correlation term is needed because no inter-gene independence is
assumed.

What it does not give you
-------------------------
The strata are assumed independent and approximately normal, and the t interval ignores the
sampling error of each stratum's own effect estimates. The cohorts and the target-gene set
were fixed after earlier results were seen, so the p-value is a descriptive summary of this
cohort, not a calibrated error rate for a fresh hypothesis. Report it as such.

Independence is also an empirical claim about the cohort, not a formality: two strata that
re-sequence the same BioSamples are not two observations. `--dependent-group` declares such
a set so that a sensitivity cohort can collapse it.
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
from scipy.stats import t as student_t


def control_expression(samplesheet, matrix_dir, studies, tissue, restrict):
    """Median control log2 CPM per gene, averaged over the reference studies.

    Only control libraries enter, so the matching variable cannot be a function of the
    drought response it is used to control for.
    """
    cols, used = [], []
    for study in studies:
        path = os.path.join(matrix_dir, f"{study}_counts.csv.gz")
        if not os.path.exists(path):
            raise SystemExit(f"missing count matrix for reference study {study}: {path}")
        mat = pd.read_csv(path, index_col=0)
        q = samplesheet[(samplesheet.study == study)
                        & (samplesheet.tissue_h == tissue)
                        & (samplesheet.treatment_class == "control")
                        & (samplesheet.sample_id.isin(mat.columns))]
        if study in restrict:
            q = q[q.timepoint.isin(restrict[study])]
        ctl = q.sample_id.tolist()
        if len(ctl) < 3:
            raise SystemExit(f"{study}: only {len(ctl)} control libraries, need >= 3")
        X = mat[ctl].astype(float)
        lib = X.sum(axis=0)
        if not (lib > 0).all():
            raise SystemExit(f"{study}: a control library has zero counts")
        cols.append(np.log2(X.div(lib, axis=1) * 1e6 + 1).median(axis=1).rename(study))
        used.extend(dict(study=study, sample_id=s, library_sum=float(lib[s])) for s in ctl)
    expr = pd.concat(cols, axis=1, join="inner")
    expr["expression"] = expr.mean(axis=1)
    return expr.dropna(), pd.DataFrame(used)


def matched_weights(expr, pool, targets, n_neighbours):
    """Frequency weights over background genes, matched to targets on control expression.

    Sampling is with replacement across targets: a background gene that is the nearest
    neighbour of two targets is counted twice. Ties break on gene id so the weights do not
    depend on row order.
    """
    rows, weights = [], pd.Series(0.0, index=pool)
    ref = expr.loc[pool, "expression"]
    for target in targets:
        d = (ref - expr.loc[target, "expression"]).abs()
        nearest = d.sort_values(kind="stable").iloc[:n_neighbours]
        for gene, dist in nearest.items():
            weights.loc[gene] += 1.0 / (len(targets) * n_neighbours)
            rows.append(dict(n_neighbours=n_neighbours, target=target, gene=gene,
                             match_distance=float(dist)))
    if not np.isclose(weights.sum(), 1.0):
        raise SystemExit(f"matching weights sum to {weights.sum()}, expected 1")
    return weights, pd.DataFrame(rows)


def summarise(diffs, label, n_neighbours):
    """One-sample t summary of the per-stratum differences."""
    k = len(diffs)
    mean = float(diffs.mean())
    se = float(diffs.std(ddof=1) / np.sqrt(k))
    crit = float(student_t.ppf(0.975, k - 1))
    return dict(n_neighbours=n_neighbours, cohort=label, k=k, delta=mean,
                ci_low=mean - crit * se, ci_high=mean + crit * se,
                p_greater=float(student_t.sf(mean / se, k - 1)),
                n_positive=int((diffs > 0).sum()))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--long", required=True,
                    help="per_study_effects_long.csv.gz from the meta-analysis")
    ap.add_argument("--sets", required=True, help="config/geneset_sets.csv")
    ap.add_argument("--panel", required=True, help="config/rfo_panel_full.csv")
    ap.add_argument("--samplesheet", required=True, help="config/samplesheet.csv")
    ap.add_argument("--matrix-dir", required=True, help="directory of per-study count matrices")
    ap.add_argument("--target-set", default="rfo_biosynthetic",
                    help="gene set in --sets whose mean is the target")
    ap.add_argument("--tissue", default="leaf/shoot")
    ap.add_argument("--reference-studies", default="",
                    help="comma-separated strata whose controls define the matching variable; "
                         "default is every stratum in --long that has a count matrix")
    ap.add_argument("--restrict-timepoints", default="",
                    help="STUDY:tp1|tp2,STUDY:tp3 - limit a reference study's controls")
    ap.add_argument("--dependent-group", default="",
                    help="comma-separated strata known to share biological samples; a "
                         "sensitivity cohort keeps only the first of them")
    ap.add_argument("--neighbours", type=int, default=20)
    ap.add_argument("--neighbour-sensitivity", default="10,50,300")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    long = pd.read_csv(args.long)
    if long.duplicated(["study", "gene"]).any():
        raise SystemExit("--long has duplicate study/gene rows")
    strata = sorted(long.study.unique())

    sets = pd.read_csv(args.sets)
    panel = set(pd.read_csv(args.panel).gene)
    samplesheet = pd.read_csv(args.samplesheet)

    declared = sets[sets["set"] == args.target_set].gene.tolist()
    # A gene enters the target mean only if every stratum estimated it: an unbalanced
    # target mean would change composition from stratum to stratum.
    seen = long.groupby("gene").study.nunique()
    complete = set(seen[seen == len(strata)].index)
    targets = [g for g in declared if g in complete]
    dropped = [g for g in declared if g not in complete]
    if len(targets) < 3:
        raise SystemExit(f"only {len(targets)} of {len(declared)} target genes are complete")

    # Which declared member is estimable where. Reported in full because a set mean whose
    # membership changes between strata is not the same quantity from stratum to stratum,
    # and a reader cannot judge that from the count of tested members alone.
    seen_in = (long[long.gene.isin(declared)]
               .assign(measured=True)
               .pivot(index="gene", columns="study", values="measured")
               .reindex(index=declared, columns=strata)
               .notna())
    measurability = seen_in.astype(int).reset_index()
    measurability.insert(1, "n_strata_measured", seen_in.sum(axis=1).values)
    measurability.insert(2, "in_target_set", [g in targets for g in declared])

    refs = ([s.strip() for s in args.reference_studies.split(",") if s.strip()]
            or [s for s in strata
                if os.path.exists(os.path.join(args.matrix_dir, f"{s}_counts.csv.gz"))])
    restrict = {}
    for item in filter(None, (x.strip() for x in args.restrict_timepoints.split(","))):
        study, tps = item.split(":", 1)
        restrict[study] = tps.split("|")

    expr, control_used = control_expression(samplesheet, args.matrix_dir, refs,
                                            args.tissue, restrict)
    # The background is fixed before any contrast is computed, and excludes the whole
    # RFO panel so that a target's neighbour can never be another pathway gene.
    pool = sorted((complete - panel).intersection(expr.index))
    missing = [g for g in targets if g not in expr.index]
    if missing:
        raise SystemExit(f"targets absent from the control expression table: {missing}")
    if len(pool) < 50 * len(targets):
        raise SystemExit(f"background pool of {len(pool)} is too small for matching")

    E = long.pivot(index="study", columns="gene", values="log2FoldChange").loc[strata]

    cohorts = {"full": strata}
    dep = [s.strip() for s in args.dependent_group.split(",") if s.strip()]
    present_dep = [s for s in dep if s in strata]
    if len(present_dep) > 1:
        cohorts["independent_samples"] = [s for s in strata
                                          if s not in present_dep[1:]]

    contrasts, neighbours, summary, loo = [], [], [], []
    settings = [args.neighbours] + [int(x) for x in args.neighbour_sensitivity.split(",")
                                    if x.strip()]
    for nb in settings:
        weights, rows = matched_weights(expr, pool, targets, nb)
        neighbours.append(rows)
        bg = E[pool].dot(weights)
        tg = E[targets].mean(axis=1)
        diff = tg - bg
        for s in strata:
            contrasts.append(dict(n_neighbours=nb, study=s, target_mean=float(tg[s]),
                                  background_mean=float(bg[s]), difference=float(diff[s])))
        for label, group in cohorts.items():
            summary.append(summarise(diff.loc[group], label, nb))
            if nb == args.neighbours:
                for omitted in group:
                    r = summarise(diff.loc[[g for g in group if g != omitted]], label, nb)
                    loo.append(dict(cohort=label, omitted_study=omitted, **{
                        k: v for k, v in r.items() if k not in ("cohort", "n_neighbours")}))

    pd.concat(neighbours).to_csv(os.path.join(args.out, "PUBLIC_matched_neighbours.csv"),
                                 index=False)
    pd.DataFrame(contrasts).to_csv(os.path.join(args.out, "PUBLIC_study_contrasts.csv"),
                                   index=False)
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(os.path.join(args.out, "PUBLIC_study_contrast_summary.csv"), index=False)
    pd.DataFrame(loo).to_csv(os.path.join(args.out, "PUBLIC_study_contrast_loo.csv"),
                             index=False)
    control_used.to_csv(os.path.join(args.out, "control_libraries_used.csv"), index=False)
    measurability.to_csv(os.path.join(args.out, "PUBLIC_target_measurability.csv"), index=False)
    pd.DataFrame({"gene": pool,
                  "control_expression": expr.loc[pool, "expression"].values}
                 ).to_csv(os.path.join(args.out, "eligible_background.csv"), index=False)
    with open(os.path.join(args.out, "study_contrast_spec.json"), "w") as fh:
        json.dump(dict(target_set=args.target_set, targets=targets,
                       targets_dropped_incomplete=dropped, n_strata=len(strata),
                       strata=strata, reference_studies=refs,
                       restrict_timepoints=restrict, tissue=args.tissue,
                       background_size=len(pool), neighbours=args.neighbours,
                       neighbour_sensitivity=settings, cohorts=cohorts,
                       dependent_group=present_dep, alternative="greater"), fh, indent=2)

    print(f"targets {len(targets)}/{len(declared)} complete across {len(strata)} strata"
          + (f"; dropped {dropped}" if dropped else ""))
    print(f"background {len(pool)} genes matched on controls from {len(refs)} strata")
    print(summary_df[summary_df.n_neighbours == args.neighbours].to_string(index=False))


if __name__ == "__main__":
    main()
