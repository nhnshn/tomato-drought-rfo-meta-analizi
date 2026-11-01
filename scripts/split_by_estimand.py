#!/usr/bin/env python3
"""split_by_estimand.py - separate per-stratum DE tables by what their contrast means.

A stratum whose control and drought arms share no timepoint estimates change from a
baseline sample, not a difference from a concurrent control. Pooling the two kinds
answers neither question cleanly, so they are written to separate directories and
meta-analysed separately. 03_de_design_aware.py records the class in each table.

A third class exists and is kept rather than discarded: strata where no library
carries a timepoint and no curated evidence resolves the control timing. Those
cannot be assigned to either estimand, so they are written to their own directory
and pooled apart. Dropping them silently would hide how much of the corpus rests on
metadata that does not say when the control was taken.
"""
import argparse, glob, os, shutil
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--de", required=True)
    ap.add_argument("--tissue", default=None, help="keep only strata ending in this tissue tag")
    ap.add_argument("--out-primary", required=True)
    ap.add_argument("--out-sensitivity", required=True)
    ap.add_argument("--out-unresolved", default=None,
                    help="strata whose control timing is unresolved; omit to drop them")
    args = ap.parse_args()
    dests = {"concurrent_control": args.out_primary,
             "baseline_referenced": args.out_sensitivity,
             "unresolved_timing": args.out_unresolved}
    for d in dests.values():
        if d:
            os.makedirs(d, exist_ok=True)
    counts = {k: 0 for k in dests}
    counts["other"] = 0
    for f in sorted(glob.glob(os.path.join(args.de, "PUBLIC_*_droughtVScontrol.csv"))):
        stratum = os.path.basename(f)[len("PUBLIC_"):-len("_droughtVScontrol.csv")]
        if args.tissue and not stratum.endswith(args.tissue):
            continue
        est = pd.read_csv(f, usecols=["estimand"], nrows=1).estimand.iloc[0]
        if est in dests and dests[est]:
            shutil.copy(f, dests[est])
            counts[est] += 1
        else:
            counts["other" if est not in counts else est] += 1
            print(f"  {stratum}: estimand '{est}' - no destination, left out")
    print(f"primary (concurrent control): {counts['concurrent_control']} strata; "
          f"sensitivity (baseline referenced): {counts['baseline_referenced']} strata; "
          f"unresolved timing: {counts['unresolved_timing']} strata")


if __name__ == "__main__":
    main()
