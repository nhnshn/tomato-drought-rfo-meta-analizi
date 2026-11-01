#!/usr/bin/env python3
"""Join the two background definitions of the sample-level test into one table.

camera's background is every other measured gene, so non-target panel members sit
in it; the stratum-level contrast excludes the whole panel. Running the rotation
both ways answers how much that choice matters, but only if the two runs end up
in one file that can be read side by side.
"""
import argparse
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--default", required=True, help="camera's own background")
ap.add_argument("--dropped", required=True, help="non-target panel genes removed")
ap.add_argument("--out", required=True)
a = ap.parse_args()

d = pd.read_csv(a.default);  d["background"] = "all measured genes (camera default)"
x = pd.read_csv(a.dropped);  x["background"] = "non-target panel genes removed"
out = pd.concat([d, x], ignore_index=True)
out.to_csv(a.out, index=False)

key = [c for c in ("gene_set", "cohort") if c in out.columns]
arm = out[out.gene_set == "rfo_biosynthetic"] if "gene_set" in out else out
print(f"{len(out)} rows -> {a.out}")
print(arm[key + [c for c in ("combined_p_up", "bh_within_cohort") if c in out.columns]]
      .to_string(index=False))
