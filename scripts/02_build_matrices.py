#!/usr/bin/env python3
"""Aggregate salmon transcript quantifications to ITAG4.1 gene level (tximport-style).

Produces, per study, a genes x samples count matrix using lengthScaledTPM counts
(the tximport default recommended for DESeq2 when transcript lengths differ between
samples), plus a cohort-wide TPM matrix for QC.

  python 02_build_matrices.py --quant quant --sheet samplesheet.csv \
      --tx2gene ref/tx2gene.tsv --out matrices
"""
import argparse, os
import numpy as np, pandas as pd

def load_sample(path):
    df = pd.read_csv(path, sep="\t")
    return df.set_index("Name")[["EffectiveLength", "TPM", "NumReads"]]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quant", default="quant")
    ap.add_argument("--sheet", required=True)
    ap.add_argument("--tx2gene", required=True)
    ap.add_argument("--out", default="matrices")
    args = ap.parse_args()

    sheet = pd.read_csv(args.sheet)
    t2g = pd.read_csv(args.tx2gene, sep="\t", header=None, names=["tx", "gene"]).set_index("tx")["gene"]
    os.makedirs(args.out, exist_ok=True)

    tpm, counts, efflen, kept = {}, {}, {}, []
    for s in sheet["sample_id"]:
        p = os.path.join(args.quant, s, "quant.sf")
        if not os.path.exists(p):
            print(f"[skip] no quant.sf for {s}"); continue
        d = load_sample(p)
        g = d.join(t2g, how="left").dropna(subset=["gene"])
        tpm[s] = g.groupby("gene")["TPM"].sum()
        counts[s] = g.groupby("gene")["NumReads"].sum()
        w = g["TPM"].clip(lower=1e-9)
        efflen[s] = (g.assign(wl=w * g["EffectiveLength"], w=w)
                      .groupby("gene")[["wl", "w"]].sum().pipe(lambda x: x.wl / x.w))
        kept.append(s)

    TPM = pd.DataFrame(tpm); RAW = pd.DataFrame(counts); LEN = pd.DataFrame(efflen)
    TPM, RAW, LEN = TPM.fillna(0), RAW.fillna(0), LEN.fillna(LEN.mean(axis=1).median())

    # tximport lengthScaledTPM: TPM * mean gene length across samples, rescaled to library size
    mean_len = LEN.mean(axis=1)
    lst = TPM.mul(mean_len, axis=0)
    lst = lst.div(lst.sum(axis=0).replace(0, np.nan), axis=1).mul(RAW.sum(axis=0), axis=1).fillna(0)

    # mtime=0 so that a rerun is byte-comparable; gzip otherwise stamps the clock.
    GZ = {"method": "gzip", "mtime": 0}
    TPM.to_csv(os.path.join(args.out, "cohort_gene_tpm.csv.gz"), compression=GZ)
    lst.round(3).to_csv(os.path.join(args.out, "cohort_gene_counts_lengthScaledTPM.csv.gz"),
                        compression=GZ)
    sub = sheet[sheet.sample_id.isin(kept)]
    for study, grp in sub.groupby("study"):
        cols = [c for c in grp["sample_id"] if c in lst.columns]
        lst[cols].round(3).to_csv(os.path.join(args.out, f"{study}_counts.csv.gz"),
                                  compression=GZ)
        grp.to_csv(os.path.join(args.out, f"{study}_samples.csv"), index=False)
    print(f"{lst.shape[0]} genes x {lst.shape[1]} samples; {sub.study.nunique()} studies -> {args.out}")

if __name__ == "__main__":
    main()
