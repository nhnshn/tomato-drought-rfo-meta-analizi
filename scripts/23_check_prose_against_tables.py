#!/usr/bin/env python
"""Check every headline number in the prose against the table that produces it.

The package already closes two links in the provenance chain: the SHA-256 inventory
proves an exported file is a faithful copy, and 20_check_export_currency.py proves an
exported table is byte-identical to the pipeline output behind it. Nothing checked the
third link -- that the numbers written in the prose are the numbers in those tables --
and that is exactly where the 1 November audit found 32 stale values, four of which
reversed the direction of a reported result.

A claim is a row of the claims file:

  id, document, locale, pattern, table, query, column, kind, tol

`pattern` is a regular expression over the document with ONE capture group holding the
written value. `query` is a pandas query selecting a single row of `table`; `column`
names the cell the claim refers to. `kind` is `value` (compare numerically within
`tol`), `count` (integer equality) or `rows` (the written value is the number of rows
the query returns). `locale` is `en` or `tr` and fixes how the written number is read. A pattern that matches nothing is reported as MISSING: a claim
whose anchor has been edited away is a stale claim, not a passing one.

Exit status is non-zero if any claim fails OR if a claim's document was not
supplied, so `make` stops both on a stale manuscript and on a partial document set.
"""
import argparse, os, re, sys
import pandas as pd


WORDS = {"one": 1, "two": 2, "three": 3, "four": 4, "five": 5, "six": 6, "seven": 7,
         "eight": 8, "nine": 9, "ten": 10, "eleven": 11, "twelve": 12, "thirteen": 13}


def norm(s, locale):
    """Parse a number written in the document's own numeric convention.

    The convention is declared per claim rather than guessed: "0.446" is four hundred
    and forty-six thousandths in the English draft and four hundred and forty-six in a
    Turkish one, and a checker that guesses wrong reports false failures, which is worse
    than no checker at all.
    """
    t = s.strip().lower()
    if t in WORDS:
        return float(WORDS[t])
    t = t.replace("\u2212", "-").replace("\u2009", "").replace("\u00a0", "").replace(" ", "")
    if locale == "tr":                       # 21.967,5  ->  21967.5
        t = t.replace(".", "").replace(",", ".")
    else:                                    # 5,346.2   ->  5346.2
        t = t.replace(",", "")
    return float(t)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--claims", required=True)
    ap.add_argument("--tables", required=True, help="directory holding the source tables")
    ap.add_argument("--docs", required=True, nargs="+", help="document files to check")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    # Markdown wraps lines at arbitrary columns, so a claim's anchor text is often
    # split across a newline. Collapse runs of whitespace before matching so that a
    # pattern can be written against the sentence as it reads.
    docs = {os.path.basename(p): re.sub(r"\s+", " ", open(p, encoding="utf-8").read())
            for p in a.docs}
    claims = pd.read_csv(a.claims, dtype=str).fillna("")
    cache, results = {}, []

    for _, c in claims.iterrows():
        doc = docs.get(c.document)
        if doc is None:
            # A claim whose document was not supplied is unverified, not satisfied. Treating
            # it as a pass would let `make prose` succeed on a partial document set and so
            # report a clean manuscript it never read.
            results.append(dict(id=c.id, document=c.document, verdict="NO-DOCUMENT",
                                written="", expected="",
                                detail="document not supplied; claim could not be checked"))
            continue
        if c.table not in cache:
            p = os.path.join(a.tables, c.table)
            cache[c.table] = pd.read_csv(p, sep="\t" if p.endswith(".tsv") else ",")
        T = cache[c.table]
        sub = T.query(c.query) if c.query else T

        if c.kind == "rows":
            expected = float(len(sub))
        else:
            if len(sub) != 1:
                results.append(dict(id=c.id, document=c.document, verdict="BAD-QUERY",
                                    written="", expected="",
                                    detail=f"query selected {len(sub)} rows, expected 1"))
                continue
            expected = float(sub.iloc[0][c.column])

        m = list(re.finditer(c.pattern, doc))
        if not m:
            results.append(dict(id=c.id, document=c.document, verdict="MISSING",
                                written="", expected=f"{expected:g}",
                                detail="pattern not found; the anchor text has changed"))
            continue

        tol = float(c.tol) if c.tol else 0.0
        for hit in m:
            try:
                written = norm(hit.group(1), c.locale or "en")
            except ValueError:
                results.append(dict(id=c.id, document=c.document, verdict="UNPARSED",
                                    written=hit.group(1), expected=f"{expected:g}",
                                    detail="captured text is not a number"))
                continue
            ok = (written == expected) if c.kind == "count" else abs(written - expected) <= tol
            results.append(dict(id=c.id, document=c.document,
                                verdict="OK" if ok else "STALE",
                                written=f"{written:g}", expected=f"{expected:g}",
                                detail="" if ok else f"prose says {written:g}, "
                                                     f"{c.table}:{c.column} says {expected:g}"))

    R = pd.DataFrame(results)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    R.to_csv(a.out, index=False)
    bad = R[R.verdict != "OK"]
    for d, g in R.groupby("document"):
        n_bad = (g.verdict != "OK").sum()
        print(f"{d}: {len(g)} claims checked, {n_bad} failing")
    if len(bad):
        print()
        print(bad[["id", "document", "verdict", "written", "expected"]].to_string(index=False))
    return 1 if len(bad) else 0


if __name__ == "__main__":
    sys.exit(main())
