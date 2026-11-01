#!/usr/bin/env python
"""Check that every flag the Makefile passes is one its script actually declares.

A dry run proves the rules expand; it does not prove the commands are valid. An
unknown flag only surfaces when the rule runs, which on this pipeline can be an
hour in. Run this after editing either the Makefile or a script's CLI:

  make -n all > dry.txt && python scripts/check_make_args.py dry.txt

Two directions are checked, because the `--tx2gene` defect showed one is not enough:
a flag the Makefile passes must be declared by the script, AND every flag the script
marks required must be passed. argparse prints required options unbracketed in its
usage line, which is what `required_flags` reads.

Pass several dry-run files to cover targets outside `all` -- `quant` and `reference`
sit outside it, which is exactly why the missing `--tx2gene` went unnoticed.
"""
import re
import shlex
import subprocess
import sys


def required_flags(usage):
    """Flags argparse shows outside brackets in its usage line, i.e. the required ones."""
    out, depth, i = set(), 0, 0
    while i < len(usage):
        ch = usage[i]
        if ch == "[":
            depth += 1
        elif ch == "]":
            depth = max(0, depth - 1)
        elif depth == 0 and usage.startswith("--", i):
            m = re.match(r"--[a-z0-9-]+", usage[i:])
            if m:
                out.add(m.group(0))
                i += len(m.group(0)) - 1
        i += 1
    return out


def main():
    paths = sys.argv[1:] or ["dry.txt"]
    cmds = [l.strip() for pth in paths for l in open(pth)
            if "scripts/" in l and ("python" in l or "bash" in l)]
    seen, bad, n = set(), 0, 0
    for c in cmds:
        try:
            toks = shlex.split(c)
        except ValueError:
            continue
        script = next((t for t in toks if t.startswith("scripts/")), None)
        if not script or not script.endswith(".py"):
            continue
        i = toks.index(script)
        flags = {t for t in toks[i + 1:] if t.startswith("--")}
        key = (script, tuple(sorted(flags)))
        if key in seen:
            continue
        seen.add(key)
        n += 1
        h = subprocess.run([sys.executable, script, "--help"],
                           capture_output=True, text=True)
        declared = set(re.findall(r"(--[a-z0-9-]+)", h.stdout))
        usage = h.stdout.split("usage:", 1)[-1].split("\n\n", 1)[0] if "usage:" in h.stdout else ""
        unsupplied = required_flags(usage) - flags - {"--help"}
        missing = flags - declared
        if missing:
            print(f"FAIL {script} does not declare {sorted(missing)}")
            bad += 1
        if unsupplied:
            print(f"FAIL {script} requires {sorted(unsupplied)} but the rule does not pass them")
            bad += 1
        if not missing and not unsupplied:
            print(f"ok   {script} ({len(flags)} flags)")
    print(f"\ndistinct invocations checked: {n} | mismatches: {bad}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
