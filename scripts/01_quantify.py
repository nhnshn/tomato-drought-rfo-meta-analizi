#!/usr/bin/env python3
"""Stream public tomato drought RNA-seq from ENA and quantify with salmon.

One sample at a time: resolve fastq URLs -> download -> salmon quant -> delete fastq.
Only quant.sf + meta_info.json are kept, so peak disk stays at a few GB regardless
of cohort size. Re-running skips samples whose quant.sf already exists.

Sample sheet (CSV) must have: study, sample_id, runs   (runs = ';'-separated SRR/ERR ids)
Optional columns are carried through untouched.

  python 01_quantify.py --sheet samplesheet.csv --index ref/salmon_itag41_k31 \
      --out quant --threads 16 --workers 3
"""
import argparse, concurrent.futures as cf, csv, hashlib, io, json, os, shutil, subprocess, sys, time
import urllib.request, urllib.error

ENA = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession={acc}"
       "&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,library_layout,read_count&format=tsv")
# NCBI's SRA Open Data mirror on S3. Measured from Modal: ~51 MB/s per stream versus
# ~2 MB/s from the ENA fastq endpoint, i.e. the archive - not the network - was the
# bottleneck. The object is .sra, so fasterq-dump reconstitutes the fastq locally.
AWS_SRA = "https://sra-pub-run-odp.s3.amazonaws.com/sra/{run}/{run}"

def ena_runs(acc, tries=4):
    for t in range(tries):
        try:
            with urllib.request.urlopen(ENA.format(acc=acc), timeout=120) as r:
                rows = list(csv.DictReader(io.TextIOWrapper(r, 'utf-8'), delimiter='\t'))
            if rows:
                return rows
        except Exception as e:
            sys.stderr.write(f"[ena] {acc} attempt {t+1}: {e}\n")
        time.sleep(5 * (t + 1))
    return []

def md5sum(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()

def fetch(url, dest, tries=8, max_bytes=0, md5=None, verify="gzip"):
    """Resumable download, verified against ENA's md5 when one is given.

    A transfer that ENA closes early (curl exit 18) resumes with -C -, but a resumed
    file can still end up corrupt; salmon then fails on a truncated gzip stream with
    no useful message. Verifying the md5 here turns that into a retry instead of a
    silently dropped sample.
    """
    def intact():
        """A file already on disk is good if its md5 matches ENA's, or - when ENA
        publishes none - if it is a complete gzip stream."""
        if not os.path.exists(dest) or os.path.getsize(dest) == 0:
            return False
        if md5:
            return md5sum(dest) == md5
        if verify == "gzip":
            return subprocess.run(["gzip", "-t", dest], capture_output=True).returncode == 0
        # .sra objects are not gzip; fasterq-dump refuses to decode a truncated one, so
        # size>0 here and the conversion itself is the integrity gate.
        return True

    for t in range(tries):
        if not max_bytes and intact():
            return os.path.getsize(dest)
        try:
            cmd = ["curl", "-fsSL", "--retry", "3", "--retry-delay", "5",
                   "--speed-limit", "10240", "--speed-time", "120", "-o", dest, url]
            cmd[1:1] = ["-r", f"0-{max_bytes}"] if max_bytes else ["-C", "-"]
            subprocess.run(cmd, check=True)
            if max_bytes or intact():
                return os.path.getsize(dest)
            sys.stderr.write(f"[dl] {url} integrity check failed after fetch, refetching\n")
            os.remove(dest)
            time.sleep(5 * (t + 1))
        except subprocess.CalledProcessError as e:
            # exit 22 on a resumed transfer is HTTP 416: the byte range starts past EOF,
            # i.e. the file is ALREADY complete. Verify before calling it a failure -
            # treating 416 as an error is what silently dropped ~9% of samples.
            if intact():
                return os.path.getsize(dest)
            sys.stderr.write(f"[dl] {url} attempt {t+1}: {e}\n")
            if os.path.exists(dest) and e.returncode == 22:
                os.remove(dest)          # unusable partial that cannot be resumed
            # exit 56 is ENA resetting the TLS connection under load; exponential backoff
            # (capped at 5 min) gives the archive room instead of hammering it.
            time.sleep(min(300, 15 * 2 ** t))
    raise RuntimeError(f"download failed integrity check after {tries} attempts: {url}")

def head_reads(path, n_reads):
    """Validation helper: truncate a fastq.gz to the first n_reads records."""
    out = path[:-3] + ".head.fq"
    with open(out, "wb") as fh:
        p1 = subprocess.Popen(["gunzip", "-c", path], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["head", "-n", str(4 * n_reads)], stdin=p1.stdout, stdout=fh)
        p1.stdout.close(); p2.communicate()
    os.remove(path)
    return out

def fetch_from_aws(run, work, threads, expect_reads=None):
    """Download the .sra object and convert it to fastq. Returns the fastq paths.

    fasterq-dump refuses to decode a truncated .sra, so a successful conversion is a
    stronger integrity check than a checksum on a compressed fastq; when ENA publishes a
    read count we additionally verify it matches.
    """
    sra = os.path.join(work, f"{run}.sra")
    fetch(AWS_SRA.format(run=run), sra, verify="none")
    cmd = ["fasterq-dump", "--split-files", "--threads", str(threads), "--outdir", work, sra]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    os.remove(sra)
    files = sorted(f for f in os.listdir(work) if f.startswith(run) and f.endswith(".fastq"))
    paths = [os.path.join(work, f) for f in files]
    if not paths:
        raise RuntimeError(f"fasterq-dump produced no fastq for {run}")
    if expect_reads:
        with open(paths[0]) as fh:
            n = sum(1 for _ in fh) // 4
        # ENA's read_count is per-spot for some submissions and per-mate for others.
        # paths[0] holds spots, so accept either convention; anything else is a real truncation.
        expect = int(expect_reads)
        if n != expect and not (len(paths) > 1 and n * len(paths) == expect):
            raise RuntimeError(f"{run}: {n} reads after conversion, ENA reports {expect} "
                               f"({len(paths)} mate file(s)) — not a mate-counting difference")
    return paths

def quant_sample(row, args):
    sample, study = row["sample_id"], row["study"]
    odir = os.path.join(args.out, sample)
    if os.path.exists(os.path.join(odir, "quant.sf")):
        return sample, "cached", None
    work = os.path.join(args.tmp, sample)
    os.makedirs(work, exist_ok=True)
    t0 = time.time()
    try:
        r1, r2, single = [], [], []
        for run in [r for r in row["runs"].split(";") if r]:
            if args.source == "aws":
                recs = ena_runs(run)
                expect = (recs[0].get("read_count") if recs else None)
                local = fetch_from_aws(run, work, args.threads, expect)
                if len(local) >= 2:
                    r1.append(local[0]); r2.append(local[1])
                else:
                    single.append(local[0])
                continue
            meta = [m for m in ena_runs(run) if m.get("fastq_ftp")]
            if not meta:
                raise RuntimeError(f"no ENA fastq record for {run}")
            urls = [u for u in meta[0]["fastq_ftp"].split(";") if u]
            # drop the unpaired "_3"/barcode file some submissions carry alongside _1/_2
            paired = [u for u in urls if u.endswith(("_1.fastq.gz", "_2.fastq.gz"))]
            urls = paired if len(paired) >= 2 else urls
            local = []
            for u in urls:
                dest = os.path.join(work, os.path.basename(u))
                fetch("https://" + u, dest)
                if args.head_reads:
                    dest = head_reads(dest, args.head_reads)
                local.append(dest)
            if len(local) >= 2:
                r1.append(local[0]); r2.append(local[1])
            else:
                single.append(local[0])
        cmd = ["salmon", "quant", "-i", args.index, "-l", "A", "-p", str(args.threads),
               "--seqBias", "--validateMappings", "-o", odir]
        if r1 and r2:
            cmd += ["-1"] + r1 + ["-2"] + r2
            if args.gc_bias:
                cmd += ["--gcBias"]
        elif single:
            cmd += ["-r"] + single
        else:
            raise RuntimeError("no reads resolved")
        os.makedirs(odir, exist_ok=True)
        with open(os.path.join(odir, "salmon.log"), "w") as lg:
            try:
                subprocess.run(cmd, check=True, stdout=lg, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError:
                lg.flush()
                tail = open(os.path.join(odir, "salmon.log")).read()[-400:].replace("\n", " | ")
                raise RuntimeError(f"salmon exit != 0 — log tail: {tail}")
        mi = json.load(open(os.path.join(odir, "aux_info", "meta_info.json")))
        stat = dict(sample_id=sample, study=study, percent_mapped=mi.get("percent_mapped"),
                    num_processed=mi.get("num_processed"),
                    library_types=";".join(mi.get("library_types", [])),
                    minutes=round((time.time() - t0) / 60, 1))
        return sample, "ok", stat
    except Exception as e:
        log = os.path.join(odir, "salmon.log")
        keep = open(log).read()[-2000:] if os.path.exists(log) else ""
        shutil.rmtree(odir, ignore_errors=True)
        os.makedirs(args.out, exist_ok=True)
        with open(os.path.join(args.out, "failures.tsv"), "a") as fh:
            fh.write(f"{study}\t{sample}\t{str(e)[:500]}\t{keep[-500:]}\n".replace("\n\t", "\t"))
        return sample, f"FAILED: {e}", None
    finally:
        shutil.rmtree(work, ignore_errors=True)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sheet", required=True)
    ap.add_argument("--index", required=True)
    ap.add_argument("--out", default="quant")
    ap.add_argument("--tmp", default="fastq_tmp")
    ap.add_argument("--threads", type=int, default=8, help="salmon threads per sample")
    ap.add_argument("--workers", type=int, default=2, help="samples processed concurrently")
    ap.add_argument("--source", choices=["aws", "ena"], default="aws",
                    help="aws = NCBI SRA Open Data mirror (.sra + fasterq-dump, ~25x faster); "
                         "ena = ENA fastq endpoint")
    ap.add_argument("--gc-bias", action="store_true")
    ap.add_argument("--head-reads", type=int, default=0,
                    help="VALIDATION ONLY: keep just the first N reads per file")
    ap.add_argument("--limit", type=int, default=0, help="process only the first N samples")
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.sheet)))
    if args.limit:
        rows = rows[:args.limit]
    os.makedirs(args.out, exist_ok=True); os.makedirs(args.tmp, exist_ok=True)
    log_path = os.path.join(args.out, "quant_log.tsv")
    done, stats = 0, []
    with cf.ThreadPoolExecutor(max_workers=args.workers) as ex:
        for sample, status, stat in ex.map(lambda r: quant_sample(r, args), rows):
            done += 1
            if stat:
                stats.append(stat)
            print(f"[{done}/{len(rows)}] {sample}: {status}" +
                  (f" mapped={stat['percent_mapped']:.1f}% reads={stat['num_processed']:,}"
                   f" {stat['minutes']}min" if stat else ""), flush=True)
    if stats:
        write_header = not os.path.exists(log_path)
        with open(log_path, "a", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(stats[0].keys()), delimiter="\t")
            if write_header:
                w.writeheader()
            w.writerows(stats)
    print(f"quant.sf written for {len(stats)} samples -> {args.out}")

if __name__ == "__main__":
    main()
