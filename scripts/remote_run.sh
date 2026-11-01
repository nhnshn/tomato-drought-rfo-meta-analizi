#!/usr/bin/env bash
# Remote shard runner: reference -> ENA streaming salmon quant -> gene matrices.
# Expects in the workdir: shard.csv, 01_quantify.py, 02_build_matrices.py
# Expects a Modal Volume mounted at /vol (persists when the sandbox dies).
#
# Durability: a sandbox killed from outside (spend cap, preemption, idle watchdog) is
# destroyed with its filesystem and ./out/ is NEVER harvested. So each finished sample is
# written to /vol as soon as it exists, and a re-run restores from /vol and skips it.
set -euo pipefail
export PATH="${MAMBA_ROOT_PREFIX:-/opt/conda}/bin:$PATH"
command -v salmon >/dev/null || { echo "salmon not on PATH ($PATH)" >&2; exit 127; }
command -v curl   >/dev/null || { echo "curl not on PATH ($PATH)" >&2; exit 127; }

THREADS=${THREADS:-2}
WORKERS=${WORKERS:-8}
VOL=${VOL:-/vol}
SOL=https://solgenomics.net/ftp/genomes/Solanum_lycopersicum/annotation/ITAG4.1_release

mkdir -p out ref quant "$VOL/quant_sf" "$VOL/logs"

echo "== restore from volume"
restored=0
for f in "$VOL"/quant_sf/*.sf.gz; do
  [ -e "$f" ] || break
  s=$(basename "$f" .sf.gz)
  if grep -q "^$s," shard.csv 2>/dev/null || grep -q ",$s," shard.csv 2>/dev/null; then
    mkdir -p "quant/$s" && gunzip -c "$f" > "quant/$s/quant.sf" && restored=$((restored+1))
  fi
done
echo "restored $restored sample(s) already finished in a previous container"

echo "== reference"
curl -fsSL -O --output-dir ref $SOL/ITAG4.1_cDNA.fasta
grep '^>' ref/ITAG4.1_cDNA.fasta | sed 's/^>//' | cut -d' ' -f1 | awk -F'.' '{print $0"\t"$1}' > ref/tx2gene.tsv
salmon index -t ref/ITAG4.1_cDNA.fasta -i ref/salmon_itag41_k31 -k 31 -p 4 2>&1 | tail -2

echo "== quantify ($(($(wc -l < shard.csv)-1)) rows)"
date -u +"start %H:%M:%S"
# push every finished quant.sf to the volume every 60 s; the volume outlives the sandbox
( while sleep 60; do
    for d in quant/*/; do s=$(basename "$d")
      if [ -f "$d/quant.sf" ] && [ ! -f "$VOL/quant_sf/$s.sf.gz" ]; then
        gzip -c "$d/quant.sf" > "$VOL/quant_sf/$s.sf.gz.part" 2>/dev/null \
          && mv "$VOL/quant_sf/$s.sf.gz.part" "$VOL/quant_sf/$s.sf.gz"
      fi
    done
    cp -f quant/quant_log.tsv "$VOL/logs/$(hostname)_quant_log.tsv" 2>/dev/null || true
    cp -f quant/failures.tsv  "$VOL/logs/$(hostname)_failures.tsv"  2>/dev/null || true
  done ) & WATCHER=$!
trap 'kill $WATCHER 2>/dev/null || true' EXIT

python 01_quantify.py --sheet shard.csv --index ref/salmon_itag41_k31 \
    --out quant --tmp fastq_tmp --threads "$THREADS" --workers "$WORKERS" --source "${SOURCE:-aws}"
date -u +"quant done %H:%M:%S"

# final flush before the matrices step
for d in quant/*/; do s=$(basename "$d")
  [ -f "$d/quant.sf" ] && [ ! -f "$VOL/quant_sf/$s.sf.gz" ] && gzip -c "$d/quant.sf" > "$VOL/quant_sf/$s.sf.gz"
done

echo "== matrices"
python 02_build_matrices.py --quant quant --sheet shard.csv --tx2gene ref/tx2gene.tsv --out matrices
cp matrices/*_counts.csv.gz matrices/*_samples.csv out/ 2>/dev/null || true
cp quant/quant_log.tsv quant/failures.tsv out/ 2>/dev/null || true
{ echo "restored_from_volume: $restored";
  echo "quant_sf_in_volume: $(ls "$VOL"/quant_sf/*.sf.gz 2>/dev/null | wc -l)";
  echo "matrix_files: $(ls out/*_counts.csv.gz 2>/dev/null | wc -l)"; } > out/summary.txt
cat out/summary.txt
