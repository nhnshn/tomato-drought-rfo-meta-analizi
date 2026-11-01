#!/usr/bin/env bash
# Build a decoy-aware salmon index for tomato ITAG4.1 (SL4.0 genome, Solyc gene IDs).
# Requires: salmon >= 1.5 (tested on 2.8.0), curl, ~16 GB RAM for the decoy index.
set -euo pipefail
REF_DIR=${1:-ref}
THREADS=${2:-16}
# DECOY=0 builds a cDNA-only index. That is what produced every result in this
# repository: the published index carries num_decoys = 0, and the default here matches
# it so a reproducer gets the same index rather than a differently-built one.
# DECOY=1 adds the SL4.0 chromosomes as decoy sequence (needs ~16 GB RAM).
DECOY=${DECOY:-0}
SOL=https://solgenomics.net/ftp/genomes/Solanum_lycopersicum

mkdir -p "$REF_DIR" && cd "$REF_DIR"
[ -f ITAG4.1_cDNA.fasta ]      || curl -fsSL -O $SOL/annotation/ITAG4.1_release/ITAG4.1_cDNA.fasta
[ -f ITAG4.1_descriptions.txt ]|| curl -fsSL -O $SOL/annotation/ITAG4.1_release/ITAG4.1_descriptions.txt

if [ "$DECOY" = "1" ]; then
  [ -f S_lycopersicum_chromosomes.4.00.fa.gz ] || curl -fsSL -O $SOL/assembly/build_4.00/S_lycopersicum_chromosomes.4.00.fa.gz
  gunzip -kf S_lycopersicum_chromosomes.4.00.fa.gz
  grep '^>' S_lycopersicum_chromosomes.4.00.fa | cut -d' ' -f1 | sed 's/^>//' > decoys.txt
  cat ITAG4.1_cDNA.fasta S_lycopersicum_chromosomes.4.00.fa > gentrome.fa
  salmon index -t gentrome.fa -d decoys.txt -i salmon_itag41_k31 -k 31 -p "$THREADS"
  rm -f gentrome.fa S_lycopersicum_chromosomes.4.00.fa
else
  salmon index -t ITAG4.1_cDNA.fasta -i salmon_itag41_k31 -k 31 -p "$THREADS"
fi

# transcript -> gene map (Solyc01g005000.3.1 -> Solyc01g005000)
grep '^>' ITAG4.1_cDNA.fasta | sed 's/^>//' | cut -d' ' -f1 \
  | awk -F'.' '{print $0"\t"$1}' > tx2gene.tsv
wc -l tx2gene.tsv
echo "index ready: $REF_DIR/salmon_itag41_k31"
