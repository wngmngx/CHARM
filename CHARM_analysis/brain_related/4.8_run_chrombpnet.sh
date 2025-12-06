#!/usr/bin/env bash
set -euo pipefail

#selected_subclasses = ['ITL23GL', 'ASC', 'CTGL', 'ITL6GL', 'PER', 'ITL5GL', 'PVGA',
#  'NPGL', 'MGL', 'ITL4GL', 'OGC', 'PTGL', 'SSTGA', 'LAMGA', 'OPC',
#  'VLMC', 'VIPGA', 'VEC', 'CLAGL', 'IOL', 'L6bGL', 'VPIA', 'STRGA',
#  'MXD', 'D2MSN', 'OBGA1', 'MSGA', 'D1MSN', 'OBNBL', 'OLFGL', 'RGL',
#  'LSXGA', 'CNUGA', 'IRGL', 'OBDOP', 'CA3GL', 'OBGL', 'DGNBL',
#  'CRC', 'OBGA2', 'ITHGL'pe f -name "$name" -print -quit)", 'CA1GL', 'GRC']

# this is just an example for one cell type
targets=(
  "combine.mm10.frag.LAMGA.tsv"
)

for name in "${targets[@]}"; do
  frag="$(find split.Li2021 -type f -name "$name" -print -quit)"
  if [ -z "$frag" ]; then
    echo "warning: $name not found under split.Li2021" >&2
    continue
  fi

  ct="$(basename "$frag" .tsv)"; ct="${ct##*.}"

  outdir="Li2021output/$ct"
  mkdir -p "$outdir"

  filt="$outdir/$(basename "$frag" .tsv).len_gt10.tsv"
  if [ ! -s "$filt" ] || [ "$frag" -nt "$filt" ]; then
    awk 'BEGIN{OFS="\t"} ($3-$2)>10' "$frag" > "$filt"
  fi

  chrombpnet pipeline \
    -ifrag "$filt" \
    -d ATAC \
    -g data/genome.fa \
    -c data/chr.len \
    -p data/top200k.narrowPeak \
    -n data/output_negatives.bed \
    -fl data/splits/fold_0.json \
    -b data/bias_models/ATAC/scATAC_dermal_fibroblast.h5 \
    -o "$outdir"
done
