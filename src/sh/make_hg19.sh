#!/usr/bin/env bash
# Fetch UCSC primary hg19 chromosomes with upstream integrity checks.
set -euo pipefail
PREFIX=https://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes
test ! -e hg19.fa
mkdir -p .work
WORK=$(mktemp -d "$PWD/.work/hg19.XXXXXXXX")
trap 'rm -rf "$WORK"' EXIT
curl -fL "$PREFIX/md5sum.txt" -o "$WORK/md5sum.txt"
for chromosome in {1..22} X Y M; do
    name="chr${chromosome}.fa.gz"
    curl -fL "$PREFIX/$name" -o "$WORK/$name"
    expected=$(awk -v name="$name" '$2==name || $2=="*"name {print $1}' "$WORK/md5sum.txt")
    test ${#expected} -eq 32
    printf '%s  %s\n' "$expected" "$WORK/$name" | md5sum -c -
    gzip -dc "$WORK/$name" >> "$WORK/hg19.fa"
    rm "$WORK/$name"
done
samtools faidx "$WORK/hg19.fa"
mv "$WORK/hg19.fa.fai" hg19.fa.fai
mv "$WORK/hg19.fa" hg19.fa
