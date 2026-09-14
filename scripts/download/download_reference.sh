#!/usr/bin/env bash
set -euo pipefail
: "${HF_REFERENCE_ROOT:?source config/server.env}"
: "${HF_PACKAGE_ROOT:?source config/server.env}"
target="${HF_REFERENCE_ROOT}/SM_V10_WBPS19"
mkdir -p "$target"
base="https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/schistosoma_mansoni/PRJEA36577"
files=(
  "schistosoma_mansoni.PRJEA36577.WBPS19.genomic.fa.gz"
  "schistosoma_mansoni.PRJEA36577.WBPS19.mRNA_transcripts.fa.gz"
  "schistosoma_mansoni.PRJEA36577.WBPS19.annotations.gff3.gz"
  "schistosoma_mansoni.PRJEA36577.WBPS19.canonical_geneset.gtf.gz"
)
for name in "${files[@]}"; do
  part="$target/$name.part"
  final="$target/$name"
  if [[ ! -s "$final" ]]; then
    curl --fail --location --retry 8 --retry-all-errors -C - -o "$part" "$base/$name"
    mv "$part" "$final"
  fi
  gzip -t "$final"
  gzip -dc "$final" > "$target/${name%.gz}.part"
  mv "$target/${name%.gz}.part" "$target/${name%.gz}"
done
(cd "$target" && sha256sum *.gz *.fa *.gff3 *.gtf | sort -k2 > reference_checksums.sha256)
cp "$HF_PACKAGE_ROOT/config/reference/reference_manifest.json" "$target/reference_manifest.pinned.json"
grep -q '^>' "$target/schistosoma_mansoni.PRJEA36577.WBPS19.genomic.fa"
grep -q '^>' "$target/schistosoma_mansoni.PRJEA36577.WBPS19.mRNA_transcripts.fa"
grep -q $'\tgene\t' "$target/schistosoma_mansoni.PRJEA36577.WBPS19.annotations.gff3"
grep -q 'gene_id' "$target/schistosoma_mansoni.PRJEA36577.WBPS19.canonical_geneset.gtf"
printf '{"status":"REFERENCE_READY","reference_id":"Schistosoma_mansoni_SM_V10_WBPS19","created_utc":"%s"}\n' "$(date -u +%FT%TZ)" > "$target/REFERENCE_READY.json"
