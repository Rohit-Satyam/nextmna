#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${1:-${ROOT}/reference/plasmoDB52}"
mkdir -p "${OUT}"
wget -O "${OUT}/PlasmoDB-52_Pfalciparum3D7_Genome.fasta" \
  https://plasmodb.org/common/downloads/release-52/Pfalciparum3D7/fasta/data/PlasmoDB-52_Pfalciparum3D7_Genome.fasta
wget -O "${OUT}/PlasmoDB-52_Pfalciparum3D7.gff" \
  https://plasmodb.org/common/downloads/release-52/Pfalciparum3D7/gff/data/PlasmoDB-52_Pfalciparum3D7.gff
samtools faidx "${OUT}/PlasmoDB-52_Pfalciparum3D7_Genome.fasta"
cut -f1,2 "${OUT}/PlasmoDB-52_Pfalciparum3D7_Genome.fasta.fai" > "${OUT}/sizes.genome"
printf 'Reference prepared in %s\n' "${OUT}"
