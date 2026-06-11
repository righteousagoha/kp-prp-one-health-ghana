#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# 04_abricate_annotation.sh
# ABRicate annotation of plasmid sequences
# Run on Galaxy (usegalaxy.org) or locally if ABRicate is installed
# ═══════════════════════════════════════════════════════════════

set -euo pipefail

PLASMID_FASTA="results/all_plasmids_combined.fasta"
ABRICATE_OUT="results/abricate"

mkdir -p ${ABRICATE_OUT}

# ── ABRicate with NCBI AMRFinderPlus database ──
abricate \
    --db ncbi \
    --mincov 80 \
    --minid 90 \
    ${PLASMID_FASTA} \
    > ${ABRICATE_OUT}/abricate_ncbi_amr.tab

# ── ABRicate with VFDB database ──
abricate \
    --db vfdb \
    --mincov 80 \
    --minid 90 \
    ${PLASMID_FASTA} \
    > ${ABRICATE_OUT}/abricate_vfdb.tab

# ── ABRicate with PlasmidFinder database ──
abricate \
    --db plasmidfinder \
    --mincov 80 \
    --minid 90 \
    ${PLASMID_FASTA} \
    > ${ABRICATE_OUT}/abricate_plasmidfinder.tab

echo "ABRicate annotation complete."
echo "Results:"
echo "  AMR genes: $(tail -n +2 ${ABRICATE_OUT}/abricate_ncbi_amr.tab | wc -l) hits"
echo "  Virulence genes: $(tail -n +2 ${ABRICATE_OUT}/abricate_vfdb.tab | wc -l) hits"
echo "  Replicon types: $(tail -n +2 ${ABRICATE_OUT}/abricate_plasmidfinder.tab | wc -l) hits"
