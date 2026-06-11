#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# 05_blast_spacer_plasmid.sh
# BLAST CRISPR spacers against reconstructed plasmid sequences
# Run on Galaxy (usegalaxy.org) or locally with BLAST+ installed
# ═══════════════════════════════════════════════════════════════

set -euo pipefail

SPACERS="results/crispr/compiled_results/all_spacers.fasta"
PLASMIDS="results/all_plasmids_combined.fasta"
BLAST_OUT="results/blast"

mkdir -p ${BLAST_OUT}

# ── Step 1: Build BLAST database from plasmid sequences ──
makeblastdb \
    -in ${PLASMIDS} \
    -dbtype nucl \
    -out ${BLAST_OUT}/plasmid_db \
    -title "K. pneumoniae plasmid database"

# ── Step 2: Run BLASTn-short ──
# blastn-short is optimised for short query sequences (spacers ~32bp)
blastn \
    -task blastn-short \
    -query ${SPACERS} \
    -db ${BLAST_OUT}/plasmid_db \
    -out ${BLAST_OUT}/spacer_vs_plasmid.tab \
    -evalue 1 \
    -word_size 7 \
    -max_target_seqs 5 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

# ── Step 3: Summary statistics ──
total=$(wc -l < ${BLAST_OUT}/spacer_vs_plasmid.tab)
significant=$(awk '$3 >= 95 && $4 >= 24' ${BLAST_OUT}/spacer_vs_plasmid.tab | wc -l)
perfect=$(awk '$3 >= 99 && $4 >= 30' ${BLAST_OUT}/spacer_vs_plasmid.tab | wc -l)

echo ""
echo "Spacer-plasmid BLAST complete."
echo "  Total hits: ${total}"
echo "  Significant (>=95% id, >=24bp): ${significant}"
echo "  Near-perfect (>=99% id, >=30bp): ${perfect}"
echo "  Output: ${BLAST_OUT}/spacer_vs_plasmid.tab"
