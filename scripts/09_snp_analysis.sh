#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# 09_snp_analysis.sh
# Core SNP analysis to validate clonal expansion across
# One Health compartments using Snippy
# ═══════════════════════════════════════════════════════════════

set -euo pipefail

ASSEMBLIES="data/assemblies"
SNP_OUT="results/snp_analysis"
METADATA="data/sample_metadata.csv"

mkdir -p ${SNP_OUT}

# ── Cross-niche STs to analyse ──
# ST147: Clinical (Kumasi) + Environmental (Accra)
# ST101: Clinical (Accra) + Animal (Accra)
# ST15:  Clinical (Kumasi, Accra)
# ST17:  Clinical + Environmental
# ST39:  Clinical (Kumasi) + Environmental (Accra)
# ST1427: Clinical + Environmental

CROSS_NICHE_STS=("ST147" "ST101" "ST15" "ST17" "ST39" "ST1427")

for ST in "${CROSS_NICHE_STS[@]}"; do
    echo "═══════════════════════════════════════"
    echo "Processing ${ST}"
    echo "═══════════════════════════════════════"
    
    ST_DIR="${SNP_OUT}/${ST}"
    mkdir -p ${ST_DIR}/snippy_out
    
    # Identify isolates belonging to this ST from metadata
    # (Assumes metadata CSV has columns: isolate_id, ST, source)
    ISOLATES=$(awk -F',' -v st="${ST}" '$2 == st {print $1}' ${METADATA})
    
    # Select reference: use the first isolate or a public reference
    REF=$(echo ${ISOLATES} | awk '{print $1}')
    REF_FASTA="${ASSEMBLIES}/${REF}.fasta"
    
    echo "Reference: ${REF}"
    echo "Isolates: ${ISOLATES}"
    
    # Run Snippy for each isolate against the reference
    for ISOLATE in ${ISOLATES}; do
        if [ "${ISOLATE}" != "${REF}" ]; then
            echo "  Snippy: ${ISOLATE} vs ${REF}"
            snippy \
                --outdir ${ST_DIR}/snippy_out/${ISOLATE} \
                --ref ${REF_FASTA} \
                --ctgs ${ASSEMBLIES}/${ISOLATE}.fasta \
                --cpus 4 \
                --force
        fi
    done
    
    # Generate core SNP alignment
    snippy-core \
        --ref ${REF_FASTA} \
        --prefix ${ST_DIR}/core \
        ${ST_DIR}/snippy_out/*
    
    # Calculate pairwise SNP distances
    snp-dists ${ST_DIR}/core.full.aln > ${ST_DIR}/${ST}_snp_distances.tsv
    
    echo "  SNP distances saved to ${ST_DIR}/${ST}_snp_distances.tsv"
    echo ""
done

echo "SNP analysis complete for all cross-niche STs."
echo ""
echo "Interpretation guide:"
echo "  <21 SNPs: Same transmission cluster (David et al., 2019)"
echo "  <50 SNPs: Likely recent clonal relationship"
echo "  >100 SNPs: Distinct clonal lineages within the same ST"
echo ""
echo "K. pneumoniae mutation rate: ~4.2e-07 subs/site/year"
echo "  = ~2-3 SNPs/genome/year (Windels et al., 2025)"
echo "  21 SNPs ≈ 7-10 years of divergence"
