#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# 03_crispr_analysis.sh
# CRISPR-Cas system, anti-CRISPR, and self-targeting analysis
# using CRISPRimmunity
# ═══════════════════════════════════════════════════════════════

set -euo pipefail

ASSEMBLIES="data/assemblies"
CRISPR_OUT="results/crispr"

mkdir -p ${CRISPR_OUT}/{per_isolate,compiled_results}

# ── Step 1: Run CRISPRimmunity per isolate ──
# CRISPRimmunity integrates:
#   - PILER-CR, CRISPRCasFinder, CRT for CRISPR array detection
#   - HmmScan (428 HMM profiles) for Cas protein classification
#   - AcRanker for anti-CRISPR candidate prediction
#   - Self-targeting spacer analysis

for assembly in ${ASSEMBLIES}/*.fasta; do
    sample=$(basename ${assembly} .fasta)
    echo "CRISPRimmunity: ${sample}"
    
    # Run CRISPRimmunity
    # (Assumes CRISPRimmunity is installed and in PATH)
    python3 CRISPRimmunity.py \
        --input ${assembly} \
        --output ${CRISPR_OUT}/per_isolate/${sample} \
        --threads 4
done

# ── Step 2: Compile results across all isolates ──

echo "Compiling CRISPR results..."

# Compile Acr results
echo -e "Sample\tacr_candidate_id\tacr_candidate_location\tacr_prot_direction\tacr_candidate_protein_size\tAcr_candidate_in_prophage\tdistance_with_Aca\tneighbor_aca_protein\thomology_with_known_aca\thomology_with_known_aca_evalue\thomology_with_known_aca_identity\thomology_with_known_aca_query_coverage\thomology_with_known_aca_hit_coverage\thomology_with_known_acr\thomology_with_known_acr_evalue\thomology_with_known_acr_identity\thomology_with_known_acr_query_coverage\thomology_with_known_acr_hit_coverage\tdistance_with_HTH\tHTH_protein\tdomain_evalue\tdomain_identity\tdomain_query_coverage\tdomain_hit_coverage\tgenome_crispr_type\tgenome_self-targeting_number\tanti_contig_protospacer_number\tprotospacer_in_prophage\tself-targeting_spacer\tself-targeting_protospacer\tself-targeting_crispr_type" \
    > ${CRISPR_OUT}/compiled_results/all_acr.tsv

for sample_dir in ${CRISPR_OUT}/per_isolate/*; do
    sample=$(basename ${sample_dir})
    acr_file="${sample_dir}/acr_result.tab"
    if [ -f "${acr_file}" ] && [ -s "${acr_file}" ]; then
        tail -n +2 "${acr_file}" | while read line; do
            echo -e "${line}\t${sample}"
        done >> ${CRISPR_OUT}/compiled_results/all_acr.tsv
    fi
done

# Compile self-targeting results
echo -e "Sample\tself_target_data" \
    > ${CRISPR_OUT}/compiled_results/all_self_targeting.tsv

for sample_dir in ${CRISPR_OUT}/per_isolate/*; do
    sample=$(basename ${sample_dir})
    st_file="${sample_dir}/self_target_result.tab"
    if [ -f "${st_file}" ] && [ -s "${st_file}" ]; then
        tail -n +2 "${st_file}" | while read line; do
            echo -e "${sample}\t${line}"
        done >> ${CRISPR_OUT}/compiled_results/all_self_targeting.tsv
    fi
done

# ── Step 3: Extract all spacer sequences ──

echo "Extracting spacer sequences..."

> ${CRISPR_OUT}/compiled_results/all_spacers.fasta

for sample_dir in ${CRISPR_OUT}/per_isolate/*; do
    sample=$(basename ${sample_dir})
    spc_file=$(find ${sample_dir} -name "*_merge.spc" 2>/dev/null | head -1)
    if [ -f "${spc_file}" ]; then
        while read line; do
            if [[ ${line} == ">"* ]]; then
                header=$(echo ${line} | sed 's/^>//')
                echo ">${sample}|${header}"
            else
                echo "${line}"
            fi
        done < "${spc_file}" >> ${CRISPR_OUT}/compiled_results/all_spacers.fasta
    fi
done

# Count results
n_spacers=$(grep -c "^>" ${CRISPR_OUT}/compiled_results/all_spacers.fasta || echo 0)
n_acr=$(tail -n +2 ${CRISPR_OUT}/compiled_results/all_acr.tsv | wc -l)

echo ""
echo "CRISPRimmunity pipeline complete."
echo "  Total spacers extracted: ${n_spacers}"
echo "  Total Acr candidates: ${n_acr}"
echo "  Output directory: ${CRISPR_OUT}/compiled_results/"
