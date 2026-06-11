#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# 02_plasmid_reconstruction.sh
# Plasmid reconstruction, typing, and MGE detection using MOB-suite
# ═══════════════════════════════════════════════════════════════

set -euo pipefail

ASSEMBLIES="data/assemblies"
MOB_OUT="results/mob_suite"

mkdir -p ${MOB_OUT}/{mob_recon,mob_typer}

# ── Step 1: MOB-recon — Plasmid reconstruction ──
# Classifies contigs as chromosome or plasmid
# Outputs: chromosome.fasta, plasmid_*.fasta, contig_report.txt, mge.report.txt

for assembly in ${ASSEMBLIES}/*.fasta; do
    sample=$(basename ${assembly} .fasta)
    echo "MOB-recon: ${sample}"
    
    mob_recon \
        --infile ${assembly} \
        --outdir ${MOB_OUT}/mob_recon/${sample} \
        --num_threads 4 \
        --force
done

# ── Step 2: MOB-typer — Plasmid characterisation ──
# Assigns replicon types, relaxase types, MPF types,
# predicted mobility, host range, and MASH nearest neighbour

for sample_dir in ${MOB_OUT}/mob_recon/*; do
    sample=$(basename ${sample_dir})
    
    for plasmid_fasta in ${sample_dir}/plasmid_*.fasta; do
        if [ -f "${plasmid_fasta}" ]; then
            plasmid_name=$(basename ${plasmid_fasta} .fasta)
            
            mob_typer \
                --infile ${plasmid_fasta} \
                --outdir ${MOB_OUT}/mob_typer/${sample} \
                --num_threads 4
        fi
    done
done

# ── Step 3: Compile master reports ──

# Master mobtyper report
echo "Compiling master MOB-typer report..."
head -1 $(find ${MOB_OUT}/mob_typer -name "mobtyper_results.txt" | head -1) \
    > ${MOB_OUT}/MASTER_mobtyper.csv

for f in $(find ${MOB_OUT}/mob_typer -name "mobtyper_results.txt"); do
    sample_dir=$(dirname ${f})
    sample=$(basename ${sample_dir})
    tail -n +2 ${f} | while read line; do
        echo -e "${line}\t${sample}"
    done
done >> ${MOB_OUT}/MASTER_mobtyper.csv

# Master contig report
echo "Compiling master contig report..."
head -1 $(find ${MOB_OUT}/mob_recon -name "contig_report.txt" | head -1) \
    > ${MOB_OUT}/MASTER_contig_report.csv

for f in $(find ${MOB_OUT}/mob_recon -name "contig_report.txt"); do
    tail -n +2 ${f}
done >> ${MOB_OUT}/MASTER_contig_report.csv

# Master MGE report
echo "Compiling master MGE report..."
head -1 $(find ${MOB_OUT}/mob_recon -name "mge.report.txt" | head -1) \
    > ${MOB_OUT}/MASTER_mge_report.csv

for f in $(find ${MOB_OUT}/mob_recon -name "mge.report.txt"); do
    sample_dir=$(dirname ${f})
    sample=$(basename ${sample_dir})
    # Remove suffix for clean isolate ID
    isolate_id=$(echo ${sample} | sed 's/_S[0-9]*_L[0-9]*//')
    tail -n +2 ${f} | while read line; do
        echo -e "${line}\t${isolate_id}"
    done
done >> ${MOB_OUT}/MASTER_mge_report.csv

# ── Step 4: Combine all plasmid FASTAs ──
# Create combined FASTA with headers: >IsolateID|ClusterID|contig_info
echo "Creating combined plasmid FASTA..."
python3 - << 'PYEOF'
import os, re, glob

outfile = open("results/all_plasmids_combined.fasta", "w")
count = 0

for fasta in sorted(glob.glob("results/mob_suite/mob_recon/*/plasmid_*.fasta")):
    parts = fasta.split(os.sep)
    sample_dir = parts[-2]
    isolate_id = re.sub(r'_S\d+_L\d+$', '', sample_dir)
    
    fname = os.path.basename(fasta)
    cluster_id = re.sub(r'^plasmid_', '', fname).replace('.fasta', '')
    
    with open(fasta) as f:
        for line in f:
            if line.startswith('>'):
                contig_info = line[1:].strip()
                outfile.write(f">{isolate_id}|{cluster_id}|{contig_info}\n")
                count += 1
            else:
                outfile.write(line)

outfile.close()
print(f"Combined {count} plasmid sequences into all_plasmids_combined.fasta")
PYEOF

echo "MOB-suite pipeline complete."
echo "Outputs:"
echo "  MASTER_mobtyper.csv — plasmid typing results"
echo "  MASTER_contig_report.csv — contig classification"
echo "  MASTER_mge_report.csv — IS elements and transposons"
echo "  all_plasmids_combined.fasta — combined plasmid sequences"
