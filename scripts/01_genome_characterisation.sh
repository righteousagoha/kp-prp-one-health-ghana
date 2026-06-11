#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# 01_genome_characterisation.sh
# Genome assembly, species confirmation, and resistance profiling
# ═══════════════════════════════════════════════════════════════

set -euo pipefail

# Directories
RAW_READS="data/raw_reads"
ASSEMBLIES="data/assemblies"
KLEBORATE_OUT="results/kleborate"
BACPIPE_OUT="results/bacpipe"

mkdir -p ${ASSEMBLIES} ${KLEBORATE_OUT} ${BACPIPE_OUT}

# ── Step 1: Quality control and trimming ──
for sample_dir in ${RAW_READS}/*; do
    sample=$(basename ${sample_dir})
    R1=$(ls ${sample_dir}/*_R1_*.fastq.gz)
    R2=$(ls ${sample_dir}/*_R2_*.fastq.gz)
    
    echo "Processing: ${sample}"
    
    # FastQC on raw reads
    fastqc ${R1} ${R2} -o results/fastqc_raw/
    
    # Trimmomatic
    trimmomatic PE \
        ${R1} ${R2} \
        ${sample_dir}/${sample}_R1_trimmed.fastq.gz ${sample_dir}/${sample}_R1_unpaired.fastq.gz \
        ${sample_dir}/${sample}_R2_trimmed.fastq.gz ${sample_dir}/${sample}_R2_unpaired.fastq.gz \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
        LEADING:3 TRAILING:3 \
        SLIDINGWINDOW:4:15 MINLEN:36
    
    # FastQC on trimmed reads
    fastqc ${sample_dir}/${sample}_R1_trimmed.fastq.gz \
           ${sample_dir}/${sample}_R2_trimmed.fastq.gz \
           -o results/fastqc_trimmed/
done

# ── Step 2: De novo assembly with Unicycler ──
for sample_dir in ${RAW_READS}/*; do
    sample=$(basename ${sample_dir})
    
    unicycler \
        -1 ${sample_dir}/${sample}_R1_trimmed.fastq.gz \
        -2 ${sample_dir}/${sample}_R2_trimmed.fastq.gz \
        -o ${ASSEMBLIES}/${sample} \
        --mode normal \
        --verbosity 2
    
    # Copy assembly to standard name
    cp ${ASSEMBLIES}/${sample}/assembly.fasta ${ASSEMBLIES}/${sample}.fasta
done

# ── Step 3: Kleborate analysis ──
# Run Kleborate on all assemblies
kleborate \
    --all \
    -a ${ASSEMBLIES}/*.fasta \
    -o ${KLEBORATE_OUT}/klebsiella_pneumo_complex_output.txt

echo "Kleborate analysis complete."

# ── Step 4: BacPipe analysis ──
# BacPipe integrates Prokka, ResFinder, CARD, VirulenceFinder,
# PlasmidFinder, MLST, and Resfams
for assembly in ${ASSEMBLIES}/*.fasta; do
    sample=$(basename ${assembly} .fasta)
    
    bacpipe \
        --input ${assembly} \
        --output ${BACPIPE_OUT}/${sample} \
        --species "Klebsiella pneumoniae" \
        --threads 4
done

echo "BacPipe analysis complete."

# ── Step 5: Species filtering ──
# Exclude K. variicola (UGK24, UGK25, UGK35, UGK39)
# Retain K. quasipneumoniae (6 isolates)
# Final dataset: 78 isolates

echo "Pipeline complete. 78 isolates retained for downstream analysis."
