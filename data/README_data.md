# Data Files

## Input data (not included in repository; available from NCBI/Zenodo)

- `raw_reads/` — Raw Illumina MiSeq FASTQ files (BioProject: TBD)
- `assemblies/` — Unicycler de novo assemblies (78 isolates)

## Compiled results (included or available on request)

### MOB-suite outputs
- `MASTER_mobtyper.csv` — Plasmid typing for 370 plasmids (replicon types, mobility, host range, MASH neighbours)
- `MASTER_contig_report.csv` — Contig classification (chromosome vs plasmid) for 29,061 contigs
- `MASTER_mge_report.csv` — Mobile genetic element annotations (IS elements, transposons)

### Kleborate outputs
- `Kleborate_Results.xlsx` — Species, MLST, AMR, virulence, K/O typing
- `klebsiella_pneumo_complex_output.txt` — Full Kleborate text output

### CRISPRimmunity outputs
- `all_spacers.fasta` — 1,207 CRISPR spacer sequences
- `all_acr.tsv` — Anti-CRISPR candidate results
- `all_self_targeting.tsv` — Self-targeting spacer results

### ABRicate outputs (from Galaxy)
- `abricate_ncbi_amr.tab` — 347 AMR gene hits on plasmids
- `abricate_vfdb.tab` — 21 virulence gene hits on plasmids
- `abricate_plasmidfinder.tab` — 274 replicon type hits

### BLAST outputs
- `spacer_vs_plasmid.tab` — 7,212 spacer-plasmid BLAST hits

### Derived data
- `plasmid_AMR_cargo.csv` — AMR genes per plasmid cluster
- `acr_genomic_location.csv` — Chromosome vs plasmid location of Acr genes
- `all_plasmids_combined.fasta` — 2,600 concatenated plasmid sequences

## Sample metadata

| Source | Location | Code | n |
|--------|----------|------|---|
| Human (Clinical) | KATH, Kumasi | KATH-* | 11 |
| Human (Clinical) | KBTH, Accra | KBTH-* | 7 |
| Human (Clinical) | HTH, Ho | HTH-* | 2 |
| Animal (Poultry) | Accra | PF-* | 17 |
| Animal (Cattle) | Accra | CTL-*, GT-* | 4 |
| Environment (Wastewater) | Korle Lagoon, Accra | UGK* | 35 |
| Environment (Wastewater) | Accra | WT-* | 2 |
| **Total** | | | **78** |

## Excluded isolates
- UGK24, UGK25, UGK35, UGK39 (*K. variicola*)
