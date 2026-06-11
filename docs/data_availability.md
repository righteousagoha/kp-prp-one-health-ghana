# Data Availability Statement

## Template for manuscript

Raw sequencing reads have been deposited in the NCBI Sequence Read Archive under
BioProject PRJNA_XXXXXX. Assembled genome sequences are available from the same
BioProject. All analysis code, including bioinformatics pipelines (genome assembly,
plasmid reconstruction, CRISPR-Cas profiling, anti-CRISPR detection, spacer-plasmid
BLAST analysis, and SNP-based clonal validation) and R scripts for statistical
analysis and figure generation, is available at
https://github.com/righteousagoha/kp-prp-one-health-ghana.
Compiled datasets including MOB-suite results, CRISPRimmunity outputs, ABRicate
reports, BLAST results, and anti-CRISPR genomic location data are deposited at
Zenodo (DOI: 10.5281/zenodo.XXXXXXX).

## Steps to complete

1. NCBI BioProject:
   - Go to https://submit.ncbi.nlm.nih.gov/
   - Create a new BioProject (type: Raw sequence reads)
   - Title: "Klebsiella pneumoniae One Health Ghana WGS"
   - Organism: Klebsiella pneumoniae
   - Register BioSamples for all 78 isolates
   - Upload raw FASTQ files to SRA
   - Record the BioProject accession (PRJNA_XXXXXX)

2. GitHub repository:
   - Create at https://github.com/new
   - Name: kp-prp-one-health-ghana
   - Push code (see instructions below)

3. Zenodo archive:
   - Go to https://zenodo.org/deposit/new
   - Upload compiled data files:
     - MASTER_mobtyper.csv
     - MASTER_contig_report.csv
     - MASTER_mge_report.csv
     - all_spacers.fasta
     - all_acr.tsv
     - abricate_ncbi_amr.tab
     - abricate_vfdb.tab
     - abricate_plasmidfinder.tab
     - spacer_vs_plasmid.tab
     - acr_genomic_location.csv
   - Link to GitHub repository
   - Record the DOI
