# Methods

## Study design, setting, and sample collection

This cross-sectional study analysed *Klebsiella pneumoniae* isolates collected between September 2021 and May 2024 from three One Health compartments across Ghana. Clinical isolates (n = 20) were obtained from archived collections at three tertiary referral hospitals: Komfo Anokye Teaching Hospital (KATH, Kumasi, Ashanti Region; n = 11), Korle-Bu Teaching Hospital (KBTH, Accra, Greater Accra Region; n = 7), and Ho Teaching Hospital (HTH, Ho, Volta Region; n = 2). These facilities serve as regional referral centres with patient catchment areas extending to neighbouring countries. Animal isolates (n = 21) were collected from poultry farms (n = 17), cattle farms (n = 4), and an abattoir within the Greater Accra Region. Environmental isolates (n = 37) were recovered from the Korle Lagoon watershed (n = 35) and wastewater ponds (n = 2) in Accra. The Korle Lagoon receives untreated and partially treated effluent from Korle-Bu Teaching Hospital and surrounding communities, providing a direct hydrological link between clinical and environmental compartments.

## Bacterial isolation and identification

A total of 82 bacterial isolates were initially included. Isolates were cultured on Simmons citrate inositol agar (SCAI) following pre-enrichment in buffered peptone water. Preliminary identification was performed using biochemical tests including triple sugar iron (TSI), urease, citrate utilisation, and sulphur-indole-motility (SIM) on yellowish colonies obtained from SCAI. Molecular confirmation of *K. pneumoniae* was performed by conventional PCR targeting the capsular polysaccharide rscA gene using the KP-27F3 (5'-GGATATCTGACCAGTCGG-3') and KP-27B3 (5'-GGGTTTTGCGTAATGATCTG-3') primer pair, which amplifies a 176 bp fragment [1]. Prior to whole-genome sequencing, isolate identity was further confirmed by matrix-assisted laser desorption/ionisation time-of-flight mass spectrometry (MALDI-TOF MS; Bruker, Billerica, MA, USA) on colonies grown on chromogenic agar at 37°C for 24 hours.

## DNA extraction and whole-genome sequencing

Genomic DNA was extracted using the Qiagen DNA Microkit (Qiagen, Hilden, Germany) according to the manufacturer's protocol and quantified using a Qubit Fluorometer 4.0 (Invitrogen, Carlsbad, CA, USA). Whole-genome sequencing was performed at the Noguchi Memorial Institute for Medical Research (Accra, Ghana) on the Illumina MiSeq platform (Illumina, San Diego, CA, USA). Sequencing libraries were prepared from 500 ng starting DNA using the Illumina DNA Prep (M) Tagmentation kit. Library quality was assessed by quantitative PCR (Kapa SYBR Fast qPCR kit; Roche, Basel, Switzerland) and the Agilent 2100 Bioanalyzer (Agilent Technologies, Santa Clara, CA, USA). Libraries were normalised based on fragment sizes and concentrations, pooled, and loaded onto the MiSeq. Raw sequencing reads (FASTQ files) were quality-filtered and trimmed using FastQC v0.11.9 [2] and Trimmomatic v0.39 [3]. High-quality reads were de novo assembled using Unicycler v0.4.9 [4] in normal mode.

## Species confirmation and genomic characterisation

Assembled genomes were analysed in parallel using two complementary platforms. Kleborate v2.4.1 [5] was used for species-level identification within the *K. pneumoniae* species complex, multilocus sequence typing (MLST), detection of acquired antimicrobial resistance genes, virulence determinants, and K and O locus typing via the integrated Kaptive module [6]. BacPipe v1.2 [7], an integrated whole-genome bacterial sequencing analysis pipeline, was used for genome annotation via Prokka [8], resistance gene detection via ResFinder [9] and CARD [10], virulence factor identification via VirulenceFinder [11], plasmid replicon typing via PlasmidFinder [12], resistance gene family classification via Resfams [13], and multilocus sequence typing [14]. Results from Kleborate and BacPipe were cross-validated to ensure concordance in species assignment, MLST, and AMR gene detection. For the investigation of plasmid-mediated resistance genes, annotated reference plasmids from NCBI were additionally mapped to whole-genome sequences using Geneious Prime 2023 [15], and consensus resistance genes were reported.

Species-level assignment identified 78 isolates as *K. pneumoniae* sensu stricto (n = 72) or *K. quasipneumoniae* (n = 6), which were retained for all downstream analyses given the documented clinical relevance and plasmid-carrying capacity of the latter species. Four *K. variicola* isolates were excluded. Forty-nine unique multilocus sequence types were identified across the 78 retained isolates.

## Plasmid reconstruction and characterisation

Plasmid sequences were reconstructed from assembled genomes using MOB-suite v3.1.4 [16]. MOB-recon was used to classify each contig as chromosomal or plasmid-derived based on reference database matching. MOB-typer assigned replicon types (using the PlasmidFinder database), relaxase types, mate-pair formation (MPF) types, predicted mobility class (conjugative, mobilizable, or non-mobilizable), predicted host range at multiple taxonomic levels, and MASH-based nearest-neighbour identification against the MOB-suite reference database. Plasmids were assigned to primary clusters based on MASH distance, enabling tracking of identical or near-identical plasmid backbones across isolates and compartments.

MOB-suite was run separately for clinical and animal isolates (from genome assemblies at the initial sequencing facility) and for environmental isolates (from assemblies processed at a second facility). Results from both runs were compiled into unified master datasets: MASTER_mobtyper.csv (370 plasmid typing records), MASTER_contig_report.csv (29,061 contig classification records), and MASTER_mge_report.csv (mobile genetic element annotations). For clinical and animal isolates, per-isolate MGE reports (mge.report.txt) were additionally compiled from individual MOB-recon output directories, yielding 195 plasmid-associated IS element and transposon records. Combined with the environmental MGE data, this produced a total of 350 plasmid-associated mobile genetic element records across 47 plasmid clusters.

## Plasmid-borne AMR, virulence, and replicon annotation

All reconstructed plasmid sequences from the 78 isolates were concatenated into a single multi-FASTA file (2,600 sequences) with headers formatted as IsolateID|ClusterID|contig_info to enable downstream tracking of each sequence to its source isolate and plasmid cluster. ABRicate v1.0.1 [17] was run on the Galaxy platform (usegalaxy.org) against three databases with default thresholds of ≥80% coverage and ≥90% identity: NCBI AMRFinderPlus for antimicrobial resistance gene detection (347 hits, 56 unique genes), VFDB for virulence factor detection (21 hits), and PlasmidFinder for replicon identification (274 hits, 35 replicon families). AMR genes were classified into 12 drug classes: carbapenemase, ESBL, other β-lactamase, aminoglycoside, quinolone, tetracycline, sulfonamide, trimethoprim, chloramphenicol, rifamycin, macrolide, and other. The VFDB output was specifically examined for canonical hypervirulence markers (iuc, iro, rmpA, rmpA2, clb).

## CRISPR-Cas immune system analysis

CRISPR-Cas system architecture, anti-CRISPR (Acr) protein candidates, and self-targeting spacers were identified using CRISPRimmunity v1.0 [18]. CRISPRimmunity integrates PILER-CR [19], CRISPRCasFinder [20], and CRT [21] for CRISPR array detection; HmmScan against 428 hidden Markov model profiles for Cas protein identification and CRISPR-Cas system classification; and AcRanker [22] for anti-CRISPR candidate prediction based on genomic context, protein features, and homology to known Acr families. CRISPRimmunity was run independently for each of the 78 isolates, and results were compiled from per-isolate output files including acr_result.tab, self_target_result.tab, cas_operons.tab, and merged spacer FASTA files.

Spacer sequences were extracted from merged spacer files for all CRISPR-bearing isolates, yielding 1,207 spacers across 34 isolates (median spacer length 32 bp; range 7-335 bp). Each isolate was classified into one of four mutually exclusive immune categories based on the co-occurrence of Cas proteins, CRISPR arrays, and Acr candidates: (i) CRISPR with anti-CRISPR (neutralised immunity; n = 29); (ii) no CRISPR with Acr present (n = 44); (iii) orphan arrays without Cas proteins (n = 4); and (iv) no CRISPR defence (n = 1).

## Genomic location of anti-CRISPR genes

To determine whether anti-CRISPR genes were chromosomally or plasmid-encoded, the contig number of each Acr candidate was extracted from the CRISPRimmunity acr_candidate_id field and cross-referenced with the MOB-suite contig report (MASTER_contig_report.csv), which classifies every contig in each genome as chromosome or plasmid. This mapping was performed for all 96 Acr genes identified across the dataset with available contig assignments.

## Spacer-plasmid BLAST analysis

To assess whether CRISPR spacers target circulating plasmids, all 1,207 spacer sequences were searched against all 2,600 reconstructed plasmid sequences using NCBI BLAST+ on the Galaxy platform. A nucleotide BLAST database was constructed from the combined plasmid FASTA file using makeblastdb. BLASTn was run using the blastn-short task (optimised for short query sequences), word size 7, e-value cutoff 1, and a maximum of 5 hits per query sequence. Significant hits were defined as ≥95% nucleotide identity over ≥24 bp alignment length (corresponding to approximately 80% of the median spacer length). Near-perfect hits were defined as ≥99% identity over ≥30 bp. Self-targeting events (spacer and plasmid from the same isolate) were distinguished from cross-targeting events (spacer matching a plasmid from a different isolate). Cross-niche targeting was quantified by assigning each spacer and each plasmid to its source compartment (Human, Animal, Environment) and constructing a 3 × 3 targeting matrix of all source-to-source combinations.

## GC deviation analysis

The absolute GC deviation of each plasmid from the *K. pneumoniae* chromosomal mean (57%) was calculated as |plasmid GC% − 57%| and expressed as a percentage. Cross-niche (shared between ≥2 compartments) and niche-specific plasmids were compared within each source compartment using two-sided Wilcoxon rank-sum tests.

## Logistic regression modelling of cross-niche plasmid determinants

To identify plasmid features associated with cross-niche distribution, a logistic regression model was constructed following the analytical framework of Dewan and Uecker [23]. The binary outcome variable was cross-niche status (shared between ≥2 compartments versus niche-specific). Three predictor variables were included: mobility score (conjugative = 3, mobilizable = 2, non-mobilizable = 1), predicted host range breadth (genus = 1 through multi-phylla = 6, assigned from MOB-suite predicted_host_range_overall_rank), and absolute GC deviation from chromosome. The model was fitted to 357 plasmids with complete annotations using the glm function in R with a binomial error distribution and logit link function. Odds ratios and 95% Wald confidence intervals were derived from exponentiated model coefficients. Statistical significance was assessed at α = 0.05.

## Resistome composition analysis

Binary presence-absence matrices of antimicrobial resistance genes were constructed for all 78 isolates from Kleborate output. Jaccard dissimilarity was computed between all isolate pairs using the vegdist function in the vegan R package v2.6-4 [24]. Hierarchical clustering was performed using the Ward.D2 agglomeration method. Principal coordinates analysis (PCoA) was conducted on the Jaccard distance matrix using the cmdscale function. The effect of source compartment (Human, Animal, Environment) on overall resistome composition was tested by permutational multivariate analysis of variance (PERMANOVA) with 999 permutations using the adonis2 function in vegan. Biplot vectors for individual AMR genes were computed and overlaid on the PCoA ordination to identify genes driving compartment-level separation. Differences in AMR gene count per isolate between source compartments were assessed using the Kruskal-Wallis rank-sum test.

## Visualisation and statistical software

All statistical analyses and data visualisations were performed in R v4.5.1 [25] using the following packages: tidyverse v2.0.0 [26] for data manipulation, ggplot2 v4.0.3 [27] for plotting, patchwork v1.2.0 [28] for multi-panel figure assembly, gggenes v0.5.1 [29] for gene arrow diagrams, ggrepel v0.9.5 [30] for non-overlapping labels, circlize v0.4.16 [31] for chord diagrams, vegan v2.6-4 [24] for ecological statistics, VennDiagram v1.7.3 [32] for set overlap visualisation, and readxl v1.4.3 for Excel file import. All figures were generated as vector-format PDF files at 600 dpi resolution and assembled into composite panels in Inkscape v1.3 (https://inkscape.org). Nature-standard figure dimensions were used throughout (single column: 89 mm; double column: 183 mm).

## Data and code availability

Raw sequencing reads have been deposited in the NCBI Sequence Read Archive under BioProject [accession to be assigned]. Assembled genome sequences are available from [repository to be assigned]. All analysis code, including bioinformatics pipelines and R scripts for statistical analysis and figure generation, is available at https://github.com/[username]/kp_prp_analysis. Compiled datasets including MOB-suite results, CRISPRimmunity outputs, ABRicate reports, BLAST results, and the R workspace are deposited at [Zenodo DOI to be assigned].

## References

1. Dong D, et al. *J Antibiot* 67, 215-220 (2014).
2. Andrews S. FastQC (2010). https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
3. Bolger AM, Lohse M, Usadel B. *Bioinformatics* 30, 2114-2120 (2014).
4. Wick RR, Judd LM, Gorrie CL, Holt KE. *PLoS Comput Biol* 13, e1005595 (2017).
5. Lam MMC, et al. *Nat Commun* 12, 4188 (2021).
6. Lam MMC, et al. *Microb Genom* 8, 000800 (2022).
7. Xavier BB, et al. *BMC Genomics* 21, 11 (2020).
8. Seemann T. *Bioinformatics* 30, 2068-2069 (2014).
9. Bortolaia V, et al. *J Antimicrob Chemother* 75, 3491-3500 (2020).
10. Alcock BP, et al. *Nucleic Acids Res* 51, D605-D612 (2023).
11. Liu B, Zheng D, Zhou S, Chen L, Yang J. *Nucleic Acids Res* 50, D912-D917 (2022).
12. Carattoli A, et al. *Antimicrob Agents Chemother* 58, 3895-3903 (2014).
13. Gibson MK, Forsberg KJ, Dantas G. *ISME J* 9, 207-216 (2015).
14. Larsen MV, et al. *J Clin Microbiol* 50, 1355-1361 (2012).
15. Kearse M, et al. *Bioinformatics* 28, 1647-1649 (2012).
16. Robertson J, Nash JHE. *Microb Genom* 4, e000206 (2018).
17. Seemann T. ABRicate (2020). https://github.com/tseemann/abricate
18. Pang Z, et al. *Nucleic Acids Res* 50, W584-W590 (2022).
19. Edgar RC. *BMC Bioinformatics* 8, 18 (2007).
20. Couvin D, et al. *Nucleic Acids Res* 46, W246-W251 (2018).
21. Bland C, et al. *BMC Bioinformatics* 8, 209 (2007).
22. Eitzinger S, et al. *Nucleic Acids Res* 48, 4698-4708 (2020).
23. Dewan I, Uecker H. *Microbiology* 169, 001362 (2023).
24. Oksanen J, et al. vegan v2.6-4 (2022). https://CRAN.R-project.org/package=vegan
25. R Core Team. R v4.5.1 (2024). https://www.R-project.org/
26. Wickham H, et al. *J Open Source Softw* 4, 1686 (2019).
27. Wickham H. ggplot2. Springer-Verlag, New York (2016).
28. Pedersen TL. patchwork (2022). https://CRAN.R-project.org/package=patchwork
29. Wilkins D. gggenes (2023). https://CRAN.R-project.org/package=gggenes
30. Slowikowski K. ggrepel (2023). https://CRAN.R-project.org/package=ggrepel
31. Gu Z, Gu L, Eils R, Schlesner M, Brors B. *Bioinformatics* 30, 2811-2812 (2014).
32. Chen H, Boutros PC. *BMC Bioinformatics* 12, 35 (2011).
