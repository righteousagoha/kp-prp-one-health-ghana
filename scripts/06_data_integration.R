# ═══════════════════════════════════════════════════════════════
# 06_data_integration.R
# Data loading, merging, quality control, and sample mapping
# ═══════════════════════════════════════════════════════════════

library(tidyverse)
library(readxl)

# ── 1. Load MOB-suite data ──
mob <- read_csv("results/mob_suite/MASTER_mobtyper.csv", show_col_types = FALSE)

# Clean isolate IDs and assign sources
mob <- mob %>%
  mutate(
    isolate_id = str_remove(sample_id, "_S\\d+_L\\d+$"),
    source = case_when(
      str_detect(isolate_id, "^KATH|^KBTH|^HTH") ~ "Human",
      str_detect(isolate_id, "^PF|^CTL|^GT") ~ "Animal",
      str_detect(isolate_id, "^WT|^UGK") ~ "Environment"
    ),
    city = case_when(
      str_detect(isolate_id, "^KATH") ~ "Kumasi",
      str_detect(isolate_id, "^KBTH|^UGK|^WT|^PF|^CTL|^GT") ~ "Accra",
      str_detect(isolate_id, "^HTH") ~ "Ho"
    )
  )

# ── 2. Load Kleborate data ──
kleb <- read_xlsx("results/kleborate/Kleborate_Results.xlsx")

# Map ST to MOB data
mob <- mob %>%
  left_join(kleb %>% select(strain, ST) %>%
              mutate(isolate_id = str_remove(strain, "_S\\d+_L\\d+$")) %>%
              distinct(isolate_id, ST),
            by = "isolate_id")

# ── 3. GC deviation from chromosomal mean ──
mob <- mob %>%
  mutate(gc_deviation = abs(gc - 0.57) * 100)

# ── 4. Define cross-niche clusters ──
cluster_sharing <- mob %>%
  filter(primary_cluster_id != "-", !is.na(primary_cluster_id)) %>%
  group_by(primary_cluster_id) %>%
  summarise(
    n_sources = n_distinct(source),
    sources = paste(sort(unique(source)), collapse = ", "),
    n_plasmids = n(),
    STs = paste(sort(unique(ST)), collapse = ", "),
    rep = first(`rep_type(s)`),
    mobility = paste(unique(predicted_mobility), collapse = "/"),
    .groups = "drop"
  )

shared_clusters <- cluster_sharing %>%
  filter(n_sources >= 2) %>%
  pull(primary_cluster_id)

cat("Total plasmids:", nrow(mob), "\n")
cat("Unique clusters:", n_distinct(mob$primary_cluster_id), "\n")
cat("Cross-niche clusters:", length(shared_clusters), "\n")

# ── 5. Load CRISPR data ──
# Compile per-isolate CRISPR results
samples <- list.dirs("results/crispr/per_isolate", recursive = FALSE, full.names = FALSE)

# Load Acr results
acr_all <- map_dfr(samples, function(s) {
  f <- file.path("results/crispr/per_isolate", s, "acr_result.tab")
  if (file.exists(f) && file.size(f) > 0) {
    tryCatch({
      df <- read_tsv(f, show_col_types = FALSE)
      if (nrow(df) > 0) { df$Sample <- s; return(df) }
    }, error = function(e) NULL)
  }
  return(NULL)
})

# Load self-targeting results
st_all <- map_dfr(samples, function(s) {
  f <- file.path("results/crispr/per_isolate", s, "self_target_result.tab")
  if (file.exists(f) && file.size(f) > 0) {
    tryCatch({
      df <- read_tsv(f, show_col_types = FALSE)
      if (nrow(df) > 0) { df$Sample <- s; return(df) }
    }, error = function(e) NULL)
  }
  return(NULL)
})

# Load spacer sequences
spacer_all <- map_dfr(samples, function(s) {
  f <- list.files(file.path("results/crispr/per_isolate", s),
                  pattern = "_merge\\.spc$", full.names = TRUE)
  if (length(f) > 0 && file.exists(f[1])) {
    lines <- readLines(f[1])
    if (length(lines) > 1) {
      headers <- which(str_detect(lines, "^>"))
      return(tibble(Sample = s, header = lines[headers], spacer_seq = lines[headers + 1]))
    }
  }
  return(NULL)
})

# ── 6. Build merged CRISPR summary ──
crispr_merged <- tibble(Sample = unique(c(
  mob %>% pull(isolate_id) %>% unique(),
  samples
))) %>%
  left_join(kleb %>%
              mutate(Sample = str_remove(strain, "_S\\d+_L\\d+$")) %>%
              distinct(Sample, ST),
            by = "Sample") %>%
  mutate(
    source = case_when(
      str_detect(Sample, "^KATH|^KBTH|^HTH") ~ "Human",
      str_detect(Sample, "^PF|^CTL|^GT") ~ "Animal",
      str_detect(Sample, "^WT|^UGK") ~ "Environment"
    )
  ) %>%
  left_join(acr_all %>% count(Sample, name = "Acr_Candidates"), by = "Sample") %>%
  left_join(st_all %>% count(Sample, name = "Self_Targeting_Spacers"), by = "Sample") %>%
  left_join(spacer_all %>% count(Sample, name = "n_spacers"), by = "Sample") %>%
  mutate(across(where(is.numeric), ~replace_na(., 0)))

# Classify immune status
crispr_merged <- crispr_merged %>%
  mutate(
    immune_status = case_when(
      CRISPR_Arrays_CI > 0 & Acr_Candidates > 0 ~ "CRISPR + Anti-CRISPR (neutralised)",
      CRISPR_Arrays_CI == 0 & Orphan_Arrays > 0 ~ "Orphan arrays (no Cas)",
      CRISPR_Arrays_CI == 0 & Acr_Candidates > 0 ~ "No CRISPR + Acr",
      TRUE ~ "No CRISPR defence"
    )
  )

cat("\n=== IMMUNE STATUS ===\n")
table(crispr_merged$immune_status)

# ── 7. Load ABRicate results ──
amr <- read_tsv("results/abricate/abricate_ncbi_amr.tab", show_col_types = FALSE) %>%
  mutate(
    isolate_id = str_remove(str_extract(SEQUENCE, "^[^|]+"), "-21$|-22$"),
    cluster_id = str_extract(SEQUENCE, "(?<=\\|)[^|]+"),
    source = case_when(
      str_detect(isolate_id, "^KATH|^KBTH|^HTH") ~ "Human",
      str_detect(isolate_id, "^PF|^CTL|^GT") ~ "Animal",
      str_detect(isolate_id, "^WT|^UGK") ~ "Environment"
    )
  ) %>%
  left_join(crispr_merged %>% select(Sample, ST) %>%
              rename(isolate_id = Sample), by = "isolate_id")

vf <- read_tsv("results/abricate/abricate_vfdb.tab", show_col_types = FALSE) %>%
  mutate(
    isolate_id = str_remove(str_extract(SEQUENCE, "^[^|]+"), "-21$|-22$"),
    cluster_id = str_extract(SEQUENCE, "(?<=\\|)[^|]+")
  )

pf <- read_tsv("results/abricate/abricate_plasmidfinder.tab", show_col_types = FALSE) %>%
  mutate(
    isolate_id = str_remove(str_extract(SEQUENCE, "^[^|]+"), "-21$|-22$"),
    cluster_id = str_extract(SEQUENCE, "(?<=\\|)[^|]+")
  )

# ── 8. Load BLAST results ──
blast <- read_tsv("results/blast/spacer_vs_plasmid.tab",
                   col_names = c("qseqid", "sseqid", "pident", "length", "mismatch",
                                 "gapopen", "qstart", "qend", "sstart", "send",
                                 "evalue", "bitscore"),
                   show_col_types = FALSE)

sig_hits <- blast %>%
  filter(pident >= 95, length >= 24) %>%
  mutate(
    spacer_isolate = str_extract(qseqid, "^[^|]+"),
    plasmid_isolate = str_extract(sseqid, "^[^|]+"),
    plasmid_cluster = str_extract(sseqid, "(?<=\\|)[^|]+"),
    self_hit = spacer_isolate == plasmid_isolate,
    cross_hit = spacer_isolate != plasmid_isolate
  )

cat("\n=== BLAST RESULTS ===\n")
cat("Total hits:", nrow(blast), "\n")
cat("Significant:", nrow(sig_hits), "\n")

# ── 9. Save workspace ──
save.image("results/workspace_full.RData")
cat("\nData integration complete. Workspace saved.\n")
