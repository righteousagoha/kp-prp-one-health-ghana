# ═══════════════════════════════════════════════════════════════
# 07_resistome_analysis.R
# Resistome clustering, PCoA, PERMANOVA (Figure 2)
# ═══════════════════════════════════════════════════════════════

library(tidyverse)
library(vegan)
load("results/workspace_full.RData")

# Binary presence-absence matrix from Kleborate
amr_genes <- kleb %>%
  select(strain, matches("^bla|^qnr|^aac|^aph|^aad|^sul|^tet|^dfr|^cat|^fos|^arr|^mph|^oqx|^mcr|^arm|^flo|^msr|^erm")) %>%
  mutate(across(-strain, ~ifelse(. == "-" | is.na(.), 0, 1)))

# Jaccard distance
jac_dist <- vegdist(amr_genes %>% select(-strain) %>% as.matrix(), method = "jaccard")

# Hierarchical clustering
hc <- hclust(jac_dist, method = "ward.D2")

# PCoA
pcoa_res <- cmdscale(jac_dist, k = 2, eig = TRUE)
var_explained <- round(pcoa_res$eig[1:2] / sum(pcoa_res$eig) * 100, 1)

# PERMANOVA
source_vec <- kleb$source[match(amr_genes$strain, kleb$strain)]
perm_result <- adonis2(jac_dist ~ source_vec, permutations = 999)
cat("PERMANOVA p =", perm_result$`Pr(>F)`[1], "\n")
cat("R² =", round(perm_result$R2[1], 3), "\n")

# ═══════════════════════════════════════════════════════════════
# 08_plasmid_sharing.R
# Cross-niche analysis, GLM, GC deviation (Figures 3-4)
# ═══════════════════════════════════════════════════════════════

# GC deviation comparison
gc_test <- mob %>%
  filter(primary_cluster_id != "-") %>%
  mutate(is_shared = primary_cluster_id %in% shared_clusters)

# Wilcoxon tests by source
for (src in c("Human", "Animal", "Environment")) {
  test <- wilcox.test(
    gc_deviation ~ is_shared,
    data = gc_test %>% filter(source == src)
  )
  cat(src, "GC deviation p =", test$p.value, "\n")
}

# Logistic regression model (Dewan & Uecker framework)
glm_data <- mob %>%
  filter(primary_cluster_id != "-") %>%
  mutate(
    is_shared = as.integer(primary_cluster_id %in% shared_clusters),
    mobility_score = case_when(
      predicted_mobility == "conjugative" ~ 3,
      predicted_mobility == "mobilizable" ~ 2,
      predicted_mobility == "non-mobilizable" ~ 1
    ),
    host_breadth = case_when(
      predicted_host_range_overall_rank == "genus" ~ 1,
      predicted_host_range_overall_rank == "family" ~ 2,
      predicted_host_range_overall_rank == "order" ~ 3,
      predicted_host_range_overall_rank == "class" ~ 4,
      predicted_host_range_overall_rank == "phylum" ~ 5,
      predicted_host_range_overall_rank == "multi-phylla" ~ 6,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(mobility_score), !is.na(host_breadth))

model <- glm(is_shared ~ mobility_score + host_breadth + gc_deviation,
             data = glm_data, family = binomial)

# Extract odds ratios
or_table <- exp(cbind(OR = coef(model), confint(model)))
cat("\n=== GLM ODDS RATIOS ===\n")
print(or_table)
cat("P-values:\n")
print(summary(model)$coefficients[, "Pr(>|z|)"])

# ═══════════════════════════════════════════════════════════════
# 09_crispr_visualisation.R
# CRISPR schematic with gggenes (Figure 5i)
# ═══════════════════════════════════════════════════════════════

library(gggenes)

# Canonical I-E operon order
ie_order <- c("DEDDh", "DinG", "RT", "TnsC", "c2c9_V-U4", "csa3",
              "cas3", "cas8e", "cse2gr11", "cas7", "cas5", "csf2gr7", "cas6e")
ie_post_array <- c("cas1", "cas2")

# Build gene features with correct operon architecture
# [See 06_data_integration.R for gene_features3 construction]

# Muted Nature palette
col_cas <- "#4878A6"
col_array <- "#45A29E"
col_orphan <- "#B0B0B0"
col_acr <- "#C44E52"
col_self <- "#DD8452"
col_backbone <- "#E8E8E8"

# ═══════════════════════════════════════════════════════════════
# 10_spacer_plasmid_analysis.R
# BLAST parsing, cross-niche targeting (Figure 5ii)
# ═══════════════════════════════════════════════════════════════

# Cross-niche targeting matrix
cross_niche_hits <- sig_hits %>%
  filter(cross_hit) %>%
  mutate(
    spacer_source = case_when(
      str_detect(spacer_isolate, "^KATH|^KBTH|^HTH") ~ "Human",
      str_detect(spacer_isolate, "^PF|^CTL|^GT") ~ "Animal",
      str_detect(spacer_isolate, "^WT|^UGK") ~ "Environment"
    ),
    plasmid_source = case_when(
      str_detect(plasmid_isolate, "^KATH|^KBTH|^HTH") ~ "Human",
      str_detect(plasmid_isolate, "^PF|^CTL|^GT") ~ "Animal",
      str_detect(plasmid_isolate, "^WT|^UGK") ~ "Environment"
    )
  ) %>%
  filter(!is.na(spacer_source), !is.na(plasmid_source))

cat("=== CROSS-NICHE TARGETING MATRIX ===\n")
cross_niche_hits %>% count(spacer_source, plasmid_source) %>% print()

# Acr genomic location analysis
acr_all <- acr_all %>%
  mutate(contig_num = str_extract(acr_candidate_id, "^\\d+"))

contig <- read_csv("results/mob_suite/MASTER_contig_report.csv", show_col_types = FALSE) %>%
  mutate(contig_num = str_extract(contig_id, "^\\d+"))

# Map each Acr to chromosome vs plasmid
acr_location <- map_dfr(1:nrow(acr_all), function(i) {
  s <- acr_all$Sample[i]
  cn <- acr_all$contig_num[i]
  match <- contig %>%
    filter(str_detect(sample_id, fixed(s)), contig_num == cn) %>%
    select(sample_id, molecule_type, primary_cluster_id, contig_num)
  if (nrow(match) > 0) {
    match %>% mutate(Sample = s, acr_id = acr_all$acr_candidate_id[i])
  } else {
    tibble(Sample = s, molecule_type = "not_found", contig_num = cn,
           acr_id = acr_all$acr_candidate_id[i])
  }
})

cat("\n=== ACR GENOMIC LOCATION ===\n")
table(acr_location$molecule_type)

# ═══════════════════════════════════════════════════════════════
# 11_prp_characterisation.R
# AMR cargo, virulence, IS elements on plasmids (Figure 6)
# ═══════════════════════════════════════════════════════════════

# AMR genes per plasmid cluster
amr_by_cluster <- amr %>%
  filter(!is.na(cluster_id), cluster_id != "") %>%
  group_by(cluster_id) %>%
  summarise(
    n_amr_genes = n_distinct(GENE),
    amr_genes = paste(sort(unique(GENE)), collapse = ", "),
    resistance_classes = paste(sort(unique(RESISTANCE)), collapse = "; "),
    STs = paste(sort(unique(ST)), collapse = ", "),
    sources = paste(sort(unique(source)), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_amr_genes))

# Virulence genes
cat("\n=== HYPERVIRULENCE MARKERS ===\n")
cat("iuc:", vf %>% filter(str_detect(GENE, "iuc|iut")) %>% nrow(), "\n")
cat("iro:", vf %>% filter(str_detect(GENE, "iro")) %>% nrow(), "\n")
cat("rmpA:", vf %>% filter(str_detect(GENE, "rmp")) %>% nrow(), "\n")
cat("clb:", vf %>% filter(str_detect(GENE, "clb")) %>% nrow(), "\n")
cat("mrk (biofilm):", vf %>% filter(str_detect(GENE, "mrk")) %>% nrow(), "\n")

# Cross-niche AMR enrichment
n_shared_amr <- sum(shared_clusters %in% amr_by_cluster$cluster_id)
n_shared_total <- length(shared_clusters)
n_specific_amr <- sum(!(amr_by_cluster$cluster_id %in% shared_clusters))
n_specific_total <- n_distinct(mob$primary_cluster_id[mob$primary_cluster_id != "-"]) - n_shared_total

fisher_mat <- matrix(c(n_shared_amr, n_shared_total - n_shared_amr,
                        n_specific_amr, n_specific_total - n_specific_amr), nrow = 2)
fisher_test <- fisher.test(fisher_mat)
cat("\nFisher's test for AMR enrichment on cross-niche clusters:\n")
cat("p =", fisher_test$p.value, "\n")
cat("OR =", fisher_test$estimate, "\n")

# MGE compilation
mge_env <- read_csv("results/mob_suite/MASTER_mge_report.csv", show_col_types = FALSE)

# Compile clinical/animal MGE from per-isolate reports
base_mge <- "results/mob_suite/mob_recon"
folders <- list.files(base_mge, full.names = FALSE)

mge_clin <- map_dfr(folders, function(f) {
  mge_file <- file.path(base_mge, f, "mge.report.txt")
  if (file.exists(mge_file)) {
    tryCatch({
      df <- read_tsv(mge_file, show_col_types = FALSE, col_types = cols(.default = "c"))
      if (nrow(df) > 0) {
        df$isolate_id <- str_remove(f, "_S\\d+_L\\d+$")
        return(df)
      }
    }, error = function(e) NULL)
  }
  return(NULL)
})

cat("\nMGE records (clinical/animal):", nrow(mge_clin), "\n")
cat("MGE records (environmental):", nrow(mge_env %>% filter(molecule_type == "plasmid")), "\n")

save.image("results/workspace_full.RData")
cat("\nAll analyses complete. Workspace saved.\n")
