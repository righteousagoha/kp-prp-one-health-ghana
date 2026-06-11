# ═══════════════════════════════════════════════════════════════
# 08_acr_genomic_location.R
# Anti-CRISPR genomic location analysis
# Determines whether Acr genes are chromosomal or plasmid-borne
# ═══════════════════════════════════════════════════════════════

library(tidyverse)
library(patchwork)

load("results/workspace_full.RData")

# ── 1. Extract contig number from Acr candidate ID ──
acr_all <- acr_all %>%
  mutate(contig_num = str_extract(acr_candidate_id, "^\\d+"))

# ── 2. Load contig report ──
contig <- read_csv("results/mob_suite/MASTER_contig_report.csv",
                    show_col_types = FALSE) %>%
  mutate(contig_num = str_extract(contig_id, "^\\d+"))

# ── 3. Map each Acr gene to chromosome vs plasmid ──
acr_with_location <- map_dfr(1:nrow(acr_all), function(i) {
  s <- acr_all$Sample[i]
  cn <- acr_all$contig_num[i]
  match <- contig %>%
    filter(str_detect(sample_id, fixed(s)), contig_num == cn) %>%
    select(sample_id, molecule_type, primary_cluster_id, contig_num) %>%
    slice(1)
  if (nrow(match) > 0) {
    match %>% mutate(
      Sample = s,
      acr_id = acr_all$acr_candidate_id[i],
      in_prophage = acr_all$Acr_candidate_in_prophage[i],
      acr_homology = acr_all$homology_with_known_acr[i]
    )
  } else {
    tibble(Sample = s, molecule_type = "not_found", contig_num = cn,
           acr_id = acr_all$acr_candidate_id[i],
           in_prophage = acr_all$Acr_candidate_in_prophage[i],
           acr_homology = acr_all$homology_with_known_acr[i])
  }
})

# ── 4. Summary statistics ──
cat("=== ACR GENOMIC LOCATION ===\n")
table(acr_with_location$molecule_type)

# ── 5. Add metadata ──
acr_with_location <- acr_with_location %>%
  left_join(crispr_merged %>% select(Sample, ST, source), by = "Sample") %>%
  mutate(
    location = case_when(
      molecule_type == "chromosome" ~ "Chromosome",
      molecule_type == "plasmid" ~ "Plasmid",
      TRUE ~ "Unassigned"
    ),
    acr_family = case_when(
      str_detect(acr_homology, "AcrIIC1") ~ "AcrIIC1",
      str_detect(acr_homology, "AcrIF11") ~ "AcrIF11",
      str_detect(acr_homology, "AcrIIA7") ~ "AcrIIA7",
      TRUE ~ "Unknown"
    ),
    source = factor(source, levels = c("Human", "Animal", "Environment"))
  )

cat("\n=== ACR FAMILIES BY LOCATION ===\n")
acr_with_location %>%
  filter(location != "Unassigned") %>%
  count(acr_family, location) %>%
  print()

cat("\n=== PROPHAGE STATUS ===\n")
table(acr_with_location$in_prophage)

cat("\n=== HTH PROTEIN ASSOCIATIONS ===\n")
table(!is.na(acr_all$HTH_protein) & acr_all$HTH_protein != "")

cat("\n=== ACA PROTEIN ASSOCIATIONS ===\n")
table(!is.na(acr_all$neighbor_aca_protein))

# ── 6. Plasmid-borne Acr details ──
cat("\n=== PLASMID-BORNE ACR DETAILS ===\n")
acr_with_location %>%
  filter(location == "Plasmid") %>%
  select(Sample, ST, source, primary_cluster_id, acr_family) %>%
  print()

# ── 7. Visualisation: Supplementary Figure ──

# Panel A: Acr family by location
pA <- acr_with_location %>%
  filter(location != "Unassigned") %>%
  count(acr_family, location) %>%
  ggplot(aes(x = acr_family, y = n, fill = location)) +
  geom_col(position = "stack", width = 0.5, colour = "grey30", linewidth = 0.3) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5),
            size = 2.5, colour = "white", fontface = "bold") +
  scale_fill_manual(values = c("Chromosome" = "#4878A6", "Plasmid" = "#CC3311"),
                    name = "Location") +
  labs(x = NULL, y = "Number of Acr genes",
       title = "Anti-CRISPR protein families and genomic location") +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(face = "italic"),
        plot.title = element_text(size = 9, face = "bold"),
        axis.line = element_line(linewidth = 0.3))

# Panel B: Acr by source and family
pB <- acr_with_location %>%
  filter(location != "Unassigned") %>%
  count(source, acr_family) %>%
  ggplot(aes(x = source, y = n, fill = acr_family)) +
  geom_col(position = "stack", width = 0.5, colour = "grey30", linewidth = 0.3) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5),
            size = 2.5, colour = "white", fontface = "bold") +
  scale_fill_manual(values = c("AcrIIC1" = "#4878A6", "AcrIF11" = "#009988",
                                "AcrIIA7" = "#882255"), name = "Acr family") +
  labs(x = NULL, y = "Number of Acr genes") +
  theme_classic(base_size = 8) +
  theme(legend.text = element_text(face = "italic"),
        axis.line = element_line(linewidth = 0.3))

# Panel C: Acr per ST tile
acr_st <- acr_with_location %>%
  filter(location != "Unassigned") %>%
  distinct(ST, source, acr_family, location) %>%
  mutate(is_cross = ST %in% c("ST101", "ST147", "ST15", "ST17", "ST39", "ST1427"))

st_ord <- acr_st %>%
  distinct(ST, source, is_cross) %>%
  arrange(source, desc(is_cross), ST) %>%
  pull(ST) %>% unique()

acr_st$ST <- factor(acr_st$ST, levels = rev(st_ord))

pC <- ggplot(acr_st, aes(x = acr_family, y = ST, fill = location)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("Chromosome" = "#4878A6", "Plasmid" = "#CC3311"),
                    name = "Location") +
  facet_grid(source ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 7) +
  theme(panel.grid = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = 7, face = "bold"),
        strip.placement = "outside",
        axis.text.x = element_text(size = 6, face = "italic"),
        axis.text.y = element_text(size = 4.5))

# Combine
fig_acr <- (pA | pB) / pC +
  plot_layout(heights = c(1, 2)) +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.tag = element_text(size = 10, face = "bold")))

ggsave("figures/fig_acr_context.pdf", fig_acr,
       width = 183, height = 200, units = "mm", dpi = 600)
ggsave("figures/fig_acr_context.png", fig_acr,
       width = 183, height = 200, units = "mm", dpi = 600)

# ── 8. Export ──
write_csv(acr_with_location, "results/acr_genomic_location.csv")
save.image("results/workspace_full.RData")
cat("\nAcr genomic location analysis complete.\n")
