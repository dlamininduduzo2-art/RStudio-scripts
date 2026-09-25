# ==============================================================================
# SILVA 144 MICROBIOME COMPOSITION AND ASSOCIATION NETWORK ANALYSIS
# ==============================================================================
#
# Potato rhizosphere bacterial microbiome
#
# INPUT FILES
#   R_taxonomy_plots/feature-table.tsv
#   R_taxonomy_plots/taxonomy144/taxonomy.tsv
#
# OUTPUT DIRECTORY
#   R_taxonomy_plots/final_outputs/
#
# FIGURES
#   (A) Phylum-level relative abundance
#   (B) Genus-level relative abundance
#   (C) Top 20 genera heatmap
#   (D) Global genus-level association network
#   (E) Most connected genera in the global network
#   (F) Clanwilliam association subnetwork
#   (G) Dendron association subnetwork
#   (H) Mamusha association subnetwork
#   (I) Wesselesbron association subnetwork
#
# IMPORTANT
#   Network edges represent statistical associations, not proven direct
#   ecological interactions.
#
#   Farm subnetworks inherit statistically supported edges from the GLOBAL
#   20-sample network. Correlations are NOT re-estimated from only five
#   samples per farm.
# ==============================================================================


# ==============================================================================
# 1. SETTINGS
# ==============================================================================

FEATURE_FILE <- "R_taxonomy_plots/feature-table.tsv"

TAXONOMY_FILE <- "R_taxonomy_plots/taxonomy144/taxonomy.tsv"

OUTPUT_DIR <- "R_taxonomy_plots/final_outputs"

FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")

TABLE_DIR <- file.path(OUTPUT_DIR, "tables")

NETWORK_DIR <- file.path(OUTPUT_DIR, "networks")

FARM_NETWORK_DIR <- file.path(
  NETWORK_DIR,
  "farm_subnetworks"
)


# Taxonomic figures
TOP_PHYLA <- 10
TOP_GENERA <- 15
TOP_HEATMAP_GENERA <- 20


# Global network thresholds
PREVALENCE_THRESHOLD <- 0.20
CORRELATION_THRESHOLD <- 0.60
FDR_THRESHOLD <- 0.05
PSEUDOCOUNT <- 0.5
CORRELATION_METHOD <- "spearman"


# Farm subnetworks
# Genus must occur in at least 2 of the farm samples
FARM_MIN_SAMPLES <- 2


# Remove chloroplast / mitochondria / obvious non-bacterial sequences
REMOVE_NON_TARGET <- TRUE


# Farm order
FARM_ORDER <- c(
  "Clanwilliam",
  "Dendron",
  "Mamusha",
  "Wesselesbron"
)


# ==============================================================================
# 2. STANDARD FIGURE LABELS
# ==============================================================================

FIGURE_TITLES <- c(
  A = "(A) Phylum-level relative abundance",
  B = "(B) Genus-level relative abundance",
  C = "(C) Top 20 genera heatmap",
  D = "(D) Global genus-level association network",
  E = "(E) Most connected genera in the global network",
  F = "(F) Clanwilliam association subnetwork",
  G = "(G) Dendron association subnetwork",
  H = "(H) Mamusha association subnetwork",
  I = "(I) Wesselesbron association subnetwork"
)


FARM_FIGURE_CODES <- c(
  Clanwilliam = "F",
  Dendron = "G",
  Mamusha = "H",
  Wesselesbron = "I"
)


# ==============================================================================
# 3. COLOURS
# ==============================================================================

# Network edge colours
EDGE_POSITIVE <- "#0072B2"
EDGE_NEGATIVE <- "#D55E00"


# Neutral taxonomy colours
COLOUR_OTHER <- "#E6E6E6"
COLOUR_UNCLASSIFIED <- "#8C8C8C"


# Farm colours
FARM_COLOURS <- c(
  Clanwilliam = "#0072B2",
  Dendron = "#009E73",
  Mamusha = "#D55E00",
  Wesselesbron = "#CC79A7"
)


# ==============================================================================
# 4. PACKAGES
# ==============================================================================

required_packages <- c(
  "tidyverse",
  "igraph",
  "scales",
  "viridisLite"
)


missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]


if (length(missing_packages) > 0) {
  stop(
    paste(
      "Install these packages first:",
      paste(missing_packages, collapse = ", ")
    )
  )
}


library(tidyverse)
library(igraph)
library(scales)
library(viridisLite)


# ==============================================================================
# 5. OUTPUT DIRECTORIES
# ==============================================================================

dir.create(
  FIGURE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  TABLE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  NETWORK_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  FARM_NETWORK_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 6. STANDARD GGPLOT THEMES
# ==============================================================================

theme_microbiome <- function() {

  theme_classic(base_size = 12) +

    theme(
      plot.title = element_text(
        face = "bold",
        size = 14,
        hjust = 0
      ),

      plot.subtitle = element_text(
        size = 10.5,
        colour = "grey30",
        margin = margin(b = 10)
      ),

      axis.title = element_text(
        face = "bold"
      ),

      axis.text = element_text(
        colour = "black"
      ),

      strip.background = element_rect(
        fill = "#F2F2F2",
        colour = NA
      ),

      strip.text = element_text(
        face = "bold",
        colour = "black"
      ),

      legend.title = element_text(
        face = "bold"
      ),

      legend.position = "right"
    )
}


theme_heatmap_clean <- function() {

  theme_minimal(base_size = 12) +

    theme(
      plot.title = element_text(
        face = "bold",
        size = 14,
        hjust = 0
      ),

      plot.subtitle = element_text(
        size = 10.5,
        colour = "grey30",
        margin = margin(b = 10)
      ),

      panel.grid = element_blank(),

      strip.background = element_rect(
        fill = "#F2F2F2",
        colour = NA
      ),

      strip.text = element_text(
        face = "bold"
      ),

      axis.text = element_text(
        colour = "black"
      ),

      legend.title = element_text(
        face = "bold"
      )
    )
}


# ==============================================================================
# 7. HELPER FUNCTIONS
# ==============================================================================


# ------------------------------------------------------------------------------
# Assign farm from sample ID
# ------------------------------------------------------------------------------

farm_from_sample <- function(x) {

  case_when(
    str_detect(x, "^ND-CL") ~ "Clanwilliam",
    str_detect(x, "^ND-DEN") ~ "Dendron",
    str_detect(x, "^ND-MA") ~ "Mamusha",
    str_detect(x, "^ND-WB") ~ "Wesselesbron",
    TRUE ~ NA_character_
  )
}


# ------------------------------------------------------------------------------
# Extract SILVA taxonomy rank
#
# Supports:
# d__Bacteria; p__Pseudomonadota...
#
# and:
# D_0__Bacteria; D_1__Proteobacteria...
# ------------------------------------------------------------------------------

extract_rank <- function(x, modern_prefix, old_rank) {

  modern_pattern <- paste0(
    "(?:^|;\\s*)",
    modern_prefix,
    "__([^;]+)"
  )

  old_pattern <- paste0(
    "(?:^|;\\s*)D_",
    old_rank,
    "__([^;]+)"
  )


  modern_match <- str_match(
    x,
    modern_pattern
  )[, 2]


  old_match <- str_match(
    x,
    old_pattern
  )[, 2]


  result <- ifelse(
    !is.na(modern_match),
    modern_match,
    old_match
  )


  result <- str_trim(result)


  result[
    result == "" |
      result == "Unassigned"
  ] <- NA


  result
}


# ------------------------------------------------------------------------------
# Taxonomic colour palette
# ------------------------------------------------------------------------------

make_taxon_palette <- function(categories) {

  named_taxa <- setdiff(
    categories,
    c(
      "Other",
      "Unclassified"
    )
  )


  colours <- setNames(
    hcl.colors(
      length(named_taxa),
      palette = "Dynamic"
    ),
    named_taxa
  )


  if ("Unclassified" %in% categories) {

    colours <- c(
      colours,
      Unclassified = COLOUR_UNCLASSIFIED
    )
  }


  if ("Other" %in% categories) {

    colours <- c(
      colours,
      Other = COLOUR_OTHER
    )
  }


  colours
}


# ------------------------------------------------------------------------------
# Safe numeric rescaling
# ------------------------------------------------------------------------------

safe_rescale <- function(
    x,
    to = c(1, 10),
    constant_value = mean(to)
) {

  if (length(x) == 0) {
    return(numeric(0))
  }


  if (length(unique(x)) <= 1) {

    return(
      rep(
        constant_value,
        length(x)
      )
    )
  }


  scales::rescale(
    x,
    to = to
  )
}


# ------------------------------------------------------------------------------
# Save ggplot as both PNG and PDF
# ------------------------------------------------------------------------------

save_figure <- function(
    plot_object,
    file_stub,
    width,
    height
) {

  ggsave(
    filename = file.path(
      FIGURE_DIR,
      paste0(file_stub, ".png")
    ),
    plot = plot_object,
    width = width,
    height = height,
    dpi = 600
  )


  ggsave(
    filename = file.path(
      FIGURE_DIR,
      paste0(file_stub, ".pdf")
    ),
    plot = plot_object,
    width = width,
    height = height
  )
}


# ==============================================================================
# 8. IMPORT FEATURE TABLE
# ==============================================================================

asv <- read.delim(
  FEATURE_FILE,
  skip = 1,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)


colnames(asv)[1] <- "FeatureID"


sample_order <- colnames(asv)[-1]


sample_metadata <- tibble(
  Sample = sample_order
) %>%

  mutate(
    Farm = farm_from_sample(Sample),

    Farm = factor(
      Farm,
      levels = FARM_ORDER
    ),

    Sample = factor(
      Sample,
      levels = sample_order
    )
  )


cat("\n========================================\n")
cat("SAMPLE METADATA\n")
cat("========================================\n")

print(
  sample_metadata,
  n = Inf
)


if (any(is.na(sample_metadata$Farm))) {

  warning(
    "At least one sample could not be assigned to a farm."
  )
}


# ==============================================================================
# 9. IMPORT SILVA TAXONOMY
# ==============================================================================

tax <- read.delim(
  TAXONOMY_FILE,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)


colnames(tax)[1] <- "FeatureID"


if (!"Taxon" %in% colnames(tax)) {

  stop(
    "No column called 'Taxon' was found in taxonomy.tsv."
  )
}


tax_clean <- tax %>%

  transmute(
    FeatureID = FeatureID,

    FullTaxonomy = Taxon,

    Domain = extract_rank(
      Taxon,
      "(?:d|k)",
      0
    ),

    Phylum = extract_rank(
      Taxon,
      "p",
      1
    ),

    Class = extract_rank(
      Taxon,
      "c",
      2
    ),

    Order = extract_rank(
      Taxon,
      "o",
      3
    ),

    Family = extract_rank(
      Taxon,
      "f",
      4
    ),

    Genus = extract_rank(
      Taxon,
      "g",
      5
    ),

    Species = extract_rank(
      Taxon,
      "s",
      6
    )
  )


write.csv(
  tax_clean,
  file.path(
    TABLE_DIR,
    "SILVA144_cleaned_taxonomy.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 10. TAXONOMIC RESOLUTION
# ==============================================================================

taxonomy_resolution <- tax_clean %>%

  summarise(
    Total_ASVs = n(),

    Phylum_resolved = sum(
      !is.na(Phylum)
    ),

    Family_resolved = sum(
      !is.na(Family)
    ),

    Genus_resolved = sum(
      !is.na(Genus)
    ),

    Species_resolved = sum(
      !is.na(Species)
    )
  )


write.csv(
  taxonomy_resolution,
  file.path(
    TABLE_DIR,
    "taxonomy_resolution_summary.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 11. FEATURE TABLE TO LONG FORMAT
# ==============================================================================

asv_long <- asv %>%

  pivot_longer(
    cols = -FeatureID,
    names_to = "Sample",
    values_to = "Count"
  ) %>%

  mutate(
    Count = as.numeric(Count),
    Sample = as.character(Sample)
  )


dat <- asv_long %>%

  left_join(
    tax_clean,
    by = "FeatureID"
  ) %>%

  left_join(
    sample_metadata %>%
      mutate(
        Sample = as.character(Sample)
      ),
    by = "Sample"
  )


# ==============================================================================
# 12. REMOVE NON-TARGET FEATURES
# ==============================================================================

if (REMOVE_NON_TARGET) {

  dat <- dat %>%

    filter(
      !str_detect(
        coalesce(
          FullTaxonomy,
          ""
        ),
        regex(
          "chloroplast|mitochondria",
          ignore_case = TRUE
        )
      )
    ) %>%

    filter(
      is.na(Domain) |
        Domain == "Bacteria"
    )
}


# ==============================================================================
# 13. SAMPLE READ TOTALS
# ==============================================================================

sample_read_totals <- dat %>%

  group_by(
    Sample,
    Farm
  ) %>%

  summarise(
    Total_reads = sum(
      Count,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


write.csv(
  sample_read_totals,
  file.path(
    TABLE_DIR,
    "sample_read_totals.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 14. RELATIVE ABUNDANCE FUNCTION
# ==============================================================================

make_rank_abundance <- function(
    data,
    rank_name
) {

  rank_values <- data[[rank_name]]


  output <- data %>%

    mutate(
      Taxon = rank_values,

      Taxon = ifelse(
        is.na(Taxon) |
          Taxon == "",
        "Unclassified",
        Taxon
      )
    ) %>%

    group_by(
      Sample,
      Farm,
      Taxon
    ) %>%

    summarise(
      Count = sum(
        Count,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%

    group_by(
      Sample
    ) %>%

    mutate(
      Sample_total = sum(
        Count,
        na.rm = TRUE
      ),

      Relative_abundance = ifelse(
        Sample_total > 0,
        Count / Sample_total * 100,
        0
      )
    ) %>%

    ungroup() %>%

    select(
      -Sample_total
    )


  output
}


# ==============================================================================
# 15. COLLAPSE TOP TAXA FUNCTION
# ==============================================================================

collapse_top_taxa <- function(
    abundance_data,
    n_top
) {

  taxa_summary <- abundance_data %>%

    group_by(Taxon) %>%

    summarise(
      Mean_abundance = mean(
        Relative_abundance,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%

    arrange(
      desc(Mean_abundance)
    )


  n_keep <- min(
    n_top,
    nrow(taxa_summary)
  )


  top_taxa <- taxa_summary %>%

    slice_head(
      n = n_keep
    ) %>%

    pull(Taxon)


  collapsed <- abundance_data %>%

    mutate(
      Taxon_plot = ifelse(
        Taxon %in% top_taxa,
        Taxon,
        "Other"
      )
    ) %>%

    group_by(
      Sample,
      Farm,
      Taxon_plot
    ) %>%

    summarise(
      Relative_abundance = sum(
        Relative_abundance,
        na.rm = TRUE
      ),
      .groups = "drop"
    )


  categories <- unique(
    c(
      top_taxa,
      if ("Other" %in% collapsed$Taxon_plot) {
        "Other"
      }
    )
  )


  complete_grid <- expand_grid(
    Sample = sample_order,
    Taxon_plot = categories
  ) %>%

    left_join(
      sample_metadata %>%
        mutate(
          Sample = as.character(Sample)
        ),
      by = "Sample"
    )


  collapsed <- complete_grid %>%

    left_join(
      collapsed,
      by = c(
        "Sample",
        "Farm",
        "Taxon_plot"
      )
    ) %>%

    mutate(
      Relative_abundance = replace_na(
        Relative_abundance,
        0
      )
    ) %>%

    group_by(Sample) %>%

    mutate(
      Plot_total = sum(
        Relative_abundance
      ),

      Relative_abundance = ifelse(
        Plot_total > 0,
        Relative_abundance / Plot_total * 100,
        0
      )
    ) %>%

    ungroup() %>%

    select(
      -Plot_total
    ) %>%

    mutate(
      Sample = factor(
        Sample,
        levels = sample_order
      ),

      Farm = factor(
        Farm,
        levels = FARM_ORDER
      ),

      Taxon_plot = factor(
        Taxon_plot,
        levels = categories
      )
    )


  list(
    data = collapsed,
    categories = categories,
    top_taxa = top_taxa
  )
}


# ==============================================================================
# 16. FIGURE A — PHYLUM RELATIVE ABUNDANCE
# ==============================================================================

phylum <- make_rank_abundance(
  dat,
  "Phylum"
)


phylum_result <- collapse_top_taxa(
  phylum,
  TOP_PHYLA
)


phylum_plot_data <- phylum_result$data


phylum_colours <- make_taxon_palette(
  phylum_result$categories
)


p_phylum <- ggplot(
  phylum_plot_data,
  aes(
    x = Sample,
    y = Relative_abundance,
    fill = Taxon_plot
  )
) +

  geom_col(
    width = 0.86
  ) +

  facet_grid(
    . ~ Farm,
    scales = "free_x",
    space = "free_x"
  ) +

  scale_fill_manual(
    values = phylum_colours,
    drop = FALSE
  ) +

  scale_y_continuous(
    breaks = seq(
      0,
      100,
      20
    ),
    expand = c(
      0,
      0
    )
  ) +

  coord_cartesian(
    ylim = c(
      0,
      100
    )
  ) +

  labs(
    title = FIGURE_TITLES["A"],

    subtitle =
      "SILVA 144 taxonomy; samples grouped by farm",

    x = NULL,

    y =
      "Relative abundance (%)",

    fill =
      "Phylum"
  ) +

  theme_microbiome() +

  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  )


print(p_phylum)


save_figure(
  p_phylum,
  "Fig_A_phylum_relative_abundance",
  14,
  7
)


# ==============================================================================
# 17. FIGURE B — GENUS RELATIVE ABUNDANCE
# ==============================================================================

genus <- make_rank_abundance(
  dat,
  "Genus"
)


genus_result <- collapse_top_taxa(
  genus,
  TOP_GENERA
)


genus_plot_data <- genus_result$data


genus_colours <- make_taxon_palette(
  genus_result$categories
)


p_genus <- ggplot(
  genus_plot_data,
  aes(
    x = Sample,
    y = Relative_abundance,
    fill = Taxon_plot
  )
) +

  geom_col(
    width = 0.86
  ) +

  facet_grid(
    . ~ Farm,
    scales = "free_x",
    space = "free_x"
  ) +

  scale_fill_manual(
    values = genus_colours,
    drop = FALSE
  ) +

  scale_y_continuous(
    breaks = seq(
      0,
      100,
      20
    ),
    expand = c(
      0,
      0
    )
  ) +

  coord_cartesian(
    ylim = c(
      0,
      100
    )
  ) +

  labs(
    title = FIGURE_TITLES["B"],

    subtitle =
      "SILVA 144 taxonomy; samples grouped by farm",

    x = NULL,

    y =
      "Relative abundance (%)",

    fill =
      "Genus"
  ) +

  theme_microbiome() +

  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),

    legend.text =
      element_text(
        face = "italic"
      )
  )


print(p_genus)


save_figure(
  p_genus,
  "Fig_B_genus_relative_abundance",
  15,
  7
)


# ==============================================================================
# 18. FIGURE C — TOP 20 GENERA HEATMAP
# ==============================================================================

identified_genus <- genus %>%

  filter(
    Taxon !=
      "Unclassified"
  )


genus_summary <- identified_genus %>%

  group_by(Taxon) %>%

  summarise(
    Mean_abundance = mean(
      Relative_abundance,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%

  arrange(
    desc(
      Mean_abundance
    )
  )


n_heatmap <- min(
  TOP_HEATMAP_GENERA,
  nrow(genus_summary)
)


top_heatmap_genera <- genus_summary %>%

  slice_head(
    n = n_heatmap
  ) %>%

  pull(Taxon)


heatmap_grid <- expand_grid(
  Sample = sample_order,
  Genus = top_heatmap_genera
) %>%

  left_join(
    sample_metadata %>%
      mutate(
        Sample = as.character(Sample)
      ),
    by = "Sample"
  )


heatmap_values <- identified_genus %>%

  filter(
    Taxon %in%
      top_heatmap_genera
  ) %>%

  transmute(
    Sample,
    Genus = Taxon,
    Relative_abundance
  )


heatmap_data <- heatmap_grid %>%

  left_join(
    heatmap_values,
    by = c(
      "Sample",
      "Genus"
    )
  ) %>%

  mutate(
    Relative_abundance = replace_na(
      Relative_abundance,
      0
    ),

    Log_abundance = log10(
      Relative_abundance +
        0.01
    ),

    Sample = factor(
      Sample,
      levels = sample_order
    ),

    Farm = factor(
      Farm,
      levels = FARM_ORDER
    ),

    Genus = factor(
      Genus,
      levels = rev(
        top_heatmap_genera
      )
    )
  )


p_heatmap <- ggplot(
  heatmap_data,
  aes(
    x = Sample,
    y = Genus,
    fill = Log_abundance
  )
) +

  geom_tile(
    colour = "white",
    linewidth = 0.35
  ) +

  facet_grid(
    . ~ Farm,
    scales = "free_x",
    space = "free_x"
  ) +

  scale_fill_gradientn(
    colours = viridis(
      100,
      option = "C"
    ),

    name =
      expression(
        log[10] *
          "(relative abundance + 0.01)"
      )
  ) +

  labs(
    title = FIGURE_TITLES["C"],

    subtitle =
      "Top genera ranked by mean relative abundance",

    x = NULL,
    y = NULL
  ) +

  theme_heatmap_clean() +

  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),

    axis.text.y = element_text(
      face = "italic"
    )
  )


print(p_heatmap)


save_figure(
  p_heatmap,
  "Fig_C_top20_genera_heatmap",
  14,
  8
)


# ==============================================================================
# 19. SAVE TAXONOMIC TABLES
# ==============================================================================

write.csv(
  phylum,
  file.path(
    TABLE_DIR,
    "phylum_relative_abundance.csv"
  ),
  row.names = FALSE
)


write.csv(
  genus,
  file.path(
    TABLE_DIR,
    "genus_relative_abundance.csv"
  ),
  row.names = FALSE
)


write.csv(
  heatmap_data,
  file.path(
    TABLE_DIR,
    "top20_genera_heatmap_data.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 20. PREPARE TRUE GENUS-LEVEL NETWORK DATA
# ==============================================================================

network_genus <- dat %>%

  filter(
    !is.na(Genus)
  ) %>%

  group_by(
    Farm,
    Sample,
    Genus
  ) %>%

  summarise(
    Count = sum(
      Count,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# ==============================================================================
# 21. GLOBAL PREVALENCE FILTER
# ==============================================================================

n_samples_total <- length(
  sample_order
)


prevalence <- network_genus %>%

  group_by(Genus) %>%

  summarise(
    Present = n_distinct(
      Sample[
        Count > 0
      ]
    ),

    Prevalence =
      Present /
      n_samples_total,

    .groups = "drop"
  )


keep_genera <- prevalence %>%

  filter(
    Prevalence >=
      PREVALENCE_THRESHOLD
  ) %>%

  pull(Genus)


cat(
  "\nGenera retained after prevalence filtering:",
  length(keep_genera),
  "\n"
)


if (length(keep_genera) < 2) {

  stop(
    "Fewer than two genera remained after prevalence filtering."
  )
}


write.csv(
  prevalence,
  file.path(
    TABLE_DIR,
    "network_genus_prevalence.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 22. BUILD SAMPLE × GENUS MATRIX
# ==============================================================================

network_counts <- network_genus %>%

  filter(
    Genus %in%
      keep_genera
  ) %>%

  group_by(
    Sample,
    Genus
  ) %>%

  summarise(
    Count = sum(Count),
    .groups = "drop"
  )


network_grid <- expand_grid(
  Sample = sample_order,
  Genus = keep_genera
)


network_long <- network_grid %>%

  left_join(
    network_counts,
    by = c(
      "Sample",
      "Genus"
    )
  ) %>%

  mutate(
    Count = replace_na(
      Count,
      0
    )
  )


network_matrix_df <- network_long %>%

  pivot_wider(
    names_from = Genus,
    values_from = Count,
    values_fill = 0
  )


network_sample_names <-
  network_matrix_df$Sample


network_matrix <- network_matrix_df %>%

  select(
    -Sample
  ) %>%

  as.matrix()


rownames(network_matrix) <-
  network_sample_names


# ==============================================================================
# 23. CLR TRANSFORMATION
# ==============================================================================

network_pc <-
  network_matrix +
  PSEUDOCOUNT


clr_matrix <- t(
  apply(
    network_pc,
    1,
    function(x) {
      log(x) -
        mean(
          log(x)
        )
    }
  )
)


# ==============================================================================
# 24. PAIRWISE SPEARMAN ASSOCIATIONS
# ==============================================================================

taxa_names <- colnames(
  clr_matrix
)


taxon_pairs <- combn(
  taxa_names,
  2,
  simplify = FALSE
)


safe_correlation <- function(pair) {

  x <- clr_matrix[, pair[1]]

  y <- clr_matrix[, pair[2]]


  if (
    sd(x) == 0 ||
      sd(y) == 0
  ) {

    return(
      tibble(
        Taxon1 = pair[1],
        Taxon2 = pair[2],
        rho = NA_real_,
        p = NA_real_
      )
    )
  }


  result <- suppressWarnings(
    cor.test(
      x,
      y,
      method = CORRELATION_METHOD,
      exact = FALSE
    )
  )


  tibble(
    Taxon1 = pair[1],
    Taxon2 = pair[2],
    rho = unname(
      result$estimate
    ),
    p = result$p.value
  )
}


edge_tests <- map_dfr(
  taxon_pairs,
  safe_correlation
)


edge_tests <- edge_tests %>%

  mutate(
    p_adj = p.adjust(
      p,
      method = "BH"
    )
  )


network_edges <- edge_tests %>%

  filter(
    !is.na(rho),
    !is.na(p_adj),

    abs(rho) >=
      CORRELATION_THRESHOLD,

    p_adj <
      FDR_THRESHOLD
  ) %>%

  mutate(
    Association = ifelse(
      rho > 0,
      "Positive",
      "Negative"
    ),

    Weight = abs(rho)
  )


cat(
  "Significant global edges:",
  nrow(network_edges),
  "\n"
)


write.csv(
  edge_tests,
  file.path(
    TABLE_DIR,
    "all_pairwise_associations.csv"
  ),
  row.names = FALSE
)


write.csv(
  network_edges,
  file.path(
    TABLE_DIR,
    "significant_global_network_edges.csv"
  ),
  row.names = FALSE
)


if (nrow(network_edges) == 0) {

  stop(
    "No network associations passed the selected thresholds."
  )
}


# ==============================================================================
# 25. BUILD GLOBAL IGRAPH NETWORK
# ==============================================================================

g <- graph_from_data_frame(
  network_edges,
  directed = FALSE
)


# ==============================================================================
# 26. GLOBAL COMMUNITY DETECTION
# ==============================================================================

communities <- cluster_louvain(
  g,
  weights = E(g)$Weight
)


V(g)$Module <- membership(
  communities
)


network_modularity <- modularity(
  communities
)


# ==============================================================================
# 27. GLOBAL NODE STATISTICS
# ==============================================================================

node_stats <- tibble(
  Genus = V(g)$name,

  Degree = degree(g),

  Betweenness = betweenness(
    g,
    directed = FALSE,
    normalized = TRUE
  ),

  Eigenvector = eigen_centrality(
    g,
    directed = FALSE,
    weights = E(g)$Weight
  )$vector,

  Module = V(g)$Module
) %>%

  arrange(
    desc(Degree),
    desc(Betweenness)
  )


write.csv(
  node_stats,
  file.path(
    TABLE_DIR,
    "global_network_node_centrality.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 28. GLOBAL NETWORK SUMMARY
# ==============================================================================

network_summary <- tibble(
  Metric = c(
    "Samples",
    "Genera after prevalence filtering",
    "Network nodes",
    "Network edges",
    "Network density",
    "Mean degree",
    "Modules",
    "Modularity"
  ),

  Value = c(
    n_samples_total,
    length(keep_genera),
    vcount(g),
    ecount(g),
    edge_density(
      g,
      loops = FALSE
    ),
    mean(
      degree(g)
    ),
    length(communities),
    network_modularity
  )
)


write.csv(
  network_summary,
  file.path(
    TABLE_DIR,
    "global_network_summary.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 29. GLOBAL NETWORK APPEARANCE
# ==============================================================================

module_ids <- sort(
  unique(
    V(g)$Module
  )
)


MODULE_COLOURS <- setNames(
  hcl.colors(
    length(module_ids),
    palette = "Set 2"
  ),
  as.character(module_ids)
)


V(g)$color <- unname(
  MODULE_COLOURS[
    as.character(
      V(g)$Module
    )
  ]
)


V(g)$frame.color <- "white"


V(g)$size <- 5 +
  safe_rescale(
    degree(g),
    to = c(
      2,
      12
    ),
    constant_value = 7
  )


E(g)$color <- ifelse(
  E(g)$rho > 0,
  EDGE_POSITIVE,
  EDGE_NEGATIVE
)


E(g)$width <- safe_rescale(
  abs(
    E(g)$rho
  ),
  to = c(
    0.6,
    3
  ),
  constant_value = 1.5
)


E(g)$lty <- ifelse(
  E(g)$rho > 0,
  1,
  2
)


# ==============================================================================
# 30. GLOBAL NETWORK LAYOUT
# ==============================================================================

set.seed(123)


global_layout <- layout_with_fr(
  g,
  weights = E(g)$Weight
)


top_n_nodes <- min(
  20L,
  nrow(node_stats)
)


global_label_taxa <- node_stats %>%

  slice_max(
    order_by = Degree,
    n = top_n_nodes,
    with_ties = FALSE
  ) %>%

  pull(Genus)


V(g)$label_hubs <- ifelse(
  V(g)$name %in%
    global_label_taxa,
  V(g)$name,
  NA_character_
)


# ==============================================================================
# 31. NETWORK DRAWING FUNCTION
# ==============================================================================

draw_network <- function(
    graph_object,
    layout_object,
    labels,
    title_text,
    subtitle_text,
    farm_border = NULL
) {

  if (is.null(farm_border)) {

    vertex_border <- "white"

  } else {

    vertex_border <- farm_border
  }


  plot(
    graph_object,

    layout = layout_object,

    vertex.label = labels,

    vertex.label.cex = 0.68,

    vertex.label.font = 3,

    vertex.label.color = "black",

    vertex.size =
      V(graph_object)$size,

    vertex.color =
      V(graph_object)$color,

    vertex.frame.color =
      vertex_border,

    vertex.frame.width = 1.5,

    edge.color =
      E(graph_object)$color,

    edge.width =
      E(graph_object)$width,

    edge.lty =
      E(graph_object)$lty,

    edge.curved = 0.08,

    main = ""
  )


  title(
    main = title_text,
    sub = subtitle_text,
    cex.main = 1.25,
    font.main = 2,
    cex.sub = 0.85
  )


  legend(
    "topleft",

    legend = c(
      "Positive association",
      "Negative association"
    ),

    col = c(
      EDGE_POSITIVE,
      EDGE_NEGATIVE
    ),

    lty = c(
      1,
      2
    ),

    lwd = 2,

    bty = "n",

    cex = 0.85
  )


  modules_present <- sort(
    unique(
      V(graph_object)$Module
    )
  )


  legend(
    "topright",

    legend = paste(
      "Module",
      modules_present
    ),

    pch = 21,

    pt.bg = MODULE_COLOURS[
      as.character(
        modules_present
      )
    ],

    pt.cex = 1.4,

    bty = "n",

    cex = 0.8
  )
}


# ==============================================================================
# 32. FIGURE D — GLOBAL NETWORK
# ==============================================================================

png(
  filename = file.path(
    FIGURE_DIR,
    "Fig_D_global_association_network.png"
  ),
  width = 4500,
  height = 3800,
  res = 400
)


draw_network(
  graph_object = g,
  layout_object = global_layout,
  labels = V(g)$label_hubs,
  title_text = FIGURE_TITLES["D"],
  subtitle_text = paste0(
    "|rho| >= ",
    CORRELATION_THRESHOLD,
    "; FDR < ",
    FDR_THRESHOLD,
    "; node colour = network module"
  )
)


dev.off()


pdf(
  file = file.path(
    FIGURE_DIR,
    "Fig_D_global_association_network.pdf"
  ),
  width = 14,
  height = 12
)


draw_network(
  graph_object = g,
  layout_object = global_layout,
  labels = V(g)$label_hubs,
  title_text = FIGURE_TITLES["D"],
  subtitle_text = paste0(
    "|rho| >= ",
    CORRELATION_THRESHOLD,
    "; FDR < ",
    FDR_THRESHOLD,
    "; node colour = network module"
  )
)


dev.off()


# ==============================================================================
# 33. FIGURE E — MOST CONNECTED GENERA
# ==============================================================================

top_nodes <- node_stats %>%

  slice_max(
    order_by = Degree,
    n = top_n_nodes,
    with_ties = FALSE
  ) %>%

  arrange(Degree) %>%

  mutate(
    Genus = factor(
      Genus,
      levels = Genus
    ),

    Module = factor(Module)
  )


p_degree <- ggplot(
  top_nodes,
  aes(
    x = Degree,
    y = Genus,
    fill = Module
  )
) +

  geom_col(
    width = 0.72
  ) +

  scale_fill_manual(
    values = MODULE_COLOURS
  ) +

  labs(
    title = FIGURE_TITLES["E"],

    subtitle =
      "Node degree represents the number of retained associations",

    x = "Degree",

    y = NULL,

    fill = "Module"
  ) +

  theme_microbiome() +

  theme(
    axis.text.y = element_text(
      face = "italic"
    )
  )


print(p_degree)


save_figure(
  p_degree,
  "Fig_E_top20_global_network_genera",
  10,
  8
)


# ==============================================================================
# 34. EXPORT GLOBAL NETWORK
# ==============================================================================

write_graph(
  g,
  file.path(
    NETWORK_DIR,
    "global_association_network.graphml"
  ),
  format = "graphml"
)


# ==============================================================================
# 35. FARM SUBNETWORK FUNCTION
#
# Correlations are NOT recalculated within farms.
# Each farm network is an induced subgraph of the GLOBAL network.
# ==============================================================================

make_farm_subnetwork <- function(
    farm_name
) {

  figure_code <-
    FARM_FIGURE_CODES[
      farm_name
    ]


  taxa_here <- network_genus %>%

    filter(
      Farm == farm_name,
      Count > 0,
      Genus %in%
        V(g)$name
    ) %>%

    group_by(Genus) %>%

    summarise(
      Samples_present =
        n_distinct(Sample),
      .groups = "drop"
    ) %>%

    filter(
      Samples_present >=
        FARM_MIN_SAMPLES
    ) %>%

    pull(Genus)


  taxa_here <- intersect(
    taxa_here,
    V(g)$name
  )


  if (length(taxa_here) < 2) {

    warning(
      paste(
        farm_name,
        "has fewer than two eligible genera."
      )
    )

    return(NULL)
  }


  g_farm <- induced_subgraph(
    g,
    vids = taxa_here
  )


  if (ecount(g_farm) == 0) {

    warning(
      paste(
        farm_name,
        "has no globally supported network edges after filtering."
      )
    )

    return(NULL)
  }


  farm_degree <- degree(
    g_farm
  )


  V(g_farm)$size <- 5 +
    safe_rescale(
      farm_degree,
      to = c(
        2,
        12
      ),
      constant_value = 7
    )


  V(g_farm)$color <- unname(
    MODULE_COLOURS[
      as.character(
        V(g_farm)$Module
      )
    ]
  )


  V(g_farm)$frame.color <-
    FARM_COLOURS[
      farm_name
    ]


  E(g_farm)$color <- ifelse(
    E(g_farm)$rho > 0,
    EDGE_POSITIVE,
    EDGE_NEGATIVE
  )


  E(g_farm)$width <- safe_rescale(
    abs(
      E(g_farm)$rho
    ),
    to = c(
      0.6,
      3
    ),
    constant_value = 1.5
  )


  E(g_farm)$lty <- ifelse(
    E(g_farm)$rho > 0,
    1,
    2
  )


  set.seed(123)


  farm_layout <- layout_with_fr(
    g_farm,
    weights = E(g_farm)$Weight
  )


  n_farm_labels <- min(
    15L,
    vcount(g_farm)
  )


  label_order <- names(
    sort(
      farm_degree,
      decreasing = TRUE
    )
  )


  label_taxa <- label_order[
    seq_len(
      n_farm_labels
    )
  ]


  farm_labels <- ifelse(
    V(g_farm)$name %in%
      label_taxa,
    V(g_farm)$name,
    NA_character_
  )


  figure_title <-
    FIGURE_TITLES[
      figure_code
    ]


  figure_subtitle <- paste0(
    "Global-network edges among genera detected in >= ",
    FARM_MIN_SAMPLES,
    " ",
    farm_name,
    " samples"
  )


  file_stub <- paste0(
    "Fig_",
    figure_code,
    "_",
    farm_name,
    "_association_subnetwork"
  )


  # ---------------------------------------------------------------------------
  # PNG
  # ---------------------------------------------------------------------------

  png(
    filename = file.path(
      FIGURE_DIR,
      paste0(
        file_stub,
        ".png"
      )
    ),
    width = 3800,
    height = 3400,
    res = 400
  )


  draw_network(
    graph_object = g_farm,
    layout_object = farm_layout,
    labels = farm_labels,
    title_text = figure_title,
    subtitle_text = figure_subtitle,
    farm_border = FARM_COLOURS[farm_name]
  )


  dev.off()


  # ---------------------------------------------------------------------------
  # PDF
  # ---------------------------------------------------------------------------

  pdf(
    file = file.path(
      FIGURE_DIR,
      paste0(
        file_stub,
        ".pdf"
      )
    ),
    width = 12,
    height = 11
  )


  draw_network(
    graph_object = g_farm,
    layout_object = farm_layout,
    labels = farm_labels,
    title_text = figure_title,
    subtitle_text = figure_subtitle,
    farm_border = FARM_COLOURS[farm_name]
  )


  dev.off()


  # ---------------------------------------------------------------------------
  # Node statistics
  # ---------------------------------------------------------------------------

  farm_nodes <- tibble(
    Farm = farm_name,

    Genus = V(g_farm)$name,

    Degree = degree(
      g_farm
    ),

    Betweenness = betweenness(
      g_farm,
      directed = FALSE,
      normalized = TRUE
    ),

    Module = V(g_farm)$Module
  ) %>%

    arrange(
      desc(Degree),
      desc(Betweenness)
    )


  write.csv(
    farm_nodes,
    file.path(
      TABLE_DIR,
      paste0(
        "network_nodes_",
        farm_name,
        ".csv"
      )
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # GraphML
  # ---------------------------------------------------------------------------

  write_graph(
    g_farm,
    file.path(
      FARM_NETWORK_DIR,
      paste0(
        farm_name,
        "_association_subnetwork.graphml"
      )
    ),
    format = "graphml"
  )


  # ---------------------------------------------------------------------------
  # Return summary
  # ---------------------------------------------------------------------------

  tibble(
    Farm = farm_name,

    Minimum_samples_required =
      FARM_MIN_SAMPLES,

    Nodes =
      vcount(g_farm),

    Edges =
      ecount(g_farm),

    Density =
      edge_density(
        g_farm,
        loops = FALSE
      ),

    Mean_degree =
      mean(
        degree(g_farm)
      )
  )
}


# ==============================================================================
# 36. FIGURES F–I — FARM SUBNETWORKS
# ==============================================================================

farm_network_results <- vector(
  "list",
  length(
    FARM_ORDER
  )
)


names(
  farm_network_results
) <- FARM_ORDER


for (farm_name in FARM_ORDER) {

  farm_network_results[
    [farm_name]
  ] <- list(
    make_farm_subnetwork(
      farm_name
    )
  )
}


farm_network_summaries <- bind_rows(
  farm_network_results
)


write.csv(
  farm_network_summaries,
  file.path(
    TABLE_DIR,
    "farm_subnetwork_summary.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 37. SAMPLE REPRESENTATION CHECK
# ==============================================================================

sample_check <- tibble(
  Sample = as.character(
    sample_metadata$Sample
  ),

  Farm = as.character(
    sample_metadata$Farm
  )
) %>%

  mutate(
    In_feature_table =
      Sample %in%
      sample_order,

    In_phylum_plot =
      Sample %in%
      as.character(
        unique(
          phylum_plot_data$Sample
        )
      ),

    In_genus_plot =
      Sample %in%
      as.character(
        unique(
          genus_plot_data$Sample
        )
      ),

    In_heatmap =
      Sample %in%
      as.character(
        unique(
          heatmap_data$Sample
        )
      ),

    In_global_network_input =
      Sample %in%
      unique(
        network_long$Sample
      )
  )


write.csv(
  sample_check,
  file.path(
    TABLE_DIR,
    "sample_representation_check.csv"
  ),
  row.names = FALSE
)


cat("\n========================================\n")
cat("SAMPLE REPRESENTATION CHECK\n")
cat("========================================\n")

print(
  sample_check,
  n = Inf
)


# ==============================================================================
# 38. FIGURE MANIFEST
# ==============================================================================

figure_manifest <- tibble(
  Figure = LETTERS[1:9],

  Title = unname(
    FIGURE_TITLES[
      LETTERS[1:9]
    ]
  ),

  File = c(
    "Fig_A_phylum_relative_abundance",
    "Fig_B_genus_relative_abundance",
    "Fig_C_top20_genera_heatmap",
    "Fig_D_global_association_network",
    "Fig_E_top20_global_network_genera",
    "Fig_F_Clanwilliam_association_subnetwork",
    "Fig_G_Dendron_association_subnetwork",
    "Fig_H_Mamusha_association_subnetwork",
    "Fig_I_Wesselesbron_association_subnetwork"
  )
)


write.csv(
  figure_manifest,
  file.path(
    OUTPUT_DIR,
    "figure_manifest.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 39. ANALYSIS PARAMETERS
# ==============================================================================

analysis_parameters <- tibble(
  Parameter = c(
    "SILVA release",
    "Top phyla",
    "Top genera",
    "Heatmap genera",
    "Global prevalence threshold",
    "Correlation method",
    "Absolute correlation threshold",
    "FDR threshold",
    "CLR pseudocount",
    "Minimum farm samples for subnetwork inclusion"
  ),

  Value = c(
    "144",
    TOP_PHYLA,
    TOP_GENERA,
    TOP_HEATMAP_GENERA,
    PREVALENCE_THRESHOLD,
    CORRELATION_METHOD,
    CORRELATION_THRESHOLD,
    FDR_THRESHOLD,
    PSEUDOCOUNT,
    FARM_MIN_SAMPLES
  )
)


write.csv(
  analysis_parameters,
  file.path(
    OUTPUT_DIR,
    "analysis_parameters.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 40. SAVE R SESSION INFORMATION
# ==============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    OUTPUT_DIR,
    "R_sessionInfo.txt"
  )
)


# ==============================================================================
# 41. FINISHED
# ==============================================================================

cat("\n========================================\n")
cat("ANALYSIS COMPLETE\n")
cat("========================================\n\n")


cat(
  "Figures:\n",
  FIGURE_DIR,
  "\n\n"
)


cat(
  "Tables:\n",
  TABLE_DIR,
  "\n\n"
)


cat(
  "Network files:\n",
  NETWORK_DIR,
  "\n\n"
)


cat(
  "Standardised figure series:\n",
  "(A) Phylum-level relative abundance\n",
  "(B) Genus-level relative abundance\n",
  "(C) Top 20 genera heatmap\n",
  "(D) Global genus-level association network\n",
  "(E) Most connected genera in the global network\n",
  "(F) Clanwilliam association subnetwork\n",
  "(G) Dendron association subnetwork\n",
  "(H) Mamusha association subnetwork\n",
  "(I) Wesselesbron association subnetwork\n"
)
