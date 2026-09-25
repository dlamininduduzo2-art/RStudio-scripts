# ==============================================================================
# SILVA 144 TAXONOMIC COMPOSITION + CO-OCCURRENCE NETWORK ANALYSIS
# Potato rhizosphere bacterial microbiome
#
# INPUTS:
#   R_taxonomy_plots/feature-table.tsv
#   R_taxonomy_plots/taxonomy144/taxonomy.tsv
#
# OUTPUTS:
#   R_taxonomy_plots/final_outputs/
#
# Produces:
#   1. Phylum relative-abundance plot
#   2. Genus relative-abundance plot
#   3. Top-20 genus heatmap
#   4. Genus-level co-occurrence network
#   5. Hub-taxa network
#   6. Hub degree plot
#   7. Network statistics
#   8. Node/edge CSV files
#   9. GraphML file for Cytoscape / Gephi
#
# IMPORTANT:
# Co-occurrence associations do not prove direct ecological interactions.
# ==============================================================================


# ==============================================================================
# 1. USER SETTINGS
# ==============================================================================

FEATURE_FILE <- "R_taxonomy_plots/feature-table.tsv"

TAXONOMY_FILE <- "R_taxonomy_plots/taxonomy144/taxonomy.tsv"

OUTPUT_DIR <- "R_taxonomy_plots/final_outputs"


# Number of taxa displayed
TOP_PHYLA <- 10
TOP_GENERA <- 15
TOP_HEATMAP_GENERA <- 20


# Network parameters
PREVALENCE_THRESHOLD <- 0.20
CORRELATION_THRESHOLD <- 0.60
FDR_THRESHOLD <- 0.05
PSEUDOCOUNT <- 0.5


# Remove chloroplast, mitochondria and non-bacterial sequences
REMOVE_NON_TARGET <- TRUE


# Preserve feature-table sample order by default
CUSTOM_SAMPLE_ORDER <- NULL


# ==============================================================================
# 2. PACKAGES
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

  install.packages(
    missing_packages
  )
}

library(tidyverse)
library(igraph)
library(scales)
library(viridisLite)


# ==============================================================================
# 3. CREATE OUTPUT DIRECTORY
# ==============================================================================

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 4. IMPORT FEATURE TABLE
# ==============================================================================

asv <- read.delim(
  FEATURE_FILE,
  skip = 1,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

colnames(asv)[1] <- "FeatureID"


all_samples <- colnames(asv)[-1]


if (is.null(CUSTOM_SAMPLE_ORDER)) {

  sample_order <- all_samples

} else {

  sample_order <- CUSTOM_SAMPLE_ORDER

  missing_samples <- setdiff(
    all_samples,
    sample_order
  )

  if (length(missing_samples) > 0) {

    sample_order <- c(
      sample_order,
      missing_samples
    )
  }
}


cat("\n========================================\n")
cat("FEATURE TABLE\n")
cat("========================================\n")

cat(
  "ASVs:",
  nrow(asv),
  "\n"
)

cat(
  "Samples:",
  length(all_samples),
  "\n\n"
)

print(sample_order)


# ==============================================================================
# 5. IMPORT SILVA 144 TAXONOMY
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
    paste(
      "No column named 'Taxon' found.",
      "Available columns:",
      paste(
        colnames(tax),
        collapse = ", "
      )
    )
  )
}


# ==============================================================================
# 6. SILVA TAXONOMY PARSER
#
# Handles both:
#
# d__Bacteria; p__Pseudomonadota...
#
# AND
#
# D_0__Bacteria; D_1__Proteobacteria...
#
# ==============================================================================

extract_rank <- function(
    x,
    modern_prefix,
    old_rank
) {

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


  modern <- str_match(
    x,
    modern_pattern
  )[, 2]


  old <- str_match(
    x,
    old_pattern
  )[, 2]


  result <- ifelse(
    !is.na(modern),
    modern,
    old
  )


  result <- str_trim(
    result
  )


  result[
    result == "" |
      result == "Unassigned"
  ] <- NA


  return(result)
}


tax_clean <- tax %>%

  transmute(

    FeatureID,

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


# ==============================================================================
# 7. TAXONOMY RESOLUTION SUMMARY
# ==============================================================================

taxonomy_resolution <- tax_clean %>%

  summarise(

    Total_ASVs = n(),

    Domain_resolved =
      sum(!is.na(Domain)),

    Phylum_resolved =
      sum(!is.na(Phylum)),

    Family_resolved =
      sum(!is.na(Family)),

    Genus_resolved =
      sum(!is.na(Genus)),

    Species_resolved =
      sum(!is.na(Species))
  )


cat("\n========================================\n")
cat("TAXONOMIC RESOLUTION\n")
cat("========================================\n")

print(
  taxonomy_resolution
)


write.csv(
  taxonomy_resolution,
  file.path(
    OUTPUT_DIR,
    "taxonomy_resolution_summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  tax_clean,
  file.path(
    OUTPUT_DIR,
    "SILVA144_cleaned_taxonomy.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 8. CONVERT FEATURE TABLE TO LONG FORMAT
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


# ==============================================================================
# 9. JOIN TAXONOMY TO COUNTS
# ==============================================================================

dat <- asv_long %>%

  left_join(
    tax_clean,
    by = "FeatureID"
  )


# ==============================================================================
# 10. REMOVE NON-TARGET FEATURES
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
# 11. SAMPLE READ TOTALS
# ==============================================================================

sample_counts <- dat %>%

  group_by(Sample) %>%

  summarise(

    Total_reads =
      sum(
        Count,
        na.rm = TRUE
      ),

    .groups = "drop"
  )


sample_counts <- tibble(
  Sample = sample_order
) %>%

  left_join(
    sample_counts,
    by = "Sample"
  ) %>%

  mutate(
    Total_reads =
      replace_na(
        Total_reads,
        0
      )
  )


cat("\n========================================\n")
cat("SAMPLE READ TOTALS\n")
cat("========================================\n")

print(
  sample_counts,
  n = Inf
)


write.csv(
  sample_counts,
  file.path(
    OUTPUT_DIR,
    "sample_read_totals.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 12. FUNCTION: CALCULATE RELATIVE ABUNDANCE
# ==============================================================================

make_rank_abundance <- function(
    data,
    rank_name
) {

  result <- data %>%

    mutate(

      Taxon = coalesce(
        .data[[rank_name]],
        "Unclassified"
      )
    ) %>%

    group_by(
      Sample,
      Taxon
    ) %>%

    summarise(

      Count =
        sum(
          Count,
          na.rm = TRUE
        ),

      .groups = "drop"
    ) %>%

    group_by(Sample) %>%

    mutate(

      Sample_total =
        sum(
          Count,
          na.rm = TRUE
        ),

      Relative_abundance =
        ifelse(
          Sample_total > 0,
          Count /
            Sample_total *
            100,
          0
        )
    ) %>%

    ungroup() %>%

    select(
      -Sample_total
    )


  return(result)
}


# ==============================================================================
# 13. FUNCTION: COLLAPSE TO TOP TAXA + OTHER
# ==============================================================================

collapse_top_taxa <- function(
    abundance_data,
    n_top
) {

  rank_summary <- abundance_data %>%

    group_by(Taxon) %>%

    summarise(

      Mean_abundance =
        mean(
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
    nrow(rank_summary)
  )


  top_taxa <- rank_summary %>%

    slice_head(
      n = n_keep
    ) %>%

    pull(Taxon)


  collapsed <- abundance_data %>%

    mutate(

      Taxon_plot =
        if_else(
          Taxon %in% top_taxa,
          Taxon,
          "Other"
        )
    ) %>%

    group_by(
      Sample,
      Taxon_plot
    ) %>%

    summarise(

      Relative_abundance =
        sum(
          Relative_abundance,
          na.rm = TRUE
        ),

      .groups = "drop"
    )


  categories <- c(
    top_taxa,
    if ("Other" %in% collapsed$Taxon_plot)
      "Other"
  )


  categories <- unique(
    categories
  )


  full_grid <- expand_grid(

    Sample =
      sample_order,

    Taxon_plot =
      categories
  )


  collapsed <- full_grid %>%

    left_join(

      collapsed,

      by = c(
        "Sample",
        "Taxon_plot"
      )
    ) %>%

    mutate(

      Relative_abundance =
        replace_na(
          Relative_abundance,
          0
        )
    ) %>%

    group_by(Sample) %>%

    mutate(

      Plot_total =
        sum(
          Relative_abundance,
          na.rm = TRUE
        ),

      Relative_abundance =
        ifelse(
          Plot_total > 0,
          Relative_abundance /
            Plot_total *
            100,
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

      Taxon_plot = factor(
        Taxon_plot,
        levels = categories
      )
    )


  return(
    list(
      data = collapsed,
      top_taxa = top_taxa,
      categories = categories
    )
  )
}


# ==============================================================================
# 14. FUNCTION: CREATE COLOURS
# ==============================================================================

make_taxon_colours <- function(
    categories
) {

  main_taxa <- setdiff(
    categories,
    "Other"
  )


  colours <- setNames(

    hcl.colors(
      length(main_taxa),
      palette = "Dark 3"
    ),

    main_taxa
  )


  if ("Other" %in% categories) {

    colours <- c(
      colours,
      "Other" = "grey82"
    )
  }


  return(colours)
}


# ==============================================================================
# 15. PHYLUM RELATIVE ABUNDANCE
# ==============================================================================

phylum <- make_rank_abundance(
  dat,
  "Phylum"
)


phylum_result <- collapse_top_taxa(
  phylum,
  TOP_PHYLA
)


phylum_plot_data <-
  phylum_result$data


phylum_colours <- make_taxon_colours(
  phylum_result$categories
)


phylum_totals <- phylum_plot_data %>%

  group_by(Sample) %>%

  summarise(

    Total =
      sum(
        Relative_abundance,
        na.rm = TRUE
      ),

    .groups = "drop"
  )


cat("\n========================================\n")
cat("PHYLUM PLOT TOTALS\n")
cat("========================================\n")

print(
  phylum_totals,
  n = Inf
)


# ==============================================================================
# 16. PHYLUM PLOT
# ==============================================================================

p_phylum <- ggplot(

  phylum_plot_data,

  aes(

    x = Sample,

    y = Relative_abundance,

    fill = Taxon_plot
  )

) +

  geom_col(
    width = 0.85
  ) +

  scale_fill_manual(
    values = phylum_colours
  ) +

  scale_y_continuous(

    breaks =
      seq(
        0,
        100,
        20
      ),

    expand =
      c(
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

  scale_x_discrete(
    drop = FALSE
  ) +

  labs(

    x = NULL,

    y = "Relative abundance (%)",

    fill = "Phylum"
  ) +

  theme_classic(
    base_size = 12
  ) +

  theme(

    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        colour = "black"
      ),

    axis.text.y =
      element_text(
        colour = "black"
      ),

    axis.title.y =
      element_text(
        face = "bold"
      ),

    legend.title =
      element_text(
        face = "bold"
      ),

    legend.position =
      "right"
  )


print(p_phylum)


ggsave(

  file.path(
    OUTPUT_DIR,
    "SILVA144_phylum_relative_abundance.png"
  ),

  p_phylum,

  width = 13,

  height = 7,

  dpi = 600
)


ggsave(

  file.path(
    OUTPUT_DIR,
    "SILVA144_phylum_relative_abundance.pdf"
  ),

  p_phylum,

  width = 13,

  height = 7
)


# ==============================================================================
# 17. GENUS RELATIVE ABUNDANCE
# ==============================================================================

genus <- make_rank_abundance(
  dat,
  "Genus"
)


genus_result <- collapse_top_taxa(
  genus,
  TOP_GENERA
)


genus_plot_data <-
  genus_result$data


genus_colours <- make_taxon_colours(
  genus_result$categories
)


genus_totals <- genus_plot_data %>%

  group_by(Sample) %>%

  summarise(

    Total =
      sum(
        Relative_abundance,
        na.rm = TRUE
      ),

    .groups = "drop"
  )


cat("\n========================================\n")
cat("GENUS PLOT TOTALS\n")
cat("========================================\n")

print(
  genus_totals,
  n = Inf
)


# ==============================================================================
# 18. GENUS PLOT
# ==============================================================================

p_genus <- ggplot(

  genus_plot_data,

  aes(

    x = Sample,

    y = Relative_abundance,

    fill = Taxon_plot
  )

) +

  geom_col(
    width = 0.85
  ) +

  scale_fill_manual(
    values = genus_colours
  ) +

  scale_y_continuous(

    breaks =
      seq(
        0,
        100,
        20
      ),

    expand =
      c(
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

  scale_x_discrete(
    drop = FALSE
  ) +

  labs(

    x = NULL,

    y = "Relative abundance (%)",

    fill = "Genus"
  ) +

  theme_classic(
    base_size = 12
  ) +

  theme(

    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        colour = "black"
      ),

    axis.text.y =
      element_text(
        colour = "black"
      ),

    axis.title.y =
      element_text(
        face = "bold"
      ),

    legend.title =
      element_text(
        face = "bold"
      ),

    legend.text =
      element_text(
        face = "italic"
      ),

    legend.position =
      "right"
  )


print(p_genus)


ggsave(

  file.path(
    OUTPUT_DIR,
    "SILVA144_genus_relative_abundance.png"
  ),

  p_genus,

  width = 14,

  height = 7,

  dpi = 600
)


ggsave(

  file.path(
    OUTPUT_DIR,
    "SILVA144_genus_relative_abundance.pdf"
  ),

  p_genus,

  width = 14,

  height = 7
)


# ==============================================================================
# 19. TOP-20 GENUS HEATMAP
# ==============================================================================

genus_identified <- genus %>%

  filter(
    Taxon != "Unclassified"
  )


genus_summary <- genus_identified %>%

  group_by(Taxon) %>%

  summarise(

    Mean_abundance =
      mean(
        Relative_abundance,
        na.rm = TRUE
      ),

    .groups = "drop"
  ) %>%

  arrange(
    desc(Mean_abundance)
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

  Sample =
    sample_order,

  Genus =
    top_heatmap_genera
)


heatmap_data <- genus_identified %>%

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

    heatmap_data,

    by = c(
      "Sample",
      "Genus"
    )
  ) %>%

  mutate(

    Relative_abundance =
      replace_na(
        Relative_abundance,
        0
      ),

    Log_abundance =
      log10(
        Relative_abundance +
          0.01
      ),

    Sample = factor(
      Sample,
      levels = sample_order
    ),

    Genus = factor(
      Genus,
      levels =
        rev(
          top_heatmap_genera
        )
    )
  )


# ==============================================================================
# 20. HEATMAP
# ==============================================================================

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

    linewidth = 0.4
  ) +

  scale_fill_gradientn(

    colours =
      viridis(
        100,
        option = "C"
      ),

    name =
      expression(
        log[10] *
          "(relative abundance + 0.01)"
      )
  ) +

  scale_x_discrete(
    drop = FALSE
  ) +

  labs(
    x = NULL,
    y = NULL
  ) +

  theme_minimal(
    base_size = 12
  ) +

  theme(

    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1,
        colour = "black"
      ),

    axis.text.y =
      element_text(
        face = "italic",
        colour = "black"
      ),

    panel.grid =
      element_blank(),

    legend.position =
      "right"
  )


print(p_heatmap)


ggsave(

  file.path(
    OUTPUT_DIR,
    "SILVA144_top20_genus_heatmap.png"
  ),

  p_heatmap,

  width = 13,

  height = 8,

  dpi = 600
)


ggsave(

  file.path(
    OUTPUT_DIR,
    "SILVA144_top20_genus_heatmap.pdf"
  ),

  p_heatmap,

  width = 13,

  height = 8
)


# ==============================================================================
# 21. SAVE ABUNDANCE TABLES
# ==============================================================================

write.csv(

  phylum,

  file.path(
    OUTPUT_DIR,
    "phylum_relative_abundance.csv"
  ),

  row.names = FALSE
)


write.csv(

  genus,

  file.path(
    OUTPUT_DIR,
    "genus_relative_abundance.csv"
  ),

  row.names = FALSE
)


write.csv(

  heatmap_data,

  file.path(
    OUTPUT_DIR,
    "top20_genus_heatmap_data.csv"
  ),

  row.names = FALSE
)


# ==============================================================================
# 22. PREPARE TRUE GENUS-LEVEL NETWORK DATA
#
# IMPORTANT:
# Only ASVs actually assigned to genus are included.
#
# Families such as Pseudomonadaceae will NOT be called genera.
# ==============================================================================

network_genus <- dat %>%

  filter(
    !is.na(Genus)
  ) %>%

  group_by(
    Sample,
    Genus
  ) %>%

  summarise(

    Count =
      sum(
        Count,
        na.rm = TRUE
      ),

    .groups = "drop"
  )


# ==============================================================================
# 23. PREVALENCE FILTER
# ==============================================================================

n_samples_total <-
  length(
    sample_order
  )


prevalence <- network_genus %>%

  group_by(Genus) %>%

  summarise(

    Present =
      sum(
        Count > 0
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


cat("\n========================================\n")
cat("NETWORK FILTERING\n")
cat("========================================\n")

cat(
  "Total samples:",
  n_samples_total,
  "\n"
)

cat(
  "Genera retained:",
  length(
    keep_genera
  ),
  "\n"
)


if (length(keep_genera) < 2) {

  stop(
    "Fewer than two genera remain after prevalence filtering."
  )
}


write.csv(

  prevalence,

  file.path(
    OUTPUT_DIR,
    "network_genus_prevalence.csv"
  ),

  row.names = FALSE
)


# ==============================================================================
# 24. BUILD COMPLETE SAMPLE x GENUS MATRIX
# ==============================================================================

network_grid <- expand_grid(

  Sample =
    sample_order,

  Genus =
    keep_genera
)


network_long <- network_grid %>%

  left_join(

    network_genus,

    by = c(
      "Sample",
      "Genus"
    )
  ) %>%

  mutate(

    Count =
      replace_na(
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


rownames(
  network_matrix
) <- network_sample_names


# ==============================================================================
# 25. CLR TRANSFORMATION
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
# 26. PAIRWISE SPEARMAN ASSOCIATIONS
# ==============================================================================

taxa_names <-
  colnames(
    clr_matrix
  )


taxon_pairs <- combn(

  taxa_names,

  2,

  simplify = FALSE
)


safe_spearman <- function(
    pair
) {

  x <-
    clr_matrix[
      ,
      pair[1]
    ]


  y <-
    clr_matrix[
      ,
      pair[2]
    ]


  if (
    sd(x) == 0 ||
      sd(y) == 0
  ) {

    return(

      tibble(

        Taxon1 =
          pair[1],

        Taxon2 =
          pair[2],

        rho =
          NA_real_,

        p =
          NA_real_
      )
    )
  }


  test <- suppressWarnings(

    cor.test(

      x,

      y,

      method =
        "spearman",

      exact =
        FALSE
    )
  )


  tibble(

    Taxon1 =
      pair[1],

    Taxon2 =
      pair[2],

    rho =
      unname(
        test$estimate
      ),

    p =
      test$p.value
  )
}


edge_tests <- map_dfr(

  taxon_pairs,

  safe_spearman
)


# ==============================================================================
# 27. MULTIPLE-TESTING CORRECTION
# ==============================================================================

edge_tests <- edge_tests %>%

  mutate(

    p_adj =
      p.adjust(
        p,
        method = "BH"
      )
  )


# ==============================================================================
# 28. RETAIN SIGNIFICANT NETWORK EDGES
# ==============================================================================

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

    Association =
      if_else(
        rho > 0,
        "Positive",
        "Negative"
      ),

    Weight =
      abs(rho)
  )


cat(
  "Significant edges retained:",
  nrow(
    network_edges
  ),
  "\n"
)


write.csv(

  edge_tests,

  file.path(
    OUTPUT_DIR,
    "all_pairwise_associations.csv"
  ),

  row.names = FALSE
)


write.csv(

  network_edges,

  file.path(
    OUTPUT_DIR,
    "significant_network_edges.csv"
  ),

  row.names = FALSE
)


if (nrow(network_edges) == 0) {

  stop(
    paste(
      "No network edges passed the thresholds.",
      "Inspect the data or adjust the network thresholds."
    )
  )
}


# ==============================================================================
# 29. BUILD IGRAPH NETWORK
# ==============================================================================

g <- graph_from_data_frame(

  network_edges,

  directed = FALSE
)


cat("\n========================================\n")
cat("NETWORK SUMMARY\n")
cat("========================================\n")

cat(
  "Nodes:",
  vcount(g),
  "\n"
)

cat(
  "Edges:",
  ecount(g),
  "\n"
)

cat(
  "Density:",
  edge_density(
    g,
    loops = FALSE
  ),
  "\n"
)

cat(
  "Mean degree:",
  mean(
    degree(g)
  ),
  "\n"
)


# ==============================================================================
# 30. COMMUNITY DETECTION
# ==============================================================================

communities <- cluster_louvain(

  g,

  weights =
    E(g)$Weight
)


V(g)$Module <-
  membership(
    communities
  )


network_modularity <-
  modularity(
    communities
  )


cat(
  "Modules:",
  length(
    communities
  ),
  "\n"
)

cat(
  "Modularity:",
  network_modularity,
  "\n"
)


# ==============================================================================
# 31. NODE CENTRALITY
# ==============================================================================

node_stats <- tibble(

  Genus =
    V(g)$name,

  Degree =
    degree(g),

  Betweenness =
    betweenness(
      g,
      directed = FALSE,
      normalized = TRUE
    ),

  Eigenvector =
    eigen_centrality(
      g,
      directed = FALSE,
      weights =
        E(g)$Weight
    )$vector,

  Module =
    V(g)$Module
) %>%

  arrange(

    desc(Degree),

    desc(Betweenness)
  )


cat("\nTop network taxa:\n")

print(
  head(
    node_stats,
    20
  ),
  n = 20
)


write.csv(

  node_stats,

  file.path(
    OUTPUT_DIR,
    "network_node_centrality.csv"
  ),

  row.names = FALSE
)


# ==============================================================================
# 32. NETWORK SUMMARY CSV
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

    length(
      keep_genera
    ),

    vcount(g),

    ecount(g),

    edge_density(
      g,
      loops = FALSE
    ),

    mean(
      degree(g)
    ),

    length(
      communities
    ),

    network_modularity
  )
)


write.csv(

  network_summary,

  file.path(
    OUTPUT_DIR,
    "network_summary.csv"
  ),

  row.names = FALSE
)


# ==============================================================================
# 33. NETWORK APPEARANCE
# ==============================================================================

V(g)$size <-

  5 +

  scales::rescale(

    degree(g),

    to = c(
      2,
      12
    )
  )


module_ids <- sort(
  unique(
    V(g)$Module
  )
)


module_colours <- setNames(

  hcl.colors(
    length(module_ids),
    palette = "Dark 3"
  ),

  module_ids
)


V(g)$color <-

  unname(
    module_colours[
      as.character(
        V(g)$Module
      )
    ]
  )


V(g)$frame.color <- NA


E(g)$width <- scales::rescale(

  abs(
    E(g)$rho
  ),

  to = c(
    0.5,
    3
  )
)


E(g)$color <- ifelse(

  E(g)$rho > 0,

  "#2C7FB8",

  "#D95F5F"
)


E(g)$lty <- ifelse(

  E(g)$rho > 0,

  1,

  2
)


# ==============================================================================
# 34. NETWORK LAYOUT
# ==============================================================================

set.seed(123)


network_layout <- layout_with_fr(

  g,

  weights =
    E(g)$Weight
)


cat(
  "\nNetwork layout dimensions:",
  paste(
    dim(network_layout),
    collapse = " x "
  ),
  "\n"
)


# ==============================================================================
# 35. DETERMINE NUMBER OF HUB TAXA TO LABEL
#
# FIXED: no n() inside slice_max()
# ==============================================================================

top_n_nodes <- min(

  20L,

  nrow(
    node_stats
  )
)


top_label_taxa <- node_stats %>%

  slice_max(

    order_by = Degree,

    n = top_n_nodes,

    with_ties = FALSE
  ) %>%

  pull(Genus)


V(g)$label_hubs <- ifelse(

  V(g)$name %in%
    top_label_taxa,

  V(g)$name,

  NA_character_
)


# ==============================================================================
# 36. FUNCTION TO DRAW NETWORK
#
# This avoids repeating long plot() code and prevents syntax errors.
# ==============================================================================

draw_network <- function(
    labels,
    title_text
) {

  plot(

    g,

    layout =
      network_layout,

    vertex.label =
      labels,

    vertex.label.cex =
      0.65,

    vertex.label.color =
      "black",

    vertex.size =
      V(g)$size,

    vertex.color =
      V(g)$color,

    vertex.frame.color =
      NA,

    edge.color =
      E(g)$color,

    edge.width =
      E(g)$width,

    edge.lty =
      E(g)$lty,

    edge.curved =
      0.08,

    main =
      title_text
  )


  legend(

    "topleft",

    legend = c(

      "Positive association",

      "Negative association"
    ),

    col = c(

      "#2C7FB8",

      "#D95F5F"
    ),

    lty = c(
      1,
      2
    ),

    lwd = 2,

    bty = "n"
  )


  legend(

    "topright",

    legend =
      paste(
        "Module",
        module_ids
      ),

    pch = 21,

    pt.bg =
      module_colours,

    pt.cex = 1.4,

    bty = "n"
  )
}


# ==============================================================================
# 37. SAVE FULL NETWORK PNG
# ==============================================================================

png(

  filename =
    file.path(
      OUTPUT_DIR,
      "cooccurrence_network_full.png"
    ),

  width = 4500,

  height = 3800,

  res = 400
)


draw_network(

  labels =
    V(g)$name,

  title_text =
    "Genus-level bacterial co-occurrence network"
)


dev.off()


# ==============================================================================
# 38. SAVE HUB-LABELLED NETWORK PNG
# ==============================================================================

png(

  filename =
    file.path(
      OUTPUT_DIR,
      "cooccurrence_network_hubs_labelled.png"
    ),

  width = 4500,

  height = 3800,

  res = 400
)


draw_network(

  labels =
    V(g)$label_hubs,

  title_text =
    "Bacterial co-occurrence network: highly connected genera"
)


dev.off()


# ==============================================================================
# 39. SAVE HUB-LABELLED NETWORK PDF
# ==============================================================================

pdf(

  file =
    file.path(
      OUTPUT_DIR,
      "cooccurrence_network_hubs_labelled.pdf"
    ),

  width = 14,

  height = 12
)


draw_network(

  labels =
    V(g)$label_hubs,

  title_text =
    "Bacterial co-occurrence network: highly connected genera"
)


dev.off()


# ==============================================================================
# 40. TOP HUB-TAXA BAR PLOT
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

    Module = factor(
      Module
    )
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

  labs(

    x = "Degree",

    y = NULL,

    fill = "Module",

    title =
      "Most connected genera in the co-occurrence network"
  ) +

  theme_classic(
    base_size = 12
  ) +

  theme(

    axis.text.y =
      element_text(
        face = "italic",
        colour = "black"
      ),

    axis.text.x =
      element_text(
        colour = "black"
      ),

    axis.title.x =
      element_text(
        face = "bold"
      ),

    legend.title =
      element_text(
        face = "bold"
      ),

    legend.position =
      "right"
  )


print(p_degree)


ggsave(

  file.path(
    OUTPUT_DIR,
    "top20_network_taxa_degree.png"
  ),

  p_degree,

  width = 10,

  height = 8,

  dpi = 600
)


ggsave(

  file.path(
    OUTPUT_DIR,
    "top20_network_taxa_degree.pdf"
  ),

  p_degree,

  width = 10,

  height = 8
)


# ==============================================================================
# 41. EXPORT NETWORK TO GRAPHML
#
# Can be opened in Cytoscape or Gephi.
# ==============================================================================

write_graph(

  g,

  file.path(
    OUTPUT_DIR,
    "cooccurrence_network.graphml"
  ),

  format = "graphml"
)


# ==============================================================================
# 42. FINAL SAMPLE REPRESENTATION CHECK
# ==============================================================================

sample_check <- tibble(

  Sample =
    sample_order,

  In_feature_table =
    sample_order %in%
      all_samples,

  In_phylum_plot =
    sample_order %in%
      as.character(
        unique(
          phylum_plot_data$Sample
        )
      ),

  In_genus_plot =
    sample_order %in%
      as.character(
        unique(
          genus_plot_data$Sample
        )
      ),

  In_heatmap =
    sample_order %in%
      as.character(
        unique(
          heatmap_data$Sample
        )
      ),

  In_network_input =
    sample_order %in%
      unique(
        network_long$Sample
      )
)


cat("\n========================================\n")
cat("FINAL SAMPLE CHECK\n")
cat("========================================\n")

print(
  sample_check,
  n = Inf
)


write.csv(

  sample_check,

  file.path(
    OUTPUT_DIR,
    "sample_representation_check.csv"
  ),

  row.names = FALSE
)


# ==============================================================================
# 43. CHECK FOR ANY FAILED SAMPLE REPRESENTATION
# ==============================================================================

if (
  any(
    !sample_check$In_phylum_plot
  ) ||
    any(
      !sample_check$In_genus_plot
    ) ||
    any(
      !sample_check$In_heatmap
    )
) {

  warning(
    "At least one sample is missing from one or more plots. Check sample_representation_check.csv."
  )

} else {

  cat(
    "\nAll samples are represented in the composition plots and heatmap.\n"
  )
}


# ==============================================================================
# 44. FINISHED
# ==============================================================================

cat("\n========================================\n")
cat("ANALYSIS COMPLETE\n")
cat("========================================\n")

cat(
  "\nAll outputs saved in:\n",
  OUTPUT_DIR,
  "\n\n"
)

cat(
  "Main figures:\n"
)

cat(
  "1. SILVA144_phylum_relative_abundance.png\n"
)

cat(
  "2. SILVA144_genus_relative_abundance.png\n"
)

cat(
  "3. SILVA144_top20_genus_heatmap.png\n"
)

cat(
  "4. cooccurrence_network_full.png\n"
)

cat(
  "5. cooccurrence_network_hubs_labelled.png\n"
)

cat(
  "6. top20_network_taxa_degree.png\n\n"
)

cat(
  "Network files:\n"
)

cat(
  "7. network_node_centrality.csv\n"
)

cat(
  "8. significant_network_edges.csv\n"
)

cat(
  "9. network_summary.csv\n"
)

cat(
  "10. cooccurrence_network.graphml\n"
)
