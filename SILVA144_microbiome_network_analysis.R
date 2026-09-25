# ==============================================================================
# SILVA 144 TAXONOMIC COMPOSITION AND MICROBIAL ASSOCIATION NETWORK ANALYSIS
# Potato rhizosphere bacterial microbiome
#
# INPUTS
#   R_taxonomy_plots/feature-table.tsv
#   R_taxonomy_plots/taxonomy144/taxonomy.tsv
#
# OUTPUT
#   R_taxonomy_plots/final_outputs/
#
# FIGURES
#   A. Phylum-level relative abundance
#   B. Genus-level relative abundance
#   C. Top 20 genera heatmap
#   D. Global genus-level association network
#   E. Most connected genera in the global network
#   F. Clanwilliam association subnetwork
#   G. Dendron association subnetwork
#   H. Mamusha association subnetwork
#   I. Wesselesbron association subnetwork
#
# IMPORTANT
#   Network edges represent statistical associations.
#   They do not demonstrate direct ecological interaction, cooperation,
#   competition, or causality.
#
#   Farm subnetworks are NOT independently inferred from only five samples.
#   They display the portion of the GLOBAL statistically supported network
#   represented by taxa occurring at each farm.
# ==============================================================================


# ==============================================================================
# 1. USER SETTINGS
# ==============================================================================

FEATURE_FILE <-
  "R_taxonomy_plots/feature-table.tsv"

TAXONOMY_FILE <-
  "R_taxonomy_plots/taxonomy144/taxonomy.tsv"

OUTPUT_DIR <-
  "R_taxonomy_plots/final_outputs"


FIGURE_DIR <-
  file.path(
    OUTPUT_DIR,
    "figures"
  )

TABLE_DIR <-
  file.path(
    OUTPUT_DIR,
    "tables"
  )

NETWORK_DIR <-
  file.path(
    OUTPUT_DIR,
    "networks"
  )

FARM_NETWORK_DIR <-
  file.path(
    NETWORK_DIR,
    "farm_subnetworks"
  )


# ------------------------------------------------------------------------------
# Taxonomic figures
# ------------------------------------------------------------------------------

TOP_PHYLA <- 10

TOP_GENERA <- 15

TOP_HEATMAP_GENERA <- 20


# ------------------------------------------------------------------------------
# Global network parameters
# ------------------------------------------------------------------------------

PREVALENCE_THRESHOLD <- 0.20

CORRELATION_THRESHOLD <- 0.60

FDR_THRESHOLD <- 0.05

PSEUDOCOUNT <- 0.5

CORRELATION_METHOD <- "spearman"


# ------------------------------------------------------------------------------
# Farm subnetworks
#
# A genus must occur in at least this many samples at the farm
# before being displayed in that farm's subnetwork.
# ------------------------------------------------------------------------------

FARM_MIN_SAMPLES <- 2


# ------------------------------------------------------------------------------
# Filtering
# ------------------------------------------------------------------------------

REMOVE_NON_TARGET <- TRUE


# ------------------------------------------------------------------------------
# Farm order
# ------------------------------------------------------------------------------

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

  A =
    "A. Phylum-level relative abundance",

  B =
    "B. Genus-level relative abundance",

  C =
    "C. Top 20 genera heatmap",

  D =
    "D. Global genus-level association network",

  E =
    "E. Most connected genera in the global network",

  F =
    "F. Clanwilliam association subnetwork",

  G =
    "G. Dendron association subnetwork",

  H =
    "H. Mamusha association subnetwork",

  I =
    "I. Wesselesbron association subnetwork"
)


FARM_FIGURE_CODES <- c(

  Clanwilliam =
    "F",

  Dendron =
    "G",

  Mamusha =
    "H",

  Wesselesbron =
    "I"
)


# ==============================================================================
# 3. STANDARD COLOURS
# ==============================================================================

# Positive / negative network edges
EDGE_POSITIVE <- "#0072B2"

EDGE_NEGATIVE <- "#D55E00"


# Neutral categories
COLOUR_OTHER <- "#E4E4E4"

COLOUR_UNCLASSIFIED <- "#9E9E9E"


# Farm colours
FARM_COLOURS <- c(

  Clanwilliam =
    "#0072B2",

  Dendron =
    "#009E73",

  Mamusha =
    "#D55E00",

  Wesselesbron =
    "#CC79A7"
)


# ==============================================================================
# 4. REQUIRED PACKAGES
# ==============================================================================

required_packages <- c(
  "tidyverse",
  "igraph",
  "scales",
  "viridisLite"
)


missing_packages <-
  required_packages[
    !required_packages %in%
      rownames(
        installed.packages()
      )
  ]


if (
  length(
    missing_packages
  ) > 0
) {

  stop(
    paste(
      "Install the following R packages first:",
      paste(
        missing_packages,
        collapse = ", "
      )
    )
  )
}


library(tidyverse)
library(igraph)
library(scales)
library(viridisLite)


# ==============================================================================
# 5. CREATE OUTPUT DIRECTORIES
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
# 6. STANDARD THEMES
# ==============================================================================

theme_microbiome <- function() {

  theme_classic(
    base_size = 12
  ) +

    theme(

      plot.title =
        element_text(
          face = "bold",
          size = 14,
          hjust = 0
        ),

      plot.subtitle =
        element_text(
          size = 10.5,
          colour = "grey30",
          margin = margin(
            b = 10
          )
        ),

      axis.title =
        element_text(
          face = "bold"
        ),

      axis.text =
        element_text(
          colour = "black"
        ),

      strip.background =
        element_rect(
          fill = "#F2F2F2",
          colour = NA
        ),

      strip.text =
        element_text(
          face = "bold",
          colour = "black"
        ),

      legend.title =
        element_text(
          face = "bold"
        ),

      legend.position =
        "right",

      plot.margin =
        margin(
          10,
          15,
          10,
          10
        )
    )
}


theme_heatmap <- function() {

  theme_minimal(
    base_size = 12
  ) +

    theme(

      plot.title =
        element_text(
          face = "bold",
          size = 14,
          hjust = 0
        ),

      plot.subtitle =
        element_text(
          size = 10.5,
          colour = "grey30",
          margin = margin(
            b = 10
          )
        ),

      panel.grid =
        element_blank(),

      strip.background =
        element_rect(
          fill = "#F2F2F2",
          colour = NA
        ),

      strip.text =
        element_text(
          face = "bold",
          colour = "black"
        ),

      axis.text =
        element_text(
          colour = "black"
        ),

      legend.title =
        element_text(
          face = "bold"
        )
    )
}


# ==============================================================================
# 7. HELPER FUNCTIONS
# ==============================================================================


# ------------------------------------------------------------------------------
# Farm from sample ID
# ------------------------------------------------------------------------------

farm_from_sample <- function(x) {

  case_when(

    str_detect(
      x,
      "^ND-CL"
    ) ~
      "Clanwilliam",

    str_detect(
      x,
      "^ND-DEN"
    ) ~
      "Dendron",

    str_detect(
      x,
      "^ND-MA"
    ) ~
      "Mamusha",

    str_detect(
      x,
      "^ND-WB"
    ) ~
      "Wesselesbron",

    TRUE ~
      "Unknown"
  )
}


# ------------------------------------------------------------------------------
# Taxonomy parser
# ------------------------------------------------------------------------------

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


  modern <-
    str_match(
      x,
      modern_pattern
    )[, 2]


  old <-
    str_match(
      x,
      old_pattern
    )[, 2]


  result <- ifelse(
    !is.na(modern),
    modern,
    old
  )


  result <-
    str_trim(
      result
    )


  result[
    result == "" |
      result == "Unassigned"
  ] <-
    NA


  result
}


# ------------------------------------------------------------------------------
# Colour palette for taxonomic stacked bars
# ------------------------------------------------------------------------------

make_taxon_palette <- function(
    categories
) {

  main_taxa <- setdiff(
    categories,
    c(
      "Other",
      "Unclassified"
    )
  )


  colours <- setNames(

    hcl.colors(
      length(main_taxa),
      palette = "Dynamic"
    ),

    main_taxa
  )


  if (
    "Unclassified" %in%
      categories
  ) {

    colours <- c(
      colours,
      Unclassified =
        COLOUR_UNCLASSIFIED
    )
  }


  if (
    "Other" %in%
      categories
  ) {

    colours <- c(
      colours,
      Other =
        COLOUR_OTHER
    )
  }


  colours
}


# ------------------------------------------------------------------------------
# Safe scaling
# ------------------------------------------------------------------------------

safe_rescale <- function(
    x,
    to,
    constant_value = mean(to)
) {

  if (
    length(x) == 0
  ) {

    return(
      numeric(0)
    )
  }


  if (
    length(
      unique(x)
    ) <= 1
  ) {

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
# Save ggplot consistently
# ------------------------------------------------------------------------------

save_figure <- function(
    plot_object,
    file_stub,
    width,
    height
) {

  ggsave(

    filename =
      file.path(
        FIGURE_DIR,
        paste0(
          file_stub,
          ".png"
        )
      ),

    plot =
      plot_object,

    width =
      width,

    height =
      height,

    dpi =
      600
  )


  ggsave(

    filename =
      file.path(
        FIGURE_DIR,
        paste0(
          file_stub,
          ".pdf"
        )
      ),

    plot =
      plot_object,

    width =
      width,

    height =
      height
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


colnames(
  asv
)[1] <-
  "FeatureID"


sample_order <-
  colnames(
    asv
  )[-1]


sample_metadata <- tibble(

  Sample =
    sample_order
) %>%

  mutate(

    Farm =
      farm_from_sample(
        Sample
      ),

    Farm =
      factor(
        Farm,
        levels = FARM_ORDER
      ),

    Sample =
      factor(
        Sample,
        levels = sample_order
      )
  )


cat(
  "\nSample metadata:\n"
)

print(
  sample_metadata,
  n = Inf
)


if (
  any(
    is.na(
      sample_metadata$Farm
    )
  )
) {

  warning(
    "At least one sample could not be assigned to a farm."
  )
}


# ==============================================================================
# 9. IMPORT SILVA 144 TAXONOMY
# ==============================================================================

tax <- read.delim(

  TAXONOMY_FILE,

  check.names = FALSE,

  quote = "",

  comment.char = ""
)


colnames(
  tax
)[1] <-
  "FeatureID"


if (
  !"Taxon" %in%
    colnames(
      tax
    )
) {

  stop(
    "No 'Taxon' column found in taxonomy file."
  )
}


tax_clean <- tax %>%

  transmute(

    FeatureID,

    FullTaxonomy =
      Taxon,

    Domain =
      extract_rank(
        Taxon,
        "(?:d|k)",
        0
      ),

    Phylum =
      extract_rank(
        Taxon,
        "p",
        1
      ),

    Class =
      extract_rank(
        Taxon,
        "c",
        2
      ),

    Order =
      extract_rank(
        Taxon,
        "o",
        3
      ),

    Family =
      extract_rank(
        Taxon,
        "f",
        4
      ),

    Genus =
      extract_rank(
        Taxon,
        "g",
        5
      ),

    Species =
      extract_rank(
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
# 10. FEATURE TABLE TO LONG FORMAT
# ==============================================================================

asv_long <- asv %>%

  pivot_longer(

    cols =
      -FeatureID,

    names_to =
      "Sample",

    values_to =
      "Count"
  ) %>%

  mutate(

    Count =
      as.numeric(
        Count
      ),

    Sample =
      as.character(
        Sample
      )
  )


dat <- asv_long %>%

  left_join(
    tax_clean,
    by = "FeatureID"
  ) %>%

  left_join(

    sample_metadata %>%
      mutate(
        Sample =
          as.character(
            Sample
          )
      ),

    by = "Sample"
  )


# ==============================================================================
# 11. REMOVE NON-TARGET SEQUENCES
# ==============================================================================

if (
  REMOVE_NON_TARGET
) {

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

      is.na(
        Domain
      ) |
        Domain ==
        "Bacteria"
    )
}


# ==============================================================================
# 12. TAXONOMIC RESOLUTION SUMMARY
# ==============================================================================

taxonomy_resolution <- tax_clean %>%

  summarise(

    Total_ASVs =
      n(),

    Phylum_resolved =
      sum(
        !is.na(
          Phylum
        )
      ),

    Family_resolved =
      sum(
        !is.na(
          Family
        )
      ),

    Genus_resolved =
      sum(
        !is.na(
          Genus
        )
      ),

    Species_resolved =
      sum(
        !is.na(
          Species
        )
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
# 13. FUNCTION: RELATIVE ABUNDANCE
# ==============================================================================

make_rank_abundance <- function(
    data,
    rank_name
) {

  data %>%

    mutate(

      Taxon =
        coalesce(
          .data[
            [rank_name]
          ],
          "Unclassified"
        )
    ) %>%

    group_by(
      Sample,
      Farm,
      Taxon
    ) %>%

    summarise(

      Count =
        sum(
          Count,
          na.rm = TRUE
        ),

      .groups =
        "drop"
    ) %>%

    group_by(
      Sample
    ) %>%

    mutate(

      Total =
        sum(
          Count
        ),

      Relative_abundance =
        ifelse(
          Total > 0,
          Count /
            Total *
            100,
          0
        )
    ) %>%

    ungroup() %>%

    select(
      -Total
    )
}


# ==============================================================================
# 14. FUNCTION: COLLAPSE TOP TAXA
# ==============================================================================

collapse_top_taxa <- function(
    abundance_data,
    n_top
) {

  summary_table <- abundance_data %>%

    group_by(
      Taxon
    ) %>%

    summarise(

      Mean_abundance =
        mean(
          Relative_abundance,
          na.rm = TRUE
        ),

      .groups =
        "drop"
    ) %>%

    arrange(
      desc(
        Mean_abundance
      )
    )


  n_keep <- min(
    n_top,
    nrow(
      summary_table
    )
  )


  top_taxa <- summary_table %>%

    slice_head(
      n = n_keep
    ) %>%

    pull(
      Taxon
    )


  collapsed <- abundance_data %>%

    mutate(

      Taxon_plot =
        if_else(
          Taxon %in%
            top_taxa,
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

      Relative_abundance =
        sum(
          Relative_abundance,
          na.rm = TRUE
        ),

      .groups =
        "drop"
    )


  categories <- unique(
    c(
      top_taxa,
      if (
        "Other" %in%
          collapsed$Taxon_plot
      ) {
        "Other"
      }
    )
  )


  complete_grid <- expand_grid(

    Sample =
      as.character(
        sample_metadata$Sample
      ),

    Taxon_plot =
      categories
  ) %>%

    left_join(

      sample_metadata %>%
        mutate(
          Sample =
            as.character(
              Sample
            )
        ),

      by =
        "Sample"
    )


  collapsed <- complete_grid %>%

    left_join(

      collapsed,

      by =
        c(
          "Sample",
          "Farm",
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

    group_by(
      Sample
    ) %>%

    mutate(

      Total =
        sum(
          Relative_abundance
        ),

      Relative_abundance =
        ifelse(
          Total > 0,
          Relative_abundance /
            Total *
            100,
          0
        )
    ) %>%

    ungroup() %>%

    select(
      -Total
    ) %>%

    mutate(

      Sample =
        factor(
          Sample,
          levels =
            sample_order
        ),

      Farm =
        factor(
          Farm,
          levels =
            FARM_ORDER
        ),

      Taxon_plot =
        factor(
          Taxon_plot,
          levels =
            categories
        )
    )


  list(

    data =
      collapsed,

    categories =
      categories,

    top_taxa =
      top_taxa
  )
}


# ==============================================================================
# 15. FIGURE A — PHYLUM RELATIVE ABUNDANCE
# ==============================================================================

phylum <-
  make_rank_abundance(
    dat,
    "Phylum"
  )


phylum_result <-
  collapse_top_taxa(
    phylum,
    TOP_PHYLA
  )


phylum_plot_data <-
  phylum_result$data


phylum_colours <-
  make_taxon_palette(
    phylum_result$categories
  )


p_phylum <- ggplot(

  phylum_plot_data,

  aes(

    x =
      Sample,

    y =
      Relative_abundance,

    fill =
      Taxon_plot
  )
) +

  geom_col(
    width = 0.86
  ) +

  facet_grid(

    . ~ Farm,

    scales =
      "free_x",

    space =
      "free_x"
  ) +

  scale_fill_manual(

    values =
      phylum_colours,

    drop =
      FALSE
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
    ylim =
      c(
        0,
        100
      )
  ) +

  labs(

    title =
      FIGURE_TITLES[
        ["A"]
      ],

    subtitle =
      "SILVA 144 taxonomy; samples grouped by farm",

    x =
      NULL,

    y =
      "Relative abundance (%)",

    fill =
      "Phylum"
  ) +

  theme_microbiome() +

  theme(

    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      )
  )


print(
  p_phylum
)


save_figure(

  p_phylum,

  "Fig_A_phylum_relative_abundance",

  14,

  7
)


# ==============================================================================
# 16. FIGURE B — GENUS RELATIVE ABUNDANCE
# ==============================================================================

genus <-
  make_rank_abundance(
    dat,
    "Genus"
  )


genus_result <-
  collapse_top_taxa(
    genus,
    TOP_GENERA
  )


genus_plot_data <-
  genus_result$data


genus_colours <-
  make_taxon_palette(
    genus_result$categories
  )


p_genus <- ggplot(

  genus_plot_data,

  aes(

    x =
      Sample,

    y =
      Relative_abundance,

    fill =
      Taxon_plot
  )
) +

  geom_col(
    width = 0.86
  ) +

  facet_grid(

    . ~ Farm,

    scales =
      "free_x",

    space =
      "free_x"
  ) +

  scale_fill_manual(

    values =
      genus_colours,

    drop =
      FALSE
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
    ylim =
      c(
        0,
        100
      )
  ) +

  labs(

    title =
      FIGURE_TITLES[
        ["B"]
      ],

    subtitle =
      "SILVA 144 taxonomy; samples grouped by farm",

    x =
      NULL,

    y =
      "Relative abundance (%)",

    fill =
      "Genus"
  ) +

  theme_microbiome() +

  theme(

    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),

    legend.text =
      element_text(
        face = "italic"
      )
  )


print(
  p_genus
)


save_figure(

  p_genus,

  "Fig_B_genus_relative_abundance",

  15,

  7
)


# ==============================================================================
# 17. FIGURE C — TOP 20 GENERA HEATMAP
# ==============================================================================

identified_genus <- genus %>%

  filter(
    Taxon !=
      "Unclassified"
  )


genus_summary <- identified_genus %>%

  group_by(
    Taxon
  ) %>%

  summarise(

    Mean_abundance =
      mean(
        Relative_abundance,
        na.rm = TRUE
      ),

    .groups =
      "drop"
  ) %>%

  arrange(
    desc(
      Mean_abundance
    )
  )


n_heatmap <- min(

  TOP_HEATMAP_GENERA,

  nrow(
    genus_summary
  )
)


top_heatmap_genera <-
  genus_summary %>%

  slice_head(
    n =
      n_heatmap
  ) %>%

  pull(
    Taxon
  )


heatmap_grid <- expand_grid(

  Sample =
    as.character(
      sample_metadata$Sample
    ),

  Genus =
    top_heatmap_genera
) %>%

  left_join(

    sample_metadata %>%
      mutate(
        Sample =
          as.character(
            Sample
          )
      ),

    by =
      "Sample"
  )


heatmap_data <-
  identified_genus %>%

  filter(
    Taxon %in%
      top_heatmap_genera
  ) %>%

  transmute(

    Sample,

    Genus =
      Taxon,

    Relative_abundance
  )


heatmap_data <- heatmap_grid %>%

  left_join(

    heatmap_data,

    by =
      c(
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

    Sample =
      factor(
        Sample,
        levels =
          sample_order
      ),

    Farm =
      factor(
        Farm,
        levels =
          FARM_ORDER
      ),

    Genus =
      factor(
        Genus,
        levels =
          rev(
            top_heatmap_genera
          )
      )
  )


p_heatmap <- ggplot(

  heatmap_data,

  aes(

    x =
      Sample,

    y =
      Genus,

    fill =
      Log_abundance
  )
) +

  geom_tile(

    colour =
      "white",

    linewidth =
      0.35
  ) +

  facet_grid(

    . ~ Farm,

    scales =
      "free_x",

    space =
      "free_x"
  ) +

  scale_fill_gradientn(

    colours =
      viridis(
        100,
        option = "D"
      ),

    name =
      expression(
        log[10] *
          "(relative abundance + 0.01)"
      )
  ) +

  labs(

    title =
      FIGURE_TITLES[
        ["C"]
      ],

    subtitle =
      "Top genera ranked by mean relative abundance",

    x =
      NULL,

    y =
      NULL
  ) +

  theme_heatmap() +

  theme(

    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      ),

    axis.text.y =
      element_text(
        face = "italic"
      )
  )


print(
  p_heatmap
)


save_figure(

  p_heatmap,

  "Fig_C_top20_genera_heatmap",

  14,

  8
)


# ==============================================================================
# 18. SAVE RELATIVE ABUNDANCE TABLES
# ==============================================================================

write.csv(

  phylum,

  file.path(
    TABLE_DIR,
    "phylum_relative_abundance.csv"
  ),

  row.names =
    FALSE
)


write.csv(

  genus,

  file.path(
    TABLE_DIR,
    "genus_relative_abundance.csv"
  ),

  row.names =
    FALSE
)


# ==============================================================================
# 19. PREPARE TRUE GENUS-LEVEL NETWORK DATA
# ==============================================================================

network_genus <- dat %>%

  filter(
    !is.na(
      Genus
    )
  ) %>%

  group_by(
    Farm,
    Sample,
    Genus
  ) %>%

  summarise(

    Count =
      sum(
        Count,
        na.rm = TRUE
      ),

    .groups =
      "drop"
  )


# ==============================================================================
# 20. GLOBAL PREVALENCE FILTER
# ==============================================================================

n_samples_total <-
  length(
    sample_order
  )


prevalence <- network_genus %>%

  group_by(
    Genus
  ) %>%

  summarise(

    Present =
      n_distinct(
        Sample[
          Count > 0
        ]
      ),

    Prevalence =
      Present /
        n_samples_total,

    .groups =
      "drop"
  )


keep_genera <- prevalence %>%

  filter(
    Prevalence >=
      PREVALENCE_THRESHOLD
  ) %>%

  pull(
    Genus
  )


cat(
  "\nGenera retained for global network:",
  length(
    keep_genera
  ),
  "\n"
)


if (
  length(
    keep_genera
  ) < 2
) {

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

  row.names =
    FALSE
)


# ==============================================================================
# 21. SAMPLE × GENUS MATRIX
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

    Count =
      sum(
        Count
      ),

    .groups =
      "drop"
  )


network_grid <- expand_grid(

  Sample =
    sample_order,

  Genus =
    keep_genera
)


network_long <- network_grid %>%

  left_join(

    network_counts,

    by =
      c(
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


network_matrix_df <-
  network_long %>%

  pivot_wider(

    names_from =
      Genus,

    values_from =
      Count,

    values_fill =
      0
  )


network_sample_names <-
  network_matrix_df$Sample


network_matrix <-
  network_matrix_df %>%

  select(
    -Sample
  ) %>%

  as.matrix()


rownames(
  network_matrix
) <-
  network_sample_names


# ==============================================================================
# 22. CLR TRANSFORMATION
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
# 23. PAIRWISE ASSOCIATIONS
# ==============================================================================

taxa_names <-
  colnames(
    clr_matrix
  )


taxon_pairs <- combn(

  taxa_names,

  2,

  simplify =
    FALSE
)


safe_correlation <- function(
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


  result <- suppressWarnings(

    cor.test(

      x,

      y,

      method =
        CORRELATION_METHOD,

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
        result$estimate
      ),

    p =
      result$p.value
  )
}


edge_tests <- map_dfr(

  taxon_pairs,

  safe_correlation
)


edge_tests <- edge_tests %>%

  mutate(

    p_adj =
      p.adjust(
        p,
        method =
          "BH"
      )
  )


network_edges <- edge_tests %>%

  filter(

    !is.na(
      rho
    ),

    !is.na(
      p_adj
    ),

    abs(
      rho
    ) >=
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
      abs(
        rho
      )
  )


write.csv(

  edge_tests,

  file.path(
    TABLE_DIR,
    "all_pairwise_associations.csv"
  ),

  row.names =
    FALSE
)


write.csv(

  network_edges,

  file.path(
    TABLE_DIR,
    "significant_global_network_edges.csv"
  ),

  row.names =
    FALSE
)


if (
  nrow(
    network_edges
  ) == 0
) {

  stop(
    "No associations passed the selected network thresholds."
  )
}


# ==============================================================================
# 24. BUILD GLOBAL IGRAPH NETWORK
# ==============================================================================

g <- graph_from_data_frame(

  network_edges,

  directed =
    FALSE
)


# ==============================================================================
# 25. GLOBAL COMMUNITY DETECTION
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


# ==============================================================================
# 26. GLOBAL NODE STATISTICS
# ==============================================================================

node_stats <- tibble(

  Genus =
    V(g)$name,

  Degree =
    degree(
      g
    ),

  Betweenness =
    betweenness(

      g,

      directed =
        FALSE,

      normalized =
        TRUE
    ),

  Eigenvector =
    eigen_centrality(

      g,

      directed =
        FALSE,

      weights =
        E(g)$Weight
    )$vector,

  Module =
    V(g)$Module
) %>%

  arrange(

    desc(
      Degree
    ),

    desc(
      Betweenness
    )
  )


write.csv(

  node_stats,

  file.path(
    TABLE_DIR,
    "global_network_node_centrality.csv"
  ),

  row.names =
    FALSE
)


network_summary <- tibble(

  Metric =
    c(
      "Samples",
      "Genera after prevalence filtering",
      "Network nodes",
      "Network edges",
      "Network density",
      "Mean degree",
      "Modules",
      "Modularity"
    ),

  Value =
    c(
      n_samples_total,
      length(
        keep_genera
      ),
      vcount(
        g
      ),
      ecount(
        g
      ),
      edge_density(
        g,
        loops = FALSE
      ),
      mean(
        degree(
          g
        )
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
    TABLE_DIR,
    "global_network_summary.csv"
  ),

  row.names =
    FALSE
)


# ==============================================================================
# 27. GLOBAL NETWORK COLOURS
# ==============================================================================

module_ids <-
  sort(
    unique(
      V(g)$Module
    )
  )


MODULE_COLOURS <- setNames(

  hcl.colors(
    length(
      module_ids
    ),
    palette =
      "Set 2"
  ),

  as.character(
    module_ids
  )
)


V(g)$color <-
  unname(

    MODULE_COLOURS[
      as.character(
        V(g)$Module
      )
    ]
  )


V(g)$frame.color <-
  NA


V(g)$size <-

  5 +

  safe_rescale(

    degree(
      g
    ),

    to =
      c(
        2,
        12
      ),

    constant_value =
      7
  )


E(g)$color <-
  ifelse(

    E(g)$rho > 0,

    EDGE_POSITIVE,

    EDGE_NEGATIVE
  )


E(g)$width <-
  safe_rescale(

    abs(
      E(g)$rho
    ),

    to =
      c(
        0.6,
        3
      ),

    constant_value =
      1.5
  )


E(g)$lty <-
  ifelse(

    E(g)$rho > 0,

    1,

    2
  )


# ==============================================================================
# 28. GLOBAL NETWORK LAYOUT
# ==============================================================================

set.seed(
  123
)


global_layout <- layout_with_fr(

  g,

  weights =
    E(g)$Weight
)


top_n_nodes <- min(

  20L,

  nrow(
    node_stats
  )
)


global_label_taxa <- node_stats %>%

  slice_max(

    order_by =
      Degree,

    n =
      top_n_nodes,

    with_ties =
      FALSE
  ) %>%

  pull(
    Genus
  )


V(g)$label_hubs <-
  ifelse(

    V(g)$name %in%
      global_label_taxa,

    V(g)$name,

    NA_character_
  )


# ==============================================================================
# 29. STANDARD NETWORK DRAWING FUNCTION
# ==============================================================================

draw_network <- function(

    graph_object,

    layout_object,

    labels,

    title_text,

    subtitle_text

) {

  plot(

    graph_object,

    layout =
      layout_object,

    vertex.label =
      labels,

    vertex.label.cex =
      0.68,

    vertex.label.font =
      3,

    vertex.label.color =
      "black",

    vertex.size =
      V(
        graph_object
      )$size,

    vertex.color =
      V(
        graph_object
      )$color,

    vertex.frame.color =
      NA,

    edge.color =
      E(
        graph_object
      )$color,

    edge.width =
      E(
        graph_object
      )$width,

    edge.lty =
      E(
        graph_object
      )$lty,

    edge.curved =
      0.08,

    main =
      ""
  )


  title(

    main =
      title_text,

    sub =
      subtitle_text,

    cex.main =
      1.25,

    font.main =
      2,

    cex.sub =
      0.85
  )


  legend(

    "topleft",

    legend =
      c(
        "Positive association",
        "Negative association"
      ),

    col =
      c(
        EDGE_POSITIVE,
        EDGE_NEGATIVE
      ),

    lty =
      c(
        1,
        2
      ),

    lwd =
      2,

    bty =
      "n",

    cex =
      0.85
  )


  modules_present <-
    sort(
      unique(
        V(
          graph_object
        )$Module
      )
    )


  legend(

    "topright",

    legend =
      paste(
        "Module",
        modules_present
      ),

    pch =
      21,

    pt.bg =
      MODULE_COLOURS[
        as.character(
          modules_present
        )
      ],

    pt.cex =
      1.4,

    bty =
      "n",

    cex =
      0.8
  )
}


# ==============================================================================
# 30. FIGURE D — GLOBAL NETWORK
# ==============================================================================

png(

  filename =
    file.path(
      FIGURE_DIR,
      "Fig_D_global_association_network.png"
    ),

  width =
    4500,

  height =
    3800,

  res =
    400
)


draw_network(

  g,

  global_layout,

  V(g)$label_hubs,

  FIGURE_TITLES[
    ["D"]
  ],

  paste0(
    "|rho| \u2265 ",
    CORRELATION_THRESHOLD,
    "; FDR < ",
    FDR_THRESHOLD,
    "; node colour = module"
  )
)


dev.off()


pdf(

  file =
    file.path(
      FIGURE_DIR,
      "Fig_D_global_association_network.pdf"
    ),

  width =
    14,

  height =
    12
)


draw_network(

  g,

  global_layout,

  V(g)$label_hubs,

  FIGURE_TITLES[
    ["D"]
  ],

  paste0(
    "|rho| \u2265 ",
    CORRELATION_THRESHOLD,
    "; FDR < ",
    FDR_THRESHOLD,
    "; node colour = module"
  )
)


dev.off()


# ==============================================================================
# 31. FIGURE E — GLOBAL HUB TAXA
# ==============================================================================

top_nodes <- node_stats %>%

  slice_max(

    order_by =
      Degree,

    n =
      top_n_nodes,

    with_ties =
      FALSE
  ) %>%

  arrange(
    Degree
  ) %>%

  mutate(

    Genus =
      factor(
        Genus,
        levels =
          Genus
      ),

    Module =
      factor(
        Module
      )
  )


p_degree <- ggplot(

  top_nodes,

  aes(

    x =
      Degree,

    y =
      Genus,

    fill =
      Module
  )
) +

  geom_col(
    width =
      0.72
  ) +

  scale_fill_manual(

    values =
      MODULE_COLOURS
  ) +

  labs(

    title =
      FIGURE_TITLES[
        ["E"]
      ],

    subtitle =
      "Node degree represents the number of retained associations",

    x =
      "Degree",

    y =
      NULL,

    fill =
      "Module"
  ) +

  theme_microbiome() +

  theme(

    axis.text.y =
      element_text(
        face =
          "italic"
      )
  )


print(
  p_degree
)


save_figure(

  p_degree,

  "Fig_E_top20_global_network_genera",

  10,

  8
)


# ==============================================================================
# 32. EXPORT GLOBAL GRAPHML
# ==============================================================================

write_graph(

  g,

  file.path(
    NETWORK_DIR,
    "global_association_network.graphml"
  ),

  format =
    "graphml"
)


# ==============================================================================
# 33. FARM-SPECIFIC SUBNETWORK FUNCTION
#
# These graphs inherit edges and global module membership from g.
# Correlations are NOT recalculated within five-sample farms.
# ==============================================================================

make_farm_subnetwork <- function(
    farm_name
) {

  figure_code <-
    FARM_FIGURE_CODES[
      [farm_name]
    ]


  taxa_here <- network_genus %>%

    filter(

      Farm ==
        farm_name,

      Count >
        0,

      Genus %in%
        V(g)$name
    ) %>%

    group_by(
      Genus
    ) %>%

    summarise(

      Samples_present =
        n_distinct(
          Sample
        ),

      .groups =
        "drop"
    ) %>%

    filter(

      Samples_present >=
        FARM_MIN_SAMPLES
    ) %>%

    pull(
      Genus
    )


  taxa_here <-
    intersect(
      taxa_here,
      V(g)$name
    )


  if (
    length(
      taxa_here
    ) < 2
  ) {

    warning(
      paste(
        farm_name,
        "has fewer than two eligible taxa."
      )
    )

    return(
      NULL
    )
  }


  g_farm <- induced_subgraph(

    g,

    vids =
      taxa_here
  )


  if (
    ecount(
      g_farm
    ) == 0
  ) {

    warning(
      paste(
        farm_name,
        "has no globally supported edges after filtering."
      )
    )

    return(
      NULL
    )
  }


  farm_degree <-
    degree(
      g_farm
    )


  V(g_farm)$size <-

    5 +

    safe_rescale(

      farm_degree,

      to =
        c(
          2,
          12
        ),

      constant_value =
        7
    )


  V(g_farm)$color <-

    MODULE_COLOURS[
      as.character(
        V(g_farm)$Module
      )
    ]


  V(g_farm)$frame.color <-
    NA


  E(g_farm)$color <-
    ifelse(

      E(g_farm)$rho > 0,

      EDGE_POSITIVE,

      EDGE_NEGATIVE
    )


  E(g_farm)$width <-
    safe_rescale(

      abs(
        E(g_farm)$rho
      ),

      to =
        c(
          0.6,
          3
        ),

      constant_value =
        1.5
    )


  E(g_farm)$lty <-
    ifelse(

      E(g_farm)$rho > 0,

      1,

      2
    )


  set.seed(
    123
  )


  farm_layout <-
    layout_with_fr(

      g_farm,

      weights =
        E(g_farm)$Weight
    )


  n_farm_labels <-
    min(

      15L,

      vcount(
        g_farm
      )
    )


  label_taxa <-
    names(

      sort(

        farm_degree,

        decreasing =
          TRUE
      )
    )[
      seq_len(
        n_farm_labels
      )
    ]


  farm_labels <-
    ifelse(

      V(g_farm)$name %in%
        label_taxa,

      V(g_farm)$name,

      NA_character_
    )


  figure_title <-
    FIGURE_TITLES[
      [figure_code]
    ]


  figure_subtitle <-
    paste0(
      "Edges inherited from the global network; genera detected in \u2265 ",
      FARM_MIN_SAMPLES,
      " farm samples"
    )


  file_stub <-
    paste0(
      "Fig_",
      figure_code,
      "_",
      gsub(
        "[^A-Za-z0-9]+",
        "_",
        farm_name
      ),
      "_association_subnetwork"
    )


  # ---------------------------------------------------------------------------
  # PNG
  # ---------------------------------------------------------------------------

  png(

    filename =
      file.path(
        FIGURE_DIR,
        paste0(
          file_stub,
          ".png"
        )
      ),

    width =
      3800,

    height =
      3400,

    res =
      400
  )


  draw_network(

    g_farm,

    farm_layout,

    farm_labels,

    figure_title,

    figure_subtitle
  )


  dev.off()


  # ---------------------------------------------------------------------------
  # PDF
  # ---------------------------------------------------------------------------

  pdf(

    file =
      file.path(
        FIGURE_DIR,
        paste0(
          file_stub,
          ".pdf"
        )
      ),

    width =
      12,

    height =
      11
  )


  draw_network(

    g_farm,

    farm_layout,

    farm_labels,

    figure_title,

    figure_subtitle
  )


  dev.off()


  # ---------------------------------------------------------------------------
  # Node statistics
  # ---------------------------------------------------------------------------

  farm_nodes <- tibble(

    Farm =
      farm_name,

    Genus =
      V(g_farm)$name,

    Degree =
      degree(
        g_farm
      ),

    Betweenness =
      betweenness(

        g_farm,

        directed =
          FALSE,

        normalized =
          TRUE
      ),

    Module =
      V(g_farm)$Module
  ) %>%

    arrange(

      desc(
        Degree
      ),

      desc(
        Betweenness
      )
    )


  write.csv(

    farm_nodes,

    file.path(

      TABLE_DIR,

      paste0(
        "network_nodes_",
        gsub(
          "[^A-Za-z0-9]+",
          "_",
          farm_name
        ),
        ".csv"
      )
    ),

    row.names =
      FALSE
  )


  # ---------------------------------------------------------------------------
  # GraphML
  # ---------------------------------------------------------------------------

  write_graph(

    g_farm,

    file.path(

      FARM_NETWORK_DIR,

      paste0(
        gsub(
          "[^A-Za-z0-9]+",
          "_",
          farm_name
        ),
        "_association_subnetwork.graphml"
      )
    ),

    format =
      "graphml"
  )


  # ---------------------------------------------------------------------------
  # Summary
  # ---------------------------------------------------------------------------

  tibble(

    Farm =
      farm_name,

    Samples_required =
      FARM_MIN_SAMPLES,

    Nodes =
      vcount(
        g_farm
      ),

    Edges =
      ecount(
        g_farm
      ),

    Density =
      edge_density(
        g_farm,
        loops =
          FALSE
      ),

    Mean_degree =
      mean(
        degree(
          g_farm
        )
      )
  )
}


# ==============================================================================
# 34. CREATE FARM SUBNETWORKS — FIGURES F–I
# ==============================================================================

farm_network_summaries <- map_dfr(

  FARM_ORDER,

  function(farm) {

    result <-
      make_farm_subnetwork(
        farm
      )


    if (
      is.null(
        result
      )
    ) {

      return(

        tibble(

          Farm =
            farm,

          Samples_required =
            FARM_MIN_SAMPLES,

          Nodes =
            NA_real_,

          Edges =
            NA_real_,

          Density =
            NA_real_,

          Mean_degree =
            NA_real_
        )
      )
    }


    result
  }
)


write.csv(

  farm_network_summaries,

  file.path(
    TABLE_DIR,
    "farm_subnetwork_summary.csv"
  ),

  row.names =
    FALSE
)


# ==============================================================================
# 35. SAMPLE REPRESENTATION CHECK
# ==============================================================================

sample_check <- tibble(

  Sample =
    as.character(
      sample_metadata$Sample
    ),

  Farm =
    as.character(
      sample_metadata$Farm
    ),

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

  row.names =
    FALSE
)


cat(
  "\n========================================\n"
)

cat(
  "SAMPLE REPRESENTATION CHECK\n"
)

cat(
  "========================================\n"
)


print(
  sample_check,
  n = Inf
)


# ==============================================================================
# 36. FIGURE MANIFEST
#
# Keeps figure labels and filenames standardised for GitHub and thesis writing.
# ==============================================================================

figure_manifest <- tibble(

  Figure =
    LETTERS[
      1:9
    ],

  Title =
    unname(
      FIGURE_TITLES[
        LETTERS[
          1:9
        ]
      ]
    ),

  PNG_file =
    c(

      "Fig_A_phylum_relative_abundance.png",

      "Fig_B_genus_relative_abundance.png",

      "Fig_C_top20_genera_heatmap.png",

      "Fig_D_global_association_network.png",

      "Fig_E_top20_global_network_genera.png",

      "Fig_F_Clanwilliam_association_subnetwork.png",

      "Fig_G_Dendron_association_subnetwork.png",

      "Fig_H_Mamusha_association_subnetwork.png",

      "Fig_I_Wesselesbron_association_subnetwork.png"
    ),

  Description =
    c(

      "Relative abundance of the most abundant bacterial phyla across samples grouped by farm.",

      "Relative abundance of the most abundant bacterial genera across samples grouped by farm.",

      "Heatmap of the 20 most abundant identified bacterial genera across samples.",

      "Global genus-level microbial association network inferred using all samples.",

      "Twenty genera with the highest degree in the global association network.",

      "Clanwilliam subnetwork containing globally supported edges among genera detected at Clanwilliam.",

      "Dendron subnetwork containing globally supported edges among genera detected at Dendron.",

      "Mamusha subnetwork containing globally supported edges among genera detected at Mamusha.",

      "Wesselesbron subnetwork containing globally supported edges among genera detected at Wesselesbron."
    )
)


write.csv(

  figure_manifest,

  file.path(
    OUTPUT_DIR,
    "figure_manifest.csv"
  ),

  row.names =
    FALSE
)


# ==============================================================================
# 37. ANALYSIS PARAMETERS
# ==============================================================================

analysis_parameters <- tibble(

  Parameter =
    c(
      "SILVA release",
      "Top phyla",
      "Top genera",
      "Heatmap genera",
      "Global prevalence threshold",
      "Correlation method",
      "Correlation threshold",
      "FDR threshold",
      "CLR pseudocount",
      "Minimum farm samples for subnetwork inclusion"
    ),

  Value =
    c(
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

  row.names =
    FALSE
)


# ==============================================================================
# 38. SAVE R SESSION INFORMATION
# ==============================================================================

capture.output(

  sessionInfo(),

  file =
    file.path(
      OUTPUT_DIR,
      "R_sessionInfo.txt"
    )
)


# ==============================================================================
# 39. FINISHED
# ==============================================================================

cat(
  "\n========================================\n"
)

cat(
  "ANALYSIS COMPLETE\n"
)

cat(
  "========================================\n\n"
)

cat(
  "Figures saved to:\n",
  FIGURE_DIR,
  "\n\n"
)

cat(
  "Tables saved to:\n",
  TABLE_DIR,
  "\n\n"
)

cat(
  "Networks saved to:\n",
  NETWORK_DIR,
  "\n\n"
)

cat(
  "Standardised figures:\n"
)

cat(
  "A. Phylum-level relative abundance\n"
)

cat(
  "B. Genus-level relative abundance\n"
)

cat(
  "C. Top 20 genera heatmap\n"
)

cat(
  "D. Global genus-level association network\n"
)

cat(
  "E. Most connected genera in the global network\n"
)

cat(
  "F. Clanwilliam association subnetwork\n"
)

cat(
  "G. Dendron association subnetwork\n"
)

cat(
  "H. Mamusha association subnetwork\n"
)

cat(
  "I. Wesselesbron association subnetwork\n"
)
