# ============================================================
# Configuration file for 16S rRNA amplicon analysis workflow
# ============================================================

# ----------------------------
# Input data paths
# ----------------------------

# Base directory containing QIIME2 outputs
qiime_dir <- "data/qiime2"

# QIIME2 artifacts
feature_table_qza <- file.path(qiime_dir, "table.qza")
tree_qza          <- file.path(qiime_dir, "rooted-tree.qza")
taxonomy_qza     <- file.path(qiime_dir, "taxonomy.qza")

# Sample metadata
metadata_tsv <- "data/metadata.tsv"

# ----------------------------
# Output settings
# ----------------------------

# Directory where plots and tables will be written
output_dir <- "output"

# Create output directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# ----------------------------
# Taxonomic filtering settings
# ----------------------------

# Kingdom to retain
target_kingdom <- "Bacteria"

# Remove common non-target taxa
remove_mitochondria <- TRUE
remove_chloroplast  <- TRUE

mitochondria_label <- "Mitochondria"
chloroplast_label  <- "Chloroplast"


# ----------------------------
# Controls (neg for decontam, and mocks)
# ----------------------------

# Column in metadata identifying control samples
control_column <- "sampletype"

# Values in that column that indicate negative controls
negative_control_labels <- c(
  "Negative control",
  "Blank",
  "Water",
  "PBS"
)

# Values corresponding to mock communities
mock_labels <- c(
  "Positive control",
  "Mock",
  "LogMC"
)

# If using DNA concentration for frequency-based decontam
# set column name here; otherwise leave as NULL
dna_concentration_column <- NULL

# Expected relative abundance for standard mock
mock_expected_path <- "data/mock_expected_abundance.csv"

# Expected relative abundance for logarithmic mock (optional)
log_mock_expected_path <- "data/log_mock_expected_abundance.csv"

# Name of logarithmic mock sample (if present)
log_mock_sample <- "LogMC"

# ----------------------------
# Diversity analysis settings
# ----------------------------

# Distance metric for beta diversity
beta_distance <- "wunifrac"

# Metadata variable to test - suggestion here is to explore environmental samples as recommended in the readme
met_var <- "sampletype"

# Define your main sample-type of interest, if you have included environmental controls

main_group <- "semen"

# ----------------------------
# Plotting defaults
# ----------------------------

plot_width  <- 7
plot_height <- 5
plot_dpi    <- 300
