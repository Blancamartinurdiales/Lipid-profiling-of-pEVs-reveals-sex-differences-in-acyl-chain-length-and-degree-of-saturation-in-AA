
# =========================================================
# Lipidomic heatmap analysis - Adolescent cohort
#
# Purpose:
# Generate a heatmap of normalized lipid abundances from
# adolescent participants, including sample group annotation
# and lipid superclass annotation.
#
# Input files:
#   - datos_adolescentes_wo_H.csv
#   - datos_adolescentes_wo_M.csv
#   - metadata_AA_H.xlsx
#   - metadata_AA_M.xlsx
#
# Output:
#   - Heatmap visualization of scaled lipid abundances
#
# =========================================================

setwd("C:/Users/user/Desktop/P3_Adolescentes/Zenodo_Adolescentes")
rm(list = ls())

library(tidyverse)
library(pheatmap)
library(readxl)
library(lipidr)
library(MetBrewer)

# =========================================================
# Load processed lipid abundance matrices
#
# Data correspond to lipidomics datasets after outlier
# removal. Male and female datasets are loaded separately
# and merged into a single abundance matrix.
# =========================================================

AA_H <- read.csv("00Data/datos_adolescentes_wo_H.csv")
AA_M <- read.csv("00Data/datos_adolescentes_wo_M.csv")

metadata_AA_H <- read_excel("00Data/metadata_AA_H.xlsx")
metadata_AA_M <- read_excel("00Data/metadata_AA_M.xlsx")

# =========================================================
# Merge abundance matrices and define lipid identifiers
#
# Lipid species are stored in column "X", which is used as
# row names for the abundance matrix.
# =========================================================

datos_AA <- full_join(AA_H, AA_M, by = "X")
datos_mat_AA <- datos_AA[, -1]
rownames(datos_mat_AA) <- datos_AA$X

# =========================================================
# Prepare sample metadata
#
# Add sex information, create combined group labels, and
# reorder metadata to match the sample order in the
# abundance matrix.
# =========================================================

metadata_AA_H$sex <- "H"
metadata_AA_M$sex <- "M"

metadata_AA <- bind_rows(metadata_AA_H, metadata_AA_M)

metadata_AA$grupo <- paste(
  metadata_AA$group,
  metadata_AA$sex,
  sep = "_"
)

metadata_AA <- metadata_AA[
  match(colnames(datos_mat_AA), metadata_AA$sample),
]

# =========================================================
# Create sample annotations for heatmap columns
#
# Samples are annotated according to experimental group
# and sex.
# =========================================================

annot_heatmap_col <- data.frame(
  Grupo = metadata_AA$grupo
)

rownames(annot_heatmap_col) <- metadata_AA$sample

ann_colors <- list(
  Grupo = c(
    "AA_H" = "#dd5129",
    "AA_M" = "#f29f05",
    "CONTROL_H" = "#0f7ba2",
    "CONTROL_M" = "#43b5a0"
  )
)

# =========================================================
# Z-score normalization
#
# Each lipid species is scaled independently across
# samples to emphasize relative abundance differences.
# =========================================================

datos_scaled_AA <- t(scale(t(datos_mat_AA)))

# =========================================================
# Assign lipid superclasses
#
# Lipid classes are extracted from lipid names and grouped
# into broader lipid superclasses for row annotation.
# =========================================================

lipidos_AA <- rownames(datos_mat_AA)

# Extract lipid class (text before the first space)
clase_lipidica <- sapply(strsplit(lipidos_AA, " "), `[`, 1)


# Map lipid classes to lipid superclasses
superclase_lipidica <- dplyr::case_when(
  
  clase_lipidica %in% c(
    "FA", "FAHFA", "ACar"
  ) ~ "Fatty Acyls",
  
  clase_lipidica %in% c(
    "TG", "DG", "MGDG", "AcylGlcADG"
  ) ~ "Glycerolipids",
  
  clase_lipidica %in% c(
    "PC", "PE", "PI", "PA", "LPC", "LPE",
    "EtherPC", "EtherPE", "EtherOxPC", "OxPC"
  ) ~ "Glycerophospholipids",
  
  clase_lipidica %in% c(
    "SM",
    "CerNS", "CerNDS", "CerADS", "CerAS",
    "CerAP", "CerNP",
    "HexCerNS", "HexCerNDS",
    "SHexCer"
  ) ~ "Sphingolipids",
  
  clase_lipidica %in% c("CE") ~ "Sterol Lipids",
  
  TRUE ~ "Others"
)

annot_heatmap_row <- data.frame(
  Clase = superclase_lipidica
)

rownames(annot_heatmap_row) <- lipidos_AA

# =========================================================
# Define color palette for lipid superclass annotations
# =========================================================

clases_paleta <- c(
  "#17154FFF", "#2F357CFF", "#6C5D9EFF", "#9D9CD5FF",
  "#B0799AFF", "#F6B3B0FF", "#E48171FF", "#BF3729FF",
  "#E69B00FF", "#F5BB50FF", "#ADA43BFF", "#355828FF"
)

lipid_colors <- setNames(
  rep(
    clases_paleta,
    length.out = length(unique(clase_lipidica))
  ),
  unique(clase_lipidica)
)

# Final color assignment for lipid superclasses
lipid_colors <- c(
  "Fatty Acyls" = "#17154FFF",
  "Glycerolipids" = "#ADA43BFF",
  "Glycerophospholipids" = "#6C5D9EFF",
  "Sphingolipids" = "#E48171FF",
  "Sterol Lipids" = "#E69B00FF",
  "Others" = "#BF3729FF"
)

ann_colors$Clase <- lipid_colors

# =========================================================
# Generate heatmap
#
# Samples are clustered using Euclidean distance and
# Ward's linkage method. Lipids are displayed in their
# original order. Sample groups and lipid superclasses are
# displayed as annotations.
# =========================================================

heatmap_AA <- pheatmap(
  datos_scaled_AA,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D",
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize_row = 6,
  fontsize_col = 8,
  border_col = FALSE,
  annotation_col = annot_heatmap_col,
  annotation_row = annot_heatmap_row,
  annotation_colors = ann_colors,
  angle_col = 45,
  main = "Heatmap - Adolescents intoxicated with alcohol"
)

