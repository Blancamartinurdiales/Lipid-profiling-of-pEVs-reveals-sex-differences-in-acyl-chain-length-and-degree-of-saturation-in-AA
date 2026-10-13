# =========================================================
# PCA analysis of lipidomic profiles - Adolescent cohort
#
# Purpose:
# Perform principal component analysis (PCA) on lipidomic
# abundance data after KNN imputation of missing values.
# PCA is used to evaluate sample clustering according to
# experimental group and sex.
#
# Input files:
#   - datos_adolescentes_wo_H.csv
#   - datos_adolescentes_wo_M.csv
#   - metadata_AA_H.xlsx
#   - metadata_AA_M.xlsx
#   - metadata_AA.xlsx
#   - function_pcaGenes_2.r
#
# Output:
#   - PCA score plots (PC1 vs PC2 and PC2 vs PC3)
#
# =========================================================

# Working directory
setwd("C:/Users/user/Desktop/P3_Adolescentes/Zenodo_Adolescentes")
rm(list = ls())

# =========================================================
# Load required libraries
# =========================================================

library(lipidr)
library(readxl)
library(gridExtra)
library(dplyr)
library(ggplot2)
library(DT)
library(ggrepel)
library(ggpubr)

# =========================================================
# Load custom PCA function
#
# This function performs PCA and returns scores and
# explained variance for downstream visualization.
# =========================================================

source("01Functions/function_pcaGenes_2.r")

# =========================================================
# Load lipidomics datasets and sample metadata
#
# Male and female datasets are loaded separately and merged
# into a single abundance matrix.
# =========================================================

AA_H <- read.csv("00Data/datos_adolescentes_wo_H.csv")
AA_M <- read.csv("00Data/datos_adolescentes_wo_M.csv")

metadata_AA_H <- read_excel("00Data/metadata_AA_H.xlsx")
metadata_AA_M <- read_excel("00Data/metadata_AA_M.xlsx")
metadata_AA <- read_excel("00Data/metadata_AA.xlsx")

# =========================================================
# Merge lipid abundance matrices
#
# Lipid identifiers are stored in column "X" and assigned
# as row names of the abundance matrix.
# =========================================================

datos <- full_join(AA_H, AA_M, by = "X")

datos_mat <- datos[, -1]
rownames(datos_mat) <- datos$X

# =========================================================
# Missing value imputation using K-nearest neighbors (KNN)
#
# Missing lipid abundance values are estimated using the
# impute.knn() function. This step is required before PCA.
# =========================================================

if (!require("impute")) {
  install.packages("impute")
  library(impute)
}

# Remove existing random seed to ensure reproducibility
# according to the impute package defaults
if (exists(".Random.seed")) rm(.Random.seed)

datos_imputed <- impute.knn(as.matrix(datos_mat))

# Extract imputed abundance matrix
datos_mat <- datos_imputed$data

# Store imputation seed information
datos_imputed$rng.seed
datos_imputed$rng.state

# =========================================================
# Prepare sample metadata
#
# Create combined labels including experimental group and
# sex, and ensure metadata matches sample order in the
# abundance matrix.
# =========================================================

metadata_AA$grupo <- paste(
  metadata_AA$group,
  metadata_AA$sex,
  sep = "_"
)

metadata_AA[
  match(colnames(datos_mat), metadata_AA$sample),
]

# =========================================================
# Define color palette for sample groups
# =========================================================

palette4 <- c(
  "AA_H" = "#dd5129",
  "AA_M" = "#f29f05",
  "CONTROL_H" = "#0f7ba2",
  "CONTROL_M" = "#43b5a0"
)

# =========================================================
# Principal Component Analysis (PCA)
#
# PCA is performed on the imputed lipid abundance matrix
# to identify major sources of variation across samples.
# =========================================================

mi.pca <- pcaGenes(datos_mat)

# Extract PCA scores for the first three principal
# components and calculate explained variance.
mi.pca.df <- as.data.frame(mi.pca$scores[, 1:3])

mi.pca.df$grupo <- metadata_AA$grupo
mi.pca.df$var.exp <- round(mi.pca$var.exp * 100, 2)

rownames(mi.pca.df) <- metadata_AA$sample
mi.pca.df$etiqueta <- metadata_AA$sample

# =========================================================
# PCA visualization: PC1 vs PC2
#
# Samples are colored according to group and displayed
# together with 95% confidence ellipses.
# =========================================================

PC1_PC2_AA <- ggplot(
  mi.pca.df,
  aes(
    x = V1,
    y = V2,
    colour = grupo
  )
) +
  
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  
  stat_ellipse(
    geom = "polygon",
    aes(fill = grupo),
    alpha = 0.2,
    level = 0.95,
    show.legend = FALSE
  ) +
  
  geom_point(
    shape = 20,
    alpha = 0.9,
    size = 2
  ) +
  
  scale_color_manual(values = palette4) +
  scale_fill_manual(values = palette4) +
  
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = "transparent")
  ) +
  
  xlab(
    paste0(
      "PC1: ",
      mi.pca.df$var.exp[1],
      "% explained variance"
    )
  ) +
  
  ylab(
    paste0(
      "PC2: ",
      mi.pca.df$var.exp[2],
      "% explained variance"
    )
  )

PC1_PC2_AA

# =========================================================
# PCA visualization: PC2 vs PC3
#
# Visualization of the second and third principal
# components using the same group annotations.
# =========================================================

PC2_PC3_AA <- ggplot(
  mi.pca.df,
  aes(
    x = V2,
    y = V3,
    colour = grupo
  )
) +
  
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  
  stat_ellipse(
    geom = "polygon",
    aes(fill = grupo),
    alpha = 0.2,
    level = 0.95,
    show.legend = FALSE
  ) +
  
  geom_point(
    shape = 20,
    alpha = 0.9,
    size = 2
  ) +
  
  scale_color_manual(values = palette4) +
  scale_fill_manual(values = palette4) +
  
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = "transparent")
  ) +
  
  xlab(
    paste0(
      "PC2: ",
      mi.pca.df$var.exp[2],
      "% explained variance"
    )
  ) +
  
  ylab(
    paste0(
      "PC3: ",
      mi.pca.df$var.exp[3],
      "% explained variance"
    )
  )

PC2_PC3_AA


