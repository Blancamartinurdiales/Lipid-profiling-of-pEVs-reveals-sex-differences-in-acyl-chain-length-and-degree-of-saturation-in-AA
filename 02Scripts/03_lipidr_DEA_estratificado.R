# =========================================================
# Differential lipid analysis and volcano plot visualization
# using the lipidr package
#
# Purpose:
# Perform differential lipid abundance analysis between
# alcohol-exposed and control adolescents separately in
# males and females. Generate volcano plots colored by
# lipid class and lipid superclass.
#
# Input files:
#   - datos_adolescentes_wo_H.csv
#   - datos_adolescentes_wo_M.csv
#   - metadata_AA_H.xlsx
#   - metadata_AA_M.xlsx
#
# Output:
#   - Differential expression result objects (.RData)
#   - Interactive result tables
#   - Volcano plots by lipid class
#   - Volcano plots by lipid superclass
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

# =========================================================
# Load lipid abundance data and sample metadata
#
# Male and female datasets are analyzed separately in order
# to investigate sex-specific lipidomic responses.
# =========================================================

AA_H <- read.csv("00Data/datos_adolescentes_wo_H.csv")
AA_M <- read.csv("00Data/datos_adolescentes_wo_M.csv")

metadata_AA_H <- read_excel("00Data/metadata_AA_H.xlsx")
metadata_AA_M <- read_excel("00Data/metadata_AA_M.xlsx")

# =========================================================
# Create LipidomicsExperiment objects
#
# Lipidr requires abundance matrices and metadata to be
# stored in a dedicated LipidomicsExperiment object.
# =========================================================

d_AA_H <- as_lipidomics_experiment(
  AA_H,
  logged = TRUE,
  normalized = TRUE
)

d_AA_H <- add_sample_annotation(
  d_AA_H,
  metadata_AA_H
)

assay(d_AA_H)

# =========================================================
# Differential lipid analysis: males
#
# Comparison:
#   Alcohol-exposed adolescents vs controls
# =========================================================

et_c_AA_H <- de_analysis(
  d_AA_H,
  AA - CONTROL
)

class(et_c_AA_H)

# =========================================================
# Differential lipid analysis: females
#
# Comparison:
#   Alcohol-exposed adolescents vs controls
# =========================================================

d_AA_M <- as_lipidomics_experiment(
  AA_M,
  logged = TRUE,
  normalized = TRUE
)

d_AA_M <- add_sample_annotation(
  d_AA_M,
  metadata_AA_M
)

assay(d_AA_M)

et_c_AA_M <- de_analysis(
  d_AA_M,
  AA - CONTROL
)

class(et_c_AA_M)

# =========================================================
# Generate standard volcano plots
#
# Volcano plots display differential abundance according
# to fold change and adjusted statistical significance.
# =========================================================

volcano_H <- plot_results_volcano(et_c_AA_H)
volcano_M <- plot_results_volcano(et_c_AA_M)

# =========================================================
# Save differential analysis results
# =========================================================

save(
  et_c_AA_H,
  file = "04Results/IM.RData"
)

save(
  et_c_AA_M,
  file = "04Results/IF.RData"
)

# =========================================================
# Interactive result tables
#
# Results can be explored through searchable and sortable
# HTML tables.
# =========================================================


DT::datatable(et_c_AA_H,
              caption = c("Tabla interactiva del resultado análisis AUD vs CONTROL_AA_H"),
              escape = FALSE, 
              options = list(
                pageLength = 100,
                scrollX = TRUE
              ),
              filter = 'top', #Pone una opciC3n para filtrar encima de la tabla
              colnames = c("Lipids", "Class", "Numero de Carbonos", "Insaturaciones", "istd", "MoleculeID", "Contrast", "LogFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")) %>%
  formatRound(columns = c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"), digits = 3) # Redondear a 3

DT::datatable(et_c_AA_M,
              caption = c("Tabla interactiva del resultado análisis AUD vs CONTROL_AA_M"),
              escape = FALSE, 
              options = list(
                pageLength = 100,
                scrollX = TRUE
              ),
              filter = 'top', #Pone una opciC3n para filtrar encima de la tabla
              colnames = c("Lipids", "Class", "Numero de Carbonos", "Insaturaciones", "istd", "MoleculeID", "Contrast", "LogFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")) %>%
  formatRound(columns = c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"), digits = 3) # Redondear a 3




# =========================================================
# Define lipid superclasses
#
# Individual lipid classes are grouped into broader lipid
# categories to facilitate biological interpretation.
# =========================================================

# =========================================================
# Volcano plot by lipid superclass - males
# =========================================================
et_c_AA_H <- et_c_AA_H %>%
  
  mutate(
    Superclass = case_when(
      
      # Fatty Acyls
      Class %in% c("FA", "FAHFA", "ACar") ~ "Fatty Acyls",
      
      # Glycerolipids
      Class %in% c("TG", "DG", "MGDG","AcylGlcADG") ~ "Glycerolipids",
      
      # Glycerophospholipids
      Class %in% c(
        "PC", "PE", "PI", "PA", "LPC", "LPE",
        "EtherPC", "EtherPE", "EtherOxPC", "OxPC"
      ) ~ "Glycerophospholipids",
      
      # Sphingolipids
      Class %in% c(
        "SM",
        "CerNS", "CerNDS", "CerADS", "CerAS",
        "CerAP", "CerNP",
        "HexCerNS", "HexCerNDS",
        "SHexCer"
      ) ~ "Sphingolipids",
      
      # Sterol Lipids
      Class %in% c("CE") ~ "Sterol Lipids",
      
      # Otros
      TRUE ~ "Others"
    )
  )

table(et_c_AA_H$Superclass)

# =========================================================
# Color palette for lipid superclasses
# =========================================================

superclass_palette <- c(
  "Fatty Acyls" = "#D8443CFF",
  "Glycerolipids" = "#F4C40FFF",
  "Glycerophospholipids" = "#92C051FF",
  "Sphingolipids" = "#633372FF",
  "Sterol Lipids" = "#E6A2A6FF",
  "Others" = "#FE9B00FF"
)


volcano_superclass_H <- ggplot(
  et_c_AA_H,
  aes(
    x = logFC,
    y = -log10(adj.P.Val),
    color = Superclass
  )
) +
  
  geom_point(size = 1.5, alpha = 0.8) +
  
  scale_color_manual(values = superclass_palette) +
  
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  
  theme_minimal() +
  
  labs(
    title = "IM comparation",
    x = "log2 Fold Change",
    y = "-log10 adjusted p-value",
    color = "Superclass"
  ) +
  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

volcano_superclass_H

# =========================================================
# Volcano plot by lipid superclass - females
# =========================================================
et_c_AA_M <- et_c_AA_M %>%
  
  mutate(
    Superclass = case_when(
      
      # Fatty Acyls
      Class %in% c("FA", "FAHFA", "ACar") ~ "Fatty Acyls",
      
      # Glycerolipids
      Class %in% c("TG", "DG", "MGDG","AcylGlcADG") ~ "Glycerolipids",
      
      # Glycerophospholipids
      Class %in% c(
        "PC", "PE", "PI", "PA", "LPC", "LPE",
        "EtherPC", "EtherPE", "EtherOxPC", "OxPC"
      ) ~ "Glycerophospholipids",
      
      # Sphingolipids
      Class %in% c(
        "SM",
        "CerNS", "CerNDS", "CerADS", "CerAS",
        "CerAP", "CerNP",
        "HexCerNS", "HexCerNDS",
        "SHexCer"
      ) ~ "Sphingolipids",
      
      # Sterol Lipids
      Class %in% c("CE") ~ "Sterol Lipids",
      
      # Otros
      TRUE ~ "Others"
    )
  )



table(et_c_AA_M$Superclass)


superclass_palette <- c(
  "Fatty Acyls" = "#D8443CFF",
  "Glycerolipids" = "#F4C40FFF",
  "Glycerophospholipids" = "#92C051FF",
  "Sphingolipids" = "#633372FF",
  "Sterol Lipids" = "#E6A2A6FF",
  "Others" = "#FE9B00FF"
)

volcano_superclass_M <- ggplot(
  et_c_AA_M,
  aes(
    x = logFC,
    y = -log10(adj.P.Val),
    color = Superclass
  )
) +
  
  geom_point(size = 1.5, alpha = 0.8) +
  
  scale_color_manual(values = superclass_palette) +
  
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  
  theme_minimal() +
  
  labs(
    title = "IF comparation",
    x = "log2 Fold Change",
    y = "-log10 adjusted p-value",
    color = "Superclass"
  ) +
  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

volcano_superclass_M

