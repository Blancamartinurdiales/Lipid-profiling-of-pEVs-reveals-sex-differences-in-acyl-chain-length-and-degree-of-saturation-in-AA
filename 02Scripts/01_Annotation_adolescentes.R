# =========================================================
# Lipid annotation workflow for lipidr compatibility
# Working directory: code
# =========================================================

setwd("C:/Users/user/Desktop/P3_Adolescentes/Zenodo_Adolescentes")

rm(list = ls())

library(lipidr)

# =========================================================
# Load Adolescents datasets (positive and negative modes)
# =========================================================
datos_neg <- read.table(
  "00Data/humano_adolescentes/DATOS_NORMALIZADOS_HUMANO_NEGATIVO_RESUMIDOS_SIN_OUTLIERS.txt",
  header = T,      
  sep = ","         
)
datos_neg$ET9_H <- NULL

datos_pos <- read.table(
  "00Data/humano_adolescentes/DATOS_NORMALIZADOS_HUMANO_POSITIVO_sinET9_RESUMIDOS_SIN_OUTLIERS.txt",
  header = TRUE,  
  sep = ","         
)

# Combine positive and negative mode datasets
datos_adolescentes <- rbind(datos_neg,datos_pos)

# Remove unnecessary columns (QC / artifacts)
columna_comun <- datos_adolescentes[, "X", drop = FALSE]
# =========================================================
# Split data by sex (H = men, M = women)
# =========================================================

# Male samples
datos_adolescentes_H <- cbind(
  columna_comun,
  datos_adolescentes[, grepl("H$", colnames(datos_adolescentes))]
)
# Female samples
datos_adolescentes_M <- cbind(
  columna_comun,
  datos_adolescentes[, grepl("M$", colnames(datos_adolescentes))]
)

# =========================================================
# Lipid name preprocessing for lipidr compatibility
# (duplicate handling + naming standardization)
# =========================================================

# -------------------------
# MEN (H)
# -------------------------

# Remove RT annotation if present
datos_adolescentes_H$X <- gsub("\\s*RT:\\d+\\.\\d+", "", datos_adolescentes_H$X)

# Identify duplicated lipids
datos_adolescentes_H_rep <- table(datos_adolescentes_H$X)
datos_adolescentes_H_rep[datos_adolescentes_H_rep > 1]

lipidos_adolescentes_H<-datos_adolescentes_H$X
dup_adolescentes_H <- lipidos_adolescentes_H[(duplicated(lipidos_adolescentes_H))]

# Add replicate index to duplicated lipids
conteo <- ave(lipidos_adolescentes_H, lipidos_adolescentes_H, FUN=seq_along)
total <- ave(lipidos_adolescentes_H, lipidos_adolescentes_H, FUN=length)
lipidos_final_adolescentes_H <- ifelse(
  total>1,
  paste0(lipidos_adolescentes_H," (",conteo,")"),
  lipidos_adolescentes_H)
lipidos_final_adolescentes_H

# -------------------------
# WOMEN (M)
# -------------------------

# Remove RT annotation if present
datos_adolescentes_M$X <- gsub("\\s*RT:\\d+\\.\\d+", "", datos_adolescentes_H$X)

# Identify duplicated lipids
datos_adolescentes_M_rep <- table(datos_adolescentes_M$X)
datos_adolescentes_M_rep[datos_adolescentes_M_rep > 1]

lipidos_adolescentes_M<-datos_adolescentes_M$X
dup_adolescentes_M <- lipidos_adolescentes_M[(duplicated(lipidos_adolescentes_M))]

# Add replicate index to duplicated lipids
conteo <- ave(lipidos_adolescentes_M, lipidos_adolescentes_M, FUN=seq_along)
total <- ave(lipidos_adolescentes_M, lipidos_adolescentes_M, FUN=length)
lipidos_final_adolescentes_M <- ifelse(
  total>1,
  paste0(lipidos_adolescentes_M," (",conteo,")"),
  lipidos_adolescentes_M)
lipidos_final_adolescentes_M

# =========================================================
# Validate lipid nomenclature with lipidr
# =========================================================

annotate_lipids(lipidos_final_adolescentes_H, no_match = "warn")

annotate_lipids(lipidos_final_adolescentes_M, no_match = "warn")

# =========================================================
# String cleaning for lipidr standard format (MEN)
# =========================================================

lipidos_final_adolescentes_H_1 <- gsub("(\\d+:\\d)_(\\d+:\\d)", "\\1/\\2", lipidos_final_adolescentes_H)
lipidos_final_adolescentes_H_2 <- gsub("(\\d+:\\d)_", "\\1/", lipidos_final_adolescentes_H_1)
lipidos_final_adolescentes_H_3 <-  gsub("(\\d+:\\d[a-zA-Z]*)_", "\\1/", lipidos_final_adolescentes_H_2)
lipidos_final_adolescentes_H_4 <- gsub("\\bd(\\d+:\\d[a-zA-Z]*)", "\\1", lipidos_final_adolescentes_H_3)
lipidos_final_adolescentes_H_5 <- gsub("^Cer_([A-Z]+)", "Cer\\1", lipidos_final_adolescentes_H_4)
lipidos_final_adolescentes_H_6 <- gsub("(\\d)e([ /])","\\1\\2",lipidos_final_adolescentes_H_5)
lipidos_final_adolescentes_H_7 <- gsub("e([ (])", "\\1", lipidos_final_adolescentes_H_6)
lipidos_final_adolescentes_H_8 <- gsub("^HexCer_([A-Z]+)", "HexCer\\1", lipidos_final_adolescentes_H_7)
lipidos_final_adolescentes_H_9 <- gsub(" t(\\d+)", " \\1", lipidos_final_adolescentes_H_8)
lipidos_final_adolescentes_H_10 <- gsub("\\+\\w+", "", lipidos_final_adolescentes_H_9)
lipidos_final_adolescentes_H_11 <- gsub("e$", "", lipidos_final_adolescentes_H_10)

annotate_lipids(lipidos_final_adolescentes_H_11, no_match = "warn")

# =========================================================
# String cleaning for lipidr standard format (WOMEN)
# =========================================================

lipidos_final_adolescentes_M_1 <- gsub("(\\d+:\\d)_(\\d+:\\d)", "\\1/\\2", lipidos_final_adolescentes_M)
lipidos_final_adolescentes_M_2 <- gsub("(\\d+:\\d)_", "\\1/", lipidos_final_adolescentes_M_1)
lipidos_final_adolescentes_M_3 <-  gsub("(\\d+:\\d[a-zA-Z]*)_", "\\1/", lipidos_final_adolescentes_M_2)
lipidos_final_adolescentes_M_4 <- gsub("\\bd(\\d+:\\d[a-zA-Z]*)", "\\1", lipidos_final_adolescentes_M_3)
lipidos_final_adolescentes_M_5 <- gsub("^Cer_([A-Z]+)", "Cer\\1", lipidos_final_adolescentes_M_4)
lipidos_final_adolescentes_M_6 <- gsub("(\\d)e([ /])","\\1\\2",lipidos_final_adolescentes_M_5)
lipidos_final_adolescentes_M_7 <- gsub("e([ (])", "\\1", lipidos_final_adolescentes_M_6)
lipidos_final_adolescentes_M_8 <- gsub("^HexCer_([A-Z]+)", "HexCer\\1", lipidos_final_adolescentes_M_7)
lipidos_final_adolescentes_M_9 <- gsub(" t(\\d+)", " \\1", lipidos_final_adolescentes_M_8)
lipidos_final_adolescentes_M_10 <- gsub("\\+\\w+", "", lipidos_final_adolescentes_M_9)
lipidos_final_adolescentes_M_11 <- gsub("e$", "", lipidos_final_adolescentes_M_10)

annotate_lipids(lipidos_final_adolescentes_M_11, no_match = "warn")

# =========================================================
# Assign cleaned lipid names back to dataset
# =========================================================

datos_adolescentes_H$X <- lipidos_final_adolescentes_H_11
datos_adolescentes_M$X <- lipidos_final_adolescentes_M_11
# Merge male and female datasets
datos_adolescentes <- cbind(datos_adolescentes_H,datos_adolescentes_M)

# =========================================================
# Export cleaned datasets
# =========================================================

getwd()
write.csv(datos_adolescentes_H, "00Data/datos_adolescentes_wo_H.csv", row.names = F )
write.csv(datos_adolescentes_M, "00Data/datos_adolescentes_wo_M.csv", row.names = F )
write.csv(datos_adolescentes, "00Data/datos_adolescentes_wo.csv", row.names = F )
