# =========================================================
# Lipid chain length and unsaturation trend analysis
#
# Purpose:
# Investigate trends in differential lipid abundance
# (log2 fold change) according to lipid chain length and
# degree of unsaturation across major lipid classes.
#
# Trend plots are generated separately for males and
# females using differential analysis results obtained
# from the lipidr workflow.
#
# Lipid classes analyzed:
#   - TG  (Triacylglycerols)
#   - FA  (Fatty Acyls)
#   - SM  (Sphingomyelins)
#   - Cer (Ceramides)
#   - PC  (Phosphatidylcholines)
#   - LPC (Lysophosphatidylcholines)
#   - PE  (Phosphatidylethanolamines)
#
# Input files:
#   - IF.RData
#   - IM.RData
#
# Output:
#   - Trend plots of logFC vs chain length
#   - Trend plots of logFC vs unsaturation degree
#
# =========================================================

# Working directory
setwd("C:/Users/user/Desktop/P3_Adolescentes")
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
# Load differential lipid analysis results
#
# IM = male subgroup
# IF = female subgroup
# =========================================================

load("04Results/IF.RData")
load("04Results/IM.RData")

# =========================================================
# Filter lipid classes of interest
#
# Internal standards (istd) are excluded from all analyses.
# =========================================================

# ---------------------------------------------------------
# Triacylglycerols (TG)
# ---------------------------------------------------------

class_tg_AA_M <- et_c_AA_M %>%
  filter(Class == "TG", !istd)

class_tg_AA_H <- et_c_AA_H %>%
  filter(Class == "TG", !istd)

# ---------------------------------------------------------
# Fatty Acyls (FA)
# ---------------------------------------------------------

class_fa_AA_M <- et_c_AA_M %>%
  filter(Class == "FA", !istd)

class_fa_AA_H <- et_c_AA_H %>%
  filter(Class == "FA", !istd)

# ---------------------------------------------------------
# Lysophosphatidylcholines (LPC)
# ---------------------------------------------------------

class_lpc_AA_M <- et_c_AA_M %>%
  filter(Class == "LPC", !istd)

class_lpc_AA_H <- et_c_AA_H %>%
  filter(Class == "LPC", !istd)

# ---------------------------------------------------------
# Sphingomyelins (SM)
# ---------------------------------------------------------

class_sm_AA_M <- et_c_AA_M %>%
  filter(Class == "SM", !istd)

class_sm_AA_H <- et_c_AA_H %>%
  filter(Class == "SM", !istd)

# ---------------------------------------------------------
# Phosphatidylcholines (PC)
# ---------------------------------------------------------

class_pc_AA_M <- et_c_AA_M %>%
  filter(Class == "PC", !istd)

class_pc_AA_H <- et_c_AA_H %>%
  filter(Class == "PC", !istd)

# ---------------------------------------------------------
# Phosphatidylethanolamines (PE)
# ---------------------------------------------------------

class_pe_AA_M <- et_c_AA_M %>%
  filter(Class == "PE", !istd)

class_pe_AA_H <- et_c_AA_H %>%
  filter(Class == "PE", !istd)

# ---------------------------------------------------------
# Ceramides (CerNS and CerNDS)
# ---------------------------------------------------------

class_cer_AA_M <- et_c_AA_M %>%
  filter(Class %in% c("CerNDS", "CerNS"), !istd)

class_cer_AA_H <- et_c_AA_H %>%
  filter(Class %in% c("CerNDS", "CerNS"), !istd)

# =========================================================
# Custom trend plot function
#
# This modified version of plot_trend() allows user-defined
# y-axis limits while preserving all observations.
#
# coord_cartesian() is used to zoom without removing data.
# =========================================================

plot_trend_mod <- function(
    de_results,
    annotation = c("length", "unsat"),
    y_limits = NULL
) {
  
  annotation <- match.arg(annotation)
  
  x <- rlang::sym(
    ifelse(
      annotation == "length",
      "total_cl",
      "total_cs"
    )
  )
  
  x_lab <- ifelse(
    annotation == "length",
    "Total chain length",
    "Total unsaturated bonds"
  )
  
  x_breaks <- if (annotation == "length") {
    seq(10, 80, 2)
  } else {
    seq(0, 12, 1)
  }
  
  p <- ggplot(
    de_results,
    aes(
      !!x,
      logFC,
      color = contrast
    )
  ) +
    geom_hline(
      color = "black",
      yintercept = 0,
      lty = 2
    ) +
    geom_smooth(alpha = 0.2) +
    labs(
      x = x_lab,
      y = "logFC"
    ) +
    scale_x_continuous(
      breaks = x_breaks
    ) +
    theme(
      legend.position = "bottom"
    )
  
  # Apply custom y-axis limits if specified
  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(
      ylim = y_limits
    )
  }
  
  return(p)
}

# =========================================================
# Trend plots using default axis limits
#
# Visualization of logFC trends according to:
#   - Total chain length
#   - Number of unsaturated bonds
#
# Generated separately for males and females.
# =========================================================

# ---------------------------------------------------------
# TG
# ---------------------------------------------------------

plot_trend(class_tg_AA_H, "length")
plot_trend(class_tg_AA_M, "length")

plot_trend(class_tg_AA_H, "unsat")
plot_trend(class_tg_AA_M, "unsat")

# ---------------------------------------------------------
# FA
# ---------------------------------------------------------

plot_trend(class_fa_AA_H, "length")
plot_trend(class_fa_AA_M, "length")

plot_trend(class_fa_AA_H, "unsat")
plot_trend(class_fa_AA_M, "unsat")

# ---------------------------------------------------------
# SM
# ---------------------------------------------------------

plot_trend(class_sm_AA_H, "length")
plot_trend(class_sm_AA_M, "length")

plot_trend(class_sm_AA_H, "unsat")
plot_trend(class_sm_AA_M, "unsat")

# ---------------------------------------------------------
# Ceramides
# ---------------------------------------------------------

plot_trend(class_cer_AA_H, "length")
plot_trend(class_cer_AA_M, "length")

plot_trend(class_cer_AA_H, "unsat")
plot_trend(class_cer_AA_M, "unsat")

# ---------------------------------------------------------
# PC
# ---------------------------------------------------------

plot_trend(class_pc_AA_H, "length")
plot_trend(class_pc_AA_M, "length")

plot_trend(class_pc_AA_H, "unsat")
plot_trend(class_pc_AA_M, "unsat")

# ---------------------------------------------------------
# LPC
# ---------------------------------------------------------

plot_trend(class_lpc_AA_H, "length")
plot_trend(class_lpc_AA_M, "length")

plot_trend(class_lpc_AA_H, "unsat")
plot_trend(class_lpc_AA_M, "unsat")

# ---------------------------------------------------------
# PE
# ---------------------------------------------------------

plot_trend(class_pe_AA_H, "length")
plot_trend(class_pe_AA_M, "length")

plot_trend(class_pe_AA_H, "unsat")
plot_trend(class_pe_AA_M, "unsat")

# =========================================================
# Trend plots using customized y-axis limits
#
# Axis limits are adjusted individually for each lipid
# class to improve visualization of trend patterns.
# =========================================================

# ---------------------------------------------------------
# TG
# ---------------------------------------------------------

plot_trend_mod(
  class_tg_AA_M,
  y_limits = c(-10, 10),
  "length"
)

plot_trend_mod(
  class_tg_AA_H,
  y_limits = c(-10, 10),
  "length"
)

plot_trend_mod(
  class_tg_AA_M,
  y_limits = c(-10, 10),
  "unsat"
)

plot_trend_mod(
  class_tg_AA_H,
  y_limits = c(-10, 10),
  "unsat"
)

# ---------------------------------------------------------
# FA
# ---------------------------------------------------------

plot_trend_mod(
  class_fa_AA_M,
  y_limits = c(-5, 8),
  "length"
)

plot_trend_mod(
  class_fa_AA_H,
  y_limits = c(-5, 8),
  "length"
)

plot_trend_mod(
  class_fa_AA_M,
  y_limits = c(-5, 5),
  "unsat"
)

plot_trend_mod(
  class_fa_AA_H,
  y_limits = c(-5, 5),
  "unsat"
)

# ---------------------------------------------------------
# SM
# ---------------------------------------------------------

plot_trend_mod(
  class_sm_AA_M,
  y_limits = c(-5, 10),
  "length"
)

plot_trend_mod(
  class_sm_AA_H,
  y_limits = c(-5, 10),
  "length"
)

plot_trend_mod(
  class_sm_AA_M,
  y_limits = c(-5, 10),
  "unsat"
)

plot_trend_mod(
  class_sm_AA_H,
  y_limits = c(-5, 10),
  "unsat"
)

# ---------------------------------------------------------
# Ceramides
# ---------------------------------------------------------

plot_trend_mod(
  class_cer_AA_M,
  y_limits = c(-10, 10),
  "length"
)

plot_trend_mod(
  class_cer_AA_H,
  y_limits = c(-10, 10),
  "length"
)

plot_trend_mod(
  class_cer_AA_M,
  y_limits = c(-10, 5),
  "unsat"
)

plot_trend_mod(
  class_cer_AA_H,
  y_limits = c(-10, 5),
  "unsat"
)

# ---------------------------------------------------------
# Phosphatidylcholines (PC)
# ---------------------------------------------------------

plot_trend_mod(
  class_pc_AA_M,
  y_limits = c(-4, 4),
  "length"
)

plot_trend_mod(
  class_pc_AA_H,
  y_limits = c(-4, 4),
  "length"
)

plot_trend_mod(
  class_pc_AA_M,
  y_limits = c(-8, 8),
  "unsat"
)

plot_trend_mod(
  class_pc_AA_H,
  y_limits = c(-8, 8),
  "unsat"
)

# ---------------------------------------------------------
# Lysophosphatidylcholines (LPC)
# ---------------------------------------------------------

plot_trend_mod(
  class_lpc_AA_M,
  y_limits = c(-4, 8),
  "length"
)

plot_trend_mod(
  class_lpc_AA_H,
  y_limits = c(-4, 8),
  "length"
)

plot_trend_mod(
  class_lpc_AA_M,
  y_limits = c(-3, 8),
  "unsat"
)

plot_trend_mod(
  class_lpc_AA_H,
  y_limits = c(-3, 8),
  "unsat"
)

# ---------------------------------------------------------
# Phosphatidylethanolamines (PE)
# ---------------------------------------------------------

plot_trend_mod(
  class_pe_AA_M,
  y_limits = c(-10, 10),
  "length"
)

plot_trend_mod(
  class_pe_AA_H,
  y_limits = c(-10, 10),
  "length"
)

plot_trend_mod(
  class_pe_AA_M,
  y_limits = c(-10, 10),
  "unsat"
)

plot_trend_mod(
  class_pe_AA_H,
  y_limits = c(-10, 10),
  "unsat"
)