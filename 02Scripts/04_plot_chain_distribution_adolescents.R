# =========================================================
# Lipid chain length and unsaturation analysis
# using the lipidr package
#
# Purpose:
# Explore lipid chain-length and unsaturation patterns
# within selected lipid classes following differential
# lipid abundance analysis.
#
# The script generates chain distribution plots for:
#   - Fatty Acyls (FA)
#   - Ceramides (CerNS, CerNDS)
#   - Sphingomyelins (SM)
#   - Lysophosphatidylcholines (LPC)
#   - Phosphatidylcholines (PC)
#   - Phosphatidylethanolamines (PE)
#   - Triacylglycerols (TG)
#
# Chain distributions are visualized using:
#   - Adjusted p-values (adj.P.Val)
#   - Log2 fold changes (logFC)
#
# Input files:
#   - IF.RData
#   - IM.RData
#
# Output:
#   - Chain distribution plots by lipid class
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
# Load differential lipid analysis results
#
# IM = male subgroup
# IF = female subgroup
# =========================================================

load("04Results/IF.RData")
load("04Results/IM.RData")

# =========================================================
# Fatty Acyls (FA)
#
# Evaluate chain-length and unsaturation distributions
# among fatty acyl species.
# =========================================================

class_fa_AA_H <- et_c_AA_H %>%
  filter(Class == "FA")

class_fa_AA_M <- et_c_AA_M %>%
  filter(Class == "FA")

# Distribution according to adjusted p-value
plot_chain_distribution(
  class_fa_AA_H,
  measure = "adj.P.Val"
)

plot_chain_distribution(
  class_fa_AA_M,
  measure = "adj.P.Val"
)

# Distribution according to log2 fold change
plot_chain_distribution(
  class_fa_AA_H,
  measure = "logFC"
)

plot_chain_distribution(
  class_fa_AA_M,
  measure = "logFC"
)

# =========================================================
# Ceramides (CerNS and CerNDS)
#
# Evaluate chain distributions for ceramide species.
# =========================================================

class_cer_AA_H <- et_c_AA_H %>%
  filter(Class %in% c("CerNDS", "CerNS"))

class_cer_AA_M <- et_c_AA_M %>%
  filter(Class %in% c("CerNDS", "CerNS"))

# Distribution according to adjusted p-value
plot_chain_distribution(
  class_cer_AA_H,
  measure = "adj.P.Val"
)

plot_chain_distribution(
  class_cer_AA_M,
  measure = "adj.P.Val"
)

# Distribution according to log2 fold change
plot_chain_distribution(
  class_cer_AA_H,
  measure = "logFC"
)

plot_chain_distribution(
  class_cer_AA_M,
  measure = "logFC"
)

# =========================================================
# Sphingomyelins (SM)
#
# Evaluate chain-length and unsaturation distributions
# among sphingomyelin species.
# =========================================================

class_sm_AA_H <- et_c_AA_H %>%
  filter(Class == "SM")

class_sm_AA_M <- et_c_AA_M %>%
  filter(Class == "SM")

# Distribution according to adjusted p-value
plot_chain_distribution(
  class_sm_AA_H,
  measure = "adj.P.Val"
)

plot_chain_distribution(
  class_sm_AA_M,
  measure = "adj.P.Val"
)

# Distribution according to log2 fold change
plot_chain_distribution(
  class_sm_AA_H,
  measure = "logFC"
)

plot_chain_distribution(
  class_sm_AA_M,
  measure = "logFC"
)

# =========================================================
# Lysophosphatidylcholines (LPC)
#
# Evaluate chain distributions for LPC species.
# =========================================================

class_lpc_AA_H <- et_c_AA_H %>%
  filter(Class == "LPC")

class_lpc_AA_M <- et_c_AA_M %>%
  filter(Class == "LPC")

# Distribution according to adjusted p-value
plot_chain_distribution(
  class_lpc_AA_H,
  measure = "adj.P.Val"
)

plot_chain_distribution(
  class_lpc_AA_M,
  measure = "adj.P.Val"
)

# Distribution according to log2 fold change
plot_chain_distribution(
  class_lpc_AA_H,
  measure = "logFC"
)

plot_chain_distribution(
  class_lpc_AA_M,
  measure = "logFC"
)

# =========================================================
# Phosphatidylcholines (PC)
#
# Evaluate chain distributions for PC species.
# =========================================================

class_pc_AA_H <- et_c_AA_H %>%
  filter(Class %in% c("PC"))

class_pc_AA_M <- et_c_AA_M %>%
  filter(Class %in% c("PC"))

# Distribution according to adjusted p-value
plot_chain_distribution(
  class_pc_AA_H,
  measure = "adj.P.Val"
)

plot_chain_distribution(
  class_pc_AA_M,
  measure = "adj.P.Val"
)

# Distribution according to log2 fold change
plot_chain_distribution(
  class_pc_AA_H,
  measure = "logFC"
)

plot_chain_distribution(
  class_pc_AA_M,
  measure = "logFC"
)

# =========================================================
# Phosphatidylethanolamines (PE)
#
# Evaluate chain distributions for PE species.
# =========================================================

class_pe_AA_H <- et_c_AA_H %>%
  filter(Class %in% c("PE"))

class_pe_AA_M <- et_c_AA_M %>%
  filter(Class %in% c("PE"))

# Distribution according to adjusted p-value
plot_chain_distribution(
  class_pe_AA_H,
  measure = "adj.P.Val"
)

plot_chain_distribution(
  class_pe_AA_M,
  measure = "adj.P.Val"
)

# Distribution according to log2 fold change
plot_chain_distribution(
  class_pe_AA_H,
  measure = "logFC"
)

plot_chain_distribution(
  class_pe_AA_M,
  measure = "logFC"
)

# =========================================================
# Triacylglycerols (TG)
#
# Evaluate chain-length and unsaturation distributions
# among triacylglycerol species.
# =========================================================

class_tg_AA_H <- et_c_AA_H %>%
  filter(Class == "TG")

class_tg_AA_M <- et_c_AA_M %>%
  filter(Class == "TG")

# Distribution according to adjusted p-value
plot_chain_distribution(
  class_tg_AA_H,
  measure = "adj.P.Val"
)

plot_chain_distribution(
  class_tg_AA_M,
  measure = "adj.P.Val"
)

# Distribution according to log2 fold change
plot_chain_distribution(
  class_tg_AA_H,
  measure = "logFC"
)

plot_chain_distribution(
  class_tg_AA_M,
  measure = "logFC"
)