# Lipid profiling of plasma extracellular vesicles reveals sex differences in acyl-chain length and degree of saturation in adolescents intoxicated with alcohol

## OVERVIEW

This repository contains the code and data used to characterize the structural lipid composition — acyl-chain length and degree of unsaturation — of plasma-derived extracellular vesicles (EVs) in adolescents intoxicated with alcohol compared with healthy controls. Lipidomic profiling was performed by LC-MS/MS, and analyses were stratified by sex to evaluate sex-specific patterns of lipid remodeling associated with chronic alcohol exposure.

## OBJECTIVE

To characterize the acyl-chain length and degree of saturation of lipid species present in plasma-derived EVs from AUD patients and healthy controls, and to identify sex-specific structural alterations across major lipid classes (fatty acids, sphingolipids, glycerophospholipids, and triacylglycerols).

**Analyses include:**

- Lipid annotation and harmonization of nomenclature
- Outlier detection and removal
- Exploratory data analysis (PCA and hierarchical clustering/heatmaps)
- Sex-stratified differential lipid abundance analysis (AUD vs. control, males and females independently)
- Structural analysis of total acyl-chain length and degree of unsaturation
- Visualization of chain-length/unsaturation distributions and class-level trends

## DATA

The dataset includes lipidomics measurements from plasma samples, acquired in positive- and negative-ionization mode and normalized prior to analysis:

```
humano_adolescents/
 |-DATOS_NORMALIZADOS_HUMANO_NEGATIVO_RESUMIDOS_SIN_OUTLIERS.csv
 |-DATOS_NORMALIZADOS_HUMANO_POSITIVO_RESUMIDOS_SIN_OUTLIERS.csv
```

Processed and filtered datasets (overall cohort and sex-stratified subsets), together with the corresponding sample metadata, are provided under `00Data/`.


**Script description:**

| Script | Description |
|---|---|
| `01_annotation_adolescents.R` | Lipid class annotation using RefMet and harmonization of lipid nomenclature |
| `02_exploratory_analysis_Heatmap_adolescents.R` | Hierarchical clustering and heatmap visualization of lipid abundance |
| `02_exploratory_analysis_PCA_adolescents.R` | Principal component analysis (PCA) of the lipidomic dataset |
| `03_lipidr_DAA_adolescents.R` | Sex-stratified differential abundance analysis (IM: males AUD vs. control; IF: females AUD vs. control) using the `lipidr` package |
| `04_plot_chain_distribution_adolescents.R` | Visualization of acyl-chain length and unsaturation distributions per lipid class |
| `05_plot_trend_class_adolescents.R` | LOESS-based trend analysis of log2 fold change (LFC) as a function of chain length and unsaturation, by lipid class |

Outputs of the differential abundance analysis (full and significant lipid species, males and females) are stored in `03Results/`.

## KEY METHODS

- Principal Component Analysis (PCA)
- Hierarchical clustering
- Heatmaps
- Lipid annotation using RefMet / LIPID MAPS
- Differential abundance analysis (`lipidr`, Benjamini-Hochberg-adjusted p-values)
- Structural lipid analysis (acyl-chain length and degree of unsaturation)

## PROJECT STRUCTURE

- **00Data/** — Raw and processed lipidomic datasets, together with sample metadata (overall cohort and sex-stratified subsets).
- **01Functions/** — Custom R functions used across the analysis scripts (PCA, outlier removal, hierarchical clustering visualization).
- **02Scripts/** — Numbered analysis pipeline, from lipid annotation to differential abundance and structural trend analysis.
- **03Results/** — Output objects from the differential abundance analysis, stratified by sex.

## ShinyApp

An interactive web application allowing exploration of the lipidomic dataset and differential abundance results is available at:

https://bioinfouvcipf.shinyapps.io/EVLipStrAdol/

## AUTHOR

**Blanca Martín Urdiales**
Bioinformatician & Biostatistician
PhD Candidate in Bioinformatics
Department of Physiology, School of Medicine and Dentistry, University of Valencia
