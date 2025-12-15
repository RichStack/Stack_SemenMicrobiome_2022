# 16S rRNA amplicon analysis workflow (QIIME2 → R)

This repository contains an R-based downstream analysis workflow for
16S rRNA amplicon sequencing data, developed during my PhD research.

The workflow is designed primarily for **host-associated low-microbial biomass samples**
and makes extensive use of negative controls and mock communities to
assess contamination, detection limits, and data plausibility.

It also makes use of environmental controls - i.e. samples from adjacent biological niches
in order to compare taxonomic and phylogenetic profiles of the different sample-types.

This is a research workflow rather than a polished software package,
and is shared for transparency, reuse, and reproducibility.

---

## Overview of the workflow

The analysis assumes that upstream processing has already been performed
using QIIME2. The R script performs the following major steps:

1. **Import of QIIME2 outputs**
   - Feature table
   - Taxonomy assignments (SILVA and/or Greengenes)
   - Rooted phylogenetic tree
   - Sample metadata
   - Construction of a `phyloseq` object

2. **Initial quality control**
   - Removal of non-bacterial features
   - Filtering of implausible or low-prevalence taxa
   - Basic inspection of sequencing depth and feature counts

3. **Contamination assessment**
   - Identification and removal of contaminant ASVs using `decontam`
   - Use of negative controls where available

4. **Mock community analysis**
   - Evaluation of mock community composition
   - Determination of minimum detection thresholds using:
     - Cell-based dilution series
     - DNA-based logarithmic mock standards

5. **Final filtering**
   - Removal of features below empirically determined thresholds
   - Generation of the final analysis-ready dataset
   - Analysis of alpha rarefaction and sampling depth

6. **Diversity analyses**
   - Alpha diversity
   - Beta diversity (ordination-based analyses)

7. **Taxonomic summaries**
   - Bar plots and other relative abundance summaries
   - Taxonomic composition across sample groups

8. **Comparative niche analysis**
   - Comparison of target samples to environmental or adjacent niches

9. **Phylogenetic visualisation**
   - Phylograms and tree-based representations of selected taxa

---

## Input requirements

The workflow expects the following inputs:

- QIIME2 feature table (`.qza`)
- Rooted phylogenetic tree (`.qza`)
- Taxonomy assignments (SILVA and/or Greengenes)
- Sample metadata file (TSV format)

Details on required columns in the metadata file are described in the
script comments.

---

## Repository structure

```text
analysis/
  main_analysis.R      # Main analysis script

config/
  config_example.R     # Example configuration file (paths & parameters)

data/
  README.md            # Place input data here (not tracked)

output/
  README.md            # Analysis outputs are written here

   
## Acknowledgements

This workflow draws on ideas and code patterns from:
- Callahan et al. DADA2
- Davis et al. decontam
- F1000Research microbiome analysis guidelines

Any mistakes or adaptations are my own.
