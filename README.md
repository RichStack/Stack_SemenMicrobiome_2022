# Low-biomass 16S amplicon workflow (QIIME2 → R)

This repository contains an R-based downstream analysis workflow
for low-microbial biomass 16S rRNA amplicon data.

## Input requirements
- QIIME2 feature table (`table.qza`)
- Rooted phylogenetic tree
- Taxonomy assignments (SILVA and/or Greengenes)
- Sample metadata (TSV)

## Assumptions
- Data processed with QIIME2 2024+
- DADA2 denoising
- Paired-end Illumina reads
- Low-biomass samples with negative controls

## Quick start
1. Clone repository
2. Copy `config/config_example.R` → `config/config.R`
3. Edit paths
4. Run `analysis/main_analysis.R`
   
## Acknowledgements

This workflow draws on ideas and code patterns from:
- Callahan et al. DADA2
- Davis et al. decontam
- F1000Research microbiome analysis guidelines

Any mistakes or adaptations are my own.
