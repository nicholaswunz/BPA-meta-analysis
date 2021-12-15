# Plastic pollutant BPA in aquatic organisms

This repository contains code and data needed to reproduce the article:

**Wu N. C., & Seebacher, F.** (2020) Effect of the plastic pollutant bisphenol A on the biology of aquatic organisms: A meta-analysis. *Global Change Biology*, **26**, 3821-3833. DOI:
[![DOI](https://zenodo.org/badge/DOI/10.1111/gcb.15127.svg)](https://doi.org/10.1111/gcb.15127)

**Raw data**
- meta data raw.csv    - Meta-data used for the meta-analysis.
- BPA level enviro.csv - Environmental BPA data used for the analysis.
- BPA level taxa.csv   - BPA in aquatic organisms data used for the analysis.
- world_2.csv          - Global BPA and plastic level data used for the analysis.
- phylo_cor_all.RDAT   - Phylogenetic co-variance matrix for the meta-analysis.

**R codes**
- meta_analysis.R      - Meta-analysis, meta-regressions and figure production.
- meta_phylo.R         - Build phylogeny and estimate phylogenetic relatedness among species.
- BPA levels.R         - Plotting BPA levels in taxa and the envrironment, produce global distribution map and analysis.

**Extra files**
- GCB-20-0299_SI.pdf - Supplementary file includes statistical outcomes and additional figures and descriptions from the main document.

## Abstract
Plastic pollution is a global environmental concern. In particular, the endocrine-disrupting chemical (EDC) bisphenol A (BPA) is nearly ubiquitous in aquatic environments globally, and it continues to be produced and released into the environment in large quantities. BPA disrupts hormone signalling and can thereby have far-reaching physiological and ecological consequences. However, it is not clear whether BPA has consistent effects across biological traits and phylogenetic groups. Hence, the aim of this study was to establish the current state-of-knowledge of the effect of BPA in aquatic organisms. We show that overall BPA exposure affected aquatic organisms negatively. It increased abnormalities, altered behaviour, and had negative effects on the cardiovascular system, development, growth, and survival. Early life stages were the most sensitive to BPA exposure in invertebrates and vertebrates, and invertebrates and amphibians seem to be particularly affected. These data provide a context for management efforts in the face of increasing plastic pollution. However, data availability is highly biased with respect to taxonomic groups and traits studies, and in the geographical distribution of sample collection. The latter is the case for both measurements of the biological responses and assessing pollution levels in water ways. Future research effort should be directed towards biological systems, such as studying endocrine disruption directly, and geographical areas (particularly in Africa and Asia) which we identify to be currently undersampled.

**Keywords:** BPA, dose-response, ecotoxicology, endocrine disrupting chemical (EDC), survival, reproduction
