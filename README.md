# 2021 neurosarcoid single-cell RNA-Seq + BCR/TCR sequencing

Published article: [The CSF in neurosarcoidosis contains consistent clonal expansion of CD8 T cells, but not CD4 T cells](https://pubmed.ncbi.nlm.nih.gov/35405431/)

Downloadable data, including raw UMI counts, TCR sequences, and BCR sequences: [FigShare Project](https://figshare.com/projects/2021_neurosarcoid_single-cell_RNA-Seq_BCR_TCR_sequencing/137667)

## Project description
These data were derived from individuals with active or inactive neurosarcoid (cases) and controls without neurosarcoid, but not necessarily without a diagnosed chronic disease. PBMCs were isolated from blood samples (PBMC) and cerebrospinal fluid (CSF) at the same visit. Relationship between samples, i.e. individual identifier, tissue, status, etc can be found in info/unified_project_info.tsv.

In order to recreate the analysis (which may not be **identical** at all points since a fixed seed wasn't used), download the single-cell data locally and run the R scripts in numerical order. The single-cell count matrices along with TCR and BCR sequencing results are in a single gzip tar file from the downloadable data link above. If you download this file into the data subdirectory and extract, you should have everything you need.

## Code description
This repository contains the appropriate layout for recreating the analysis performed in the listed paper. All analysis was performed in R, and the scripts are in the analysis sub-directory. The 'here' package is used extensively so that the scripts should work across systems.

### Required R packages
* cowplot
* DoubletFinder
* dplyr
* escape
* ggplot2
* grid
* here
* monocle
* patchwork
* pheatmap
* purrr
* RColorBrewer
* readr
* rlang
* scales
* sctransform
* Seurat
* SeuratWrappers
* stringr
* tibble
* tidyr
* tidyverse
* VGAM
