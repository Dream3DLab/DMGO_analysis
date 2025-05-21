[![DOI](https://zenodo.org/badge/15356014.svg)](https://zenodo.org/badge/latestdoi/15356014)

# Transcriptomic profiling of BrO, DMGO and GD2 CAR T cells


> **Note** :construction: <br>
> We are currently collecting and uploading all notebooks used in this analysis.<br>
> This note wil disappear as soon as all code has been made available via this page.

Repository accompanying the publication **De novo H3.3K27M-altered Diffuse Midline Glioma in human brainstem organoids to dissect GD2 CAR T cell function**

![](supplemental_files/graphical_abstract.jpg)

The repository is organized firstly per model and secondly per coding language

* `code/BrO/`: Code, jupyter notebooks and .Rmd files used for processing and analyzing the time course sequencing data of the Brainstem Organoid, including annotation, integration and reference mapping against the [HNOCA](https://doi.org/10.1038/s41586-024-08172-8) and the [HDBCA](https://doi.org/10.1126/science.adf1226).
* `code/DMGO/`: Code, jupyter notebooks and .Rmd files used for processing and analyzing de novo generated DMG tumors in BrO, including annotation and reference mapping against [Liu et al., 2022](https://doi.org/10.1038/s41588-022-01236-3) for DMG specific analysis and pan-tumor comparisons.
* `code/CART/`: Code, .Rmd files and R scripts used for processing and analyzing each CART sequencing batch.

For code related to the lineage tracing data please refer to the [GitHub](https://github.com/anna-alemany/TrackerSeq_BROs) of the [Alemany lab](https://www.alemany-lab.com/).

## Data availability
Raw fastq files have been deposited on ArrayExpress under accession [E-MTAB-15147]. <br>
Processed datasets (.h5ad and .rds) are available through [Zenodo](https://zenodo.org/badge/latestdoi/15356014).

## Citation
> Bessler, N.; Wezenaar A.K.L.; Ariese, H.C.R.; Honhoff, C. et al. De novo H3.3K27M-altered Diffuse Midline Glioma in human brainstem organoids to dissect GD2 CAR T cell function. Nature Cancer, 2025

## Acknowledgements

We thank [Cristian Ruiz Moreno](https://github.com/ccruizm) for all the help and feedback whenever we needed.  
Processing .Rmd scripts have been inspired by [Selin Jessa](https://github.com/sjessa).  