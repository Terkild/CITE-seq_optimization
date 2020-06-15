# CITE-seq optimization
Code and results from TotalSeqC antibody titration and pipeline benchmarking for CITE-seq experiments.

Data and Seurat V3 objects from the study is deposited at FigShare under this DOI: [10.6084/m9.figshare.c.5018987](https://doi.org/10.6084/m9.figshare.c.5018987)

This repository contains all the code used in the processing of the aligned data and data analysis (used for generating all figures) included in the manuscript:

**Pre-processing:**
* [Loading data, Demultiplexing, Preprocessing and down-sampling](Demux_Preprocess_Downsample.md) - Supplementary Figure S1
* [Load unfiltered data and determine cell-containing vs. empty droplets](Load-unfiltered-data.md) - Supplementary Figure S6

**Data analysis:**
* [Antibody concentration titration](Antibody-titration.md) - Figure 1, 2 and Supplementary Figure S2
* [Reducing staining volume](Volume-titration.md) - Figure 3 and Supplementary Figure S3
* [Reducing cell number at staining](Cell-number-titration.md) - Figure 4 and Supplementary Figure S4
* [Reducing cell number mitigates reduced staining volume](Volume-and-cell-number-titration.md) - Supplementary Figure S5
* [ADT signal in cell-containing vs. empty droplets](ADT-reads-in-cells-vs-empty-drops.md) - Figure 5 and Supplementary Figure S8
* [10X Datasets: UMI per marker plots](10X-Datasets-UMI-per-marker.md) - Supplementary Figure S7
* [Comparison of ADT counting methods](ADT-counting-methods.md) - Figure 6 and Supplementary Figure S9

We also included the [Snakefiles](Snakemake/) used with Snakemake to generate the alignment and counting data from our dataset and for the 10X datasets.