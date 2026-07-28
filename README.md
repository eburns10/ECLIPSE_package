# ECLIPSE
ECLIPSE is a novel method of scTCR-seq analysis that can predict the 
clonotypes of T cells possessing >2 or only 1 TCR chain. These populations of 
cells, which can represent >50% of cells in some samples, are not well 
suited for traditional scTCR-seq analysis. By mapping chain pairings across the 
sample and using the EM algorithm, ECLIPSE is able to predict chains lost to 
sequencing dropout and discern between cells that have >2 TCR due to true  
biology vs. sequencing artifacts. Together this leads to larger clone sizes 
and less cells without a clonotype assignment, allowing users to make best use 
of their scTCR-seq data.

### Authors
Ethan C. Burns, Zhaochen Ye, Kelly Street, and David A. Braun

### Installation
To download the development version of ECLIPSE, please run:
``` r
remotes::install_github("eburns10/ECLIPSE")
```

### Manuscript
ECLIPSE is described in detail in our pre-print on 
bioRxiv: "VDJdive and ECLIPSE enhance single-cell TCR sequencing analysis 
through the probabilistic resolution of ambiguous clonotypes" available 
[here](https://doi.org/10.64898/2026.02.18.706444)

### Issues
Please use our issues 
[page](https://github.com/eburns10/ECLIPSE_bioconductor_submission/issues) or 
email ethan.burns [@] yale.edu

### How to use
Please view the vignette in the vignettes folder


