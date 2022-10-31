# poisson.mash.paper
This repository provides R scripts to perform the analyses in the poisson mash paper (Liu et al, 2022).

### Code related to simulation study
The scripts related to the simulation study are located in the folder ./simulations.
In addition to the implementation of poisson mash which calls the R package poisson.mash.alpha,
we also provide R scripts to implement all competing methods considered in Liu et al (2022). <br />

### Code related to application to the single cell cytokines dataset
The scripts related to application to the single cell cytokines data are located in the folder ./applications.
We provide R scripts to apply poisson mash to scRNA-seq data for neutrophils, which calls the R package poisson.mash.alpha,
and analyze the results. Application of poisson mash to scRNA-seq data from other cell types is essentially the same.

### Dependency on other R packages
To implement poisson mash, the following R packages need to be installed: poisson.mash.alpha, glmpca, scran. 


## Reference

Yusha Liu, Peter Carbonetto, Michihiro Takahama, Adam Gruenbaum, Dongyue Xie, Nicolas Chevrier, and Matthew Stephens. (2022).
A flexible model for correlated count data, with application to analysis of gene expression differences in multi-condition experiments.
<https://arxiv.org/abs/2210.00697>
