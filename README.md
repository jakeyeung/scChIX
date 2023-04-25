# scChIX
Repository for scChIX

2021-05-01

This is the R package for deconvolving multiplexed histone modifications ([https://www.biorxiv.org/content/10.1101/2021.04.26.440629v1](https://www.biorxiv.org/content/10.1101/2021.04.26.440629v1)). 


2023-01-01

More has now been added for the journal version of the paper: "scChIX-seq infers dynamic relationships between histone modifications in single cells" [https://www.nature.com/articles/s41587-022-01560-3](https://www.nature.com/articles/s41587-022-01560-3).

Specifically we've added in the vignettes: 

- simulation study to show how scChIX-seq performs on ground truth data. [Simulation studies and downstream analysis of estimates versus ground truth:](https://github.com/jakeyeung/scChIX/blob/main/vignettes/scChIX-simulation.md)

- application of scChIX-seq on early organogenesis, seeing relationships between heterochromatin and active transcription. [Comparing cell type and heterochromatin relationships, and global changes in ratios of two histone modifications:](https://github.com/jakeyeung/scChIX/blob/main/vignettes/scChIX-gastrulation.md)

- Uncover dynamic relationships between histone modifications by inferring two pseudotime labels during differentiation. [Inferring dual continuous pseudotime during macrophage differentiation](https://github.com/jakeyeung/scChIX/blob/main/vignettes/scChIX-macrophagedifferentiation.md)

```
library(devtools)
devtools::install_github("jakeyeung/scChIX")
```

Follow the vignettes for an example of deconvolving the multiplexed signal. The output is a linked UMAP. In the example it is H3K27me3 UMAP on the left, H3K9me3 UMAP on the right, with lines connecting the two maps together:

![Linked UMAP example:](example_umap.png)


## Vignettes for guided analysis of downstream analysis

Here are a list of vignettes for guided analysis:

[Simulation studies and downstream analysis of estimates versus ground truth:](https://github.com/jakeyeung/scChIX/blob/main/vignettes/scChIX-simulation.md)

[Analysis of scChIX downstream:](https://github.com/jakeyeung/scChIX/blob/main/vignettes/scChIX-vignette.md)

[Comparing cell type and heterochromatin relationships, and global changes in ratios of two histone modifications:](https://github.com/jakeyeung/scChIX/blob/main/vignettes/scChIX-gastrulation.md)

[Inferring dual continuous pseudotime during macrophage differentiation](https://github.com/jakeyeung/scChIX/blob/main/vignettes/scChIX-macrophagedifferentiation.md)
