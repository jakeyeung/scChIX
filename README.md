# scChIX
Repository for scChIX

This is a simple R package for deconvolving multiplexed histone modifications ([https://www.biorxiv.org/content/10.1101/2021.04.26.440629v1](https://www.biorxiv.org/content/10.1101/2021.04.26.440629v1)). 

```
library(devtools)
devtools::install_github("jakeyeung/scChIX")
```

Follow the vignettes for an example of deconvolving the multiplexed signal. The output is a linked UMAP. In the example it is H3K27me3 UMAP on the left, H3K9me3 UMAP on the right, with lines connecting the two maps together:

![Linked UMAP example:](example_umap.png)
