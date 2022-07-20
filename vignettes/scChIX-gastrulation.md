-   [Set up gastrulation data for scChIX snakemake
    workflow](#set-up-gastrulation-data-for-scchix-snakemake-workflow)
-   [Run snakemake workflow](#run-snakemake-workflow)
-   [Uncover cell type and heterochromatin
    relationships](#uncover-cell-type-and-heterochromatin-relationships)

## Set up gastrulation data for scChIX snakemake workflow

First, we load the count matrices that are used for running the scChIX
snakemake workflow.

``` r
data(CountMatsGastrulationInputs)
# Loading objects:
#   countmats.gastru
```

Save these count matrices as `.rds` files into a
`snakemake_inputs/countmats` directory where the `Snakefile`,
`config.yaml`, and `cluster.json` files are.

``` r
# copy snakemake workflow files into directory for outputs
outdir <- "/tmp"
copycmd=paste0("cp snakemake_workflow/Snakefile snakemake_workflow/cluster.json snakemake_workflow/config.yaml snakemake_workflow/run_snakemake.gastru_K9.sh", outdir)
# system(cpycmd)  # run this to move files

# save matrices as .rds to snakemake_inputs
mkdircmd=paste0("mkdir -p ", outdir, "/snakemake_inputs/countmats")
# system(mkdircmd)  # make directory before writing mats
# run chunk to save countmats as rds into snakemake input directory
# use countmat_var_filt.${mark}.rds as name to match expected filename in Snakemake workflow
# jmarks <- names(countmats.gastru); names(jmarks)
# for (jmark in jmarks){
#   saveRDS(countmats.gastru[[jmark]], file.path(outdir, paste0("countmat_var_filt.", jmark, ".rds")))
# }
print(outdir)
#> [1] "/tmp"
```

## Run snakemake workflow

``` r
bashscript=file.path(outdir, paste0("run_snakemake.gastru_K9.sh")) # modify for your specific conda environment and HPC cluster settings
runcmd=paste0("bash ", bashscript)
# system(runcmd)  # launch snakemake workflow
```

## Uncover cell type and heterochromatin relationships

Let’s load the outputs of the snakemake workflow and look at the cell
type and heterochromatin relationships.

``` r
data(GastrulationScChIXOutputsK36K9m3)
# Loading objects:
#   act.repress.coord.lst
fits.out <- act.repress.coord.lst

data(GastrulationLouvainCelltypeAnnotations)
# Loading objects:
#   louvain.celltype.metadata
#   ctype.colcode.metadata
louv2ctype.act <- hash::hash(louvain.celltype.metadata$K36$louv.act, louvain.celltype.metadata$K36$cluster)
louv2ctype.repress <- hash::hash(louvain.celltype.metadata$K9m3$louv.repress, louvain.celltype.metadata$K9m3$cluster)
ctype2colcode <- hash::hash(ctype.colcode.metadata$cluster, ctype.colcode.metadata$colorcode)
```

We wrangle the fits and then plot the celltype to heterochromatin
relationship

``` r
# if louvains are now from clusters need eto rethink jcoord
cell.vec <- names(fits.out)
names(cell.vec) <- cell.vec
coords.dbl <- lapply(cell.vec, function(jcell){
  jfit <- fits.out[[jcell]]
  jweight <- fits.out[[jcell]]$w
  p.mat <- SoftMax(jfit$ll.mat)
  jcoord <- which(jfit$ll.mat == max(jfit$ll.mat), arr.ind = TRUE)
  jmax <- max(p.mat)
  
  jlouv.act <- rownames(p.mat)[[jcoord[[1]]]]
  jlouv.repress <- colnames(p.mat)[[jcoord[[2]]]]
  jcluster.act <- louv2ctype.act[[jlouv.act]]
  jcluster.repress <- louv2ctype.repress[[jlouv.repress]]
  colcode <- ctype2colcode[[jcluster.act]]
  
  if (grepl("_", jlouv.act)){
    jlouv.act <- strsplit(jlouv.act, split = "_")[[1]][[2]]
  }
  if (grepl("_", jlouv.repress)){
    jlouv.repress <- strsplit(jlouv.repress, split = "_")[[1]][[2]]
  }
  out.dat <- data.frame(cell = jcell, celltype.act = jcluster.act, celltype.repress = jcluster.repress, colcode = colcode, lnprob = jmax, w = jweight, stringsAsFactors = FALSE)
  return(out.dat)
}) %>%
  bind_rows()


clstrs.order.active <- c("Erythroid", "WhiteBloodCells", "Endothelial", "NeuralTubeNeuralProgs", "Neurons", "SchwannCellPrecursor", "Epithelial", "MesenchymalProgs", "Cardiomyocytes")
clstrs.order.repress <- c("Erythroid", "WhiteBloodCells", "NonBlood")

coords.dbl$celltype.act <- factor(coords.dbl$celltype.act, levels = clstrs.order.active)
coords.dbl$celltype.repress <- factor(coords.dbl$celltype.repress, levels = clstrs.order.repress)

m.grid <- ggplot(coords.dbl, aes(x = celltype.act, y = celltype.repress, color = colcode)) +
  geom_point(alpha = 0.25, position = ggforce::position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  scale_color_identity() + 
  theme(aspect.ratio=0.6, axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  ggtitle("Each dot is a double stained cell,\nX-Y shows the cluster pair it is assigned")
print(m.grid)
```

<img src="scChIX-gastrulation_files/figure-markdown_github/wrangle-scchix-1.png" width="80%" />

The scChIX output also reveals global ratios of H3K36me3 and H3K9me3 and
how they are lower in erythroids versus other cell types.

``` r
m.ratios <- ggplot(coords.dbl, aes(x = celltype.act, y = log2(w / (1 - w)), fill = colcode)) +
  geom_boxplot() + 
  theme_bw() +
  scale_fill_identity() + 
  ylab("log2 H3K36me3 to H3K9me3 ratio") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("log2 H3K36me3 to H3K9me3 ratio inferred from double-incubated cells")
print(m.ratios)
```

<img src="scChIX-gastrulation_files/figure-markdown_github/ratios-1.png" width="80%" />