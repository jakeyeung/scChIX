``` r
library(scChIX)
library(dplyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(igraph)
library(umap)
```

Model selection to infer cluster-pair for each double-incubated cells
---------------------------------------------------------------------

First we load the two pre-trained multinomial models (one for each mark)
inferred from the single-incubated data.

The model is a probability of getting a read within a region relative to
other regions.

This can be estimated from as simple as just fraction of reads in bin /
total reads for each cluster.

Here we ran Latent Dirichlet Allocation (Package: `topicmodels`) to
estimate this bin-specific probability genome-wide.

``` r
# data(ObjsForFittingSubset)
data(ModelFromSingleForFitting)
# Loading objects:
#   dat.impute.repress.lst
#   dat.impute.active

data(RawDblCountMatSubset)
# Loading objects:
#   count.mat.dbl.subset
```

Now for every double-incubated cell (here we take a subset of 90 cells
for speed), we perform model selection. The fitting infers which
cluster-pair (one from H3K27me3, one from H3K9me3) best fits the
observed double-incubated cut for each cell.

``` r
# bounds for fitting mixing fraction
wlower <- 0
wupper <- 1
# method for fitting mixing fraction
jmethod <- "Brent"
count.mat.dbl <- count.mat.dbl.subset
cell.count.raw.merged.lst <- as.list(as.data.frame(as.matrix(count.mat.dbl)))
system.time(
  act.repress.coord.lst <- lapply(cell.count.raw.merged.lst, function(cell.count.raw.merged){
    optim.out <- FitMixingWeight(cell.count.raw.merged = cell.count.raw.merged,
                                 dat.impute.repress.lst = dat.impute.repress.lst,
                                 dat.impute.active = dat.impute.active, w.init = 0.5, w.lower = wlower, w.upper = wupper, jmethod = "Brent")
    ll.mat <- GetLLMerged(optim.out$par, cell.count.raw.merged, dat.impute.repress.lst, dat.impute.active, return.mat = TRUE)
    return(list(ll.mat = ll.mat, w = optim.out$par))
  })
)
#>    user  system elapsed 
#>  73.732   0.499  74.274
```

We plot the fit for one cell as an example of the output. Here we show a
double-incubated granulocyte cell (which we know from ground truth from
FACS) is correctly predicted to come from a mixture of
granulocyte-specific H3K27me3 distribution plus a granulocyte-specific
H3K9me3 distribution.

``` r

(jcell.tmp <- names(act.repress.coord.lst)[[1]])
#> [1] "BM-EtOH-MNctrl-K27m3-K9m3-p1_220"
(jclst <- subset(annots.dat, cell == jcell.tmp)$celltype)
#> [1] "Granulocytes"

jcell.tmp <- "BM-EtOH-MNctrl-K27m3-K9m3-p1_220"
jclst <- "Granulocyte"  # ground truth

LL.mat.tmp <- act.repress.coord.lst[[jcell.tmp]]$ll.mat
LL.mat.tmp <- data.frame(H3K27me3 = rownames(LL.mat.tmp), LL.mat.tmp, stringsAsFactors = FALSE)
LL.long.tmp <- data.table::melt(LL.mat.tmp, id.vars = "H3K27me3", value.name = "LL", variable.name = "H3K9me3")

k27me3.ordering.all <- c("topic16", "topic3", "topic9", "topic25", "topic11")
k9me3.ordering.all <- c("topic12", "topic28", "topic1", "topic20")
ctype.ordering.all <- c("Granulocytes", "Bcells", "NKcells")

k27me3.clstr.labs <- list("topic16" = "Granulocytes",
                          "topic3" = "Bcells",
                          "topic9"  = "Bcells",
                          "topic25" = "NKcells",
                          "topic11" = "NKcells")
k9me3.clstr.labs <- list("topic12" = "Granulocytes",
                          "topic28" = "Bcells",
                          "topic1"  = "NKcells",
                          "topic20" = "NKcells")


LL.long.tmp$H3K27me3 <- factor(LL.long.tmp$H3K27me3, levels = k27me3.ordering.all)
LL.long.tmp$H3K9me3 <- factor(LL.long.tmp$H3K9me3, levels = k9me3.ordering.all)

m <- ggplot(LL.long.tmp, aes(x = H3K27me3, y = H3K9me3, fill = LL)) +
  geom_tile() +
  theme_bw() +
  scale_fill_viridis_c() +
  ggtitle(paste("Ground truth:", jclst, "\nCell name: ", jcell.tmp)) +
  theme(aspect.ratio= 4/8, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "bottom")
print(m)
```

<img src="scChIX-vignette_files/figure-markdown_github/visualize_one_cell-1.png" width="80%" />

For each double-incubated cell, we select the best cluster-pair inferred
from the model and summarize the assignments in a 2D grid.

``` r
# Summarize all cells  ----------------------------------------------------

cell.vec <- names(act.repress.coord.lst)
names(cell.vec) <- cell.vec
coords.dbl <- lapply(cell.vec, function(jcell){
  jfit <- act.repress.coord.lst[[jcell]]
  jweight <- act.repress.coord.lst[[jcell]]$w
  p.mat <- SoftMax(jfit$ll.mat)
  jcoord <- which(jfit$ll.mat == max(jfit$ll.mat), arr.ind = TRUE)
  jmax <- max(p.mat)

  jclstr.k27me3 <- rownames(p.mat)[[jcoord[[1]]]]
  jclstr.k9me3 <- colnames(p.mat)[[jcoord[[2]]]]

  if (grepl("_", jclstr.k27me3)){
    jclstr.k27me3 <- strsplit(jclstr.k27me3, split = "_")[[1]][[2]]
  }
  if (grepl("_", jclstr.k9me3)){
    jclstr.k9me3 <- strsplit(jclstr.k9me3, split = "_")[[1]][[2]]
  }

  out.dat <- data.frame(cell = jcell, clstr.k27me3 = jclstr.k27me3, clstr.k9me3 = jclstr.k9me3, lnprob = jmax, w = jweight, stringsAsFactors = FALSE)
  return(out.dat)
}) %>%
  bind_rows()


coords.dbl.annots <- left_join(coords.dbl, annots.dat)
#> Joining, by = "cell"
coords.dbl.annots$clstr.k27me3 <- factor(coords.dbl.annots$clstr.k27me3, levels = k27me3.ordering.all)
coords.dbl.annots$clstr.k9me3 <- factor(coords.dbl.annots$clstr.k9me3, levels = k9me3.ordering.all)

# annotate clusters
coords.dbl.annots$clstr.k27me3.annot <- sapply(as.character(coords.dbl.annots$clstr.k27me3), function(clst){
  k27me3.clstr.labs[[clst]]
})
coords.dbl.annots$clstr.k9me3.annot <- sapply(as.character(coords.dbl.annots$clstr.k9me3), function(clst){
  k9me3.clstr.labs[[clst]]
})


cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
m.grid <- ggplot(coords.dbl.annots, aes(x = clstr.k27me3.annot, y = clstr.k9me3.annot, color = celltype)) +
  geom_point(position = ggforce::position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  theme(aspect.ratio=1) +
  scale_color_manual(values = cbPalette) + xlab("H3K27me3 Cell Type Clusters") + ylab("H3K9me3 Cell Type Clusters") +
  ggtitle("Each dot is a double stained cell, colored by ground truth label.\nX-Y shows the cluster pair it is assigned")
print(m.grid)
```

<img src="scChIX-vignette_files/figure-markdown_github/visualize_all_cells-1.png" width="80%" />

Unmixing cut fragments in double-incubated cells to respective histone modification
-----------------------------------------------------------------------------------

For each double-incubated cell, we now have a model of how H3K27me3 and
H3K9me3 are mixed together to generate the observed double-incubated cut
fragments. We use this cell-specific model to assign reads in each
regions to either H3K27me3 or H3K9me3.

``` r

# Deconvolve double-incubated count matrix ------------------------------------------------------


all.cells <- coords.dbl$cell
names(all.cells) <- all.cells
col.i <- seq_len(ncol(count.mat.dbl))
names(col.i) <- colnames(count.mat.dbl)
all.x.raw <- lapply(col.i, function(i) count.mat.dbl[, i])  # https://stackoverflow.com/questions/6819804/how-to-convert-a-matrix-to-a-list-of-column-vectors-in-r/6823557

all.mixweights <- coords.dbl$w
names(all.mixweights) <- all.cells
all.clstr.k27me3 <- as.character(coords.dbl$clstr.k27me3)
all.clstr.k9me3 <- as.character(coords.dbl$clstr.k9me3)
names(all.clstr.k27me3) <- all.cells
names(all.clstr.k9me3) <- all.cells
all.p.active <- lapply(all.clstr.k27me3, function(clstr.k27me3) dat.impute.active[clstr.k27me3, ])
all.p.repress <- lapply(all.clstr.k9me3, function(clstr.k9me3) dat.impute.repress.lst[[clstr.k9me3]])


jnames.all <- lapply(X = list(all.x.raw, all.mixweights, all.p.active, all.p.repress), FUN = names)
assertthat::assert_that(all(sapply(jnames.all, identical, jnames.all[[1]])))
#> [1] TRUE
# assertthat::assert_that(all(names(all.x.raw) == names(all.mixweights)))

system.time(
  x.raw.unmixed <- lapply(all.cells, function(jcell){
    # print(jcell)
    return(UnmixRawCounts(x.raw = all.x.raw[[jcell]], mixweight = all.mixweights[[jcell]], p.active = all.p.active[[jcell]], p.repress = all.p.repress[[jcell]], random.seed = 0))
  })
)
#>    user  system elapsed 
#>   0.066   0.000   0.066
```

Projecting unmixed cells
------------------------

We take the H3K27me3 cuts deconvolved from double-incubated cells and
project them onto the H3K27me3 manifold previously learned from
single-incubated data. And similarly we project H3K9me3 cuts onto the
H3K9me3 manifold. We project using the Latent Dirichet Allocation
framework, which takes a vector of discrete counts and infers the
cell-to-topic weights (which tells you where the new cells will be
located in the UMAP), while freezing the topic-to-region weights (which
have already been learned from the single-incubated data).

We project only 10 cells onto the manifold to show how its done without
having to wait too long to finish while still running on a single core.

``` r

data(H3K27me3_LdaOutputs)  # out.lda.k27me3
data(H3K9me3_LdaOutputs)  # out.lda.k9me3

rnames <- rownames(count.mat.dbl)
x.k27me3.mat <- as.data.frame(lapply(x.raw.unmixed, function(outlst) return(outlst$x.raw.active)), row.names = rnames)
x.k9me3.mat <- as.data.frame(lapply(x.raw.unmixed, function(outlst) return(outlst$x.raw.repress)), row.names = rnames)

set.seed(0)
nsubset <- 10  # faster
rand.i1 <- sample(seq(ncol(x.k27me3.mat)), size = nsubset)
rand.i2 <- sample(seq(ncol(x.k9me3.mat)), size = nsubset)

print("Projecting to a few unmixed cells to K27me3 manifold")
#> [1] "Projecting to a few unmixed cells to K27me3 manifold"
system.time(
  out.lda.predict.k27me3 <- topicmodels::posterior(out.lda.k27me3, t(as.matrix(x.k27me3.mat[, rand.i1])))
)
#>    user  system elapsed 
#>  55.807   0.088  56.321

print("Projecting to a few unmixed cells to K9me3 manifold")
#> [1] "Projecting to a few unmixed cells to K9me3 manifold"
system.time(
  out.lda.predict.k9me3 <- topicmodels::posterior(out.lda.k27me3, t(as.matrix(x.k27me3.mat[, rand.i2])))
)
#>    user  system elapsed 
#>  78.330   0.004  78.351
```

Linking UMAPs from distinct histone modifications together
----------------------------------------------------------

We can visualize both UMAPs together now. We link the two UMAPs together
with the deconvoved double-incubated cells as connections between the
two UMAPs.

``` r
jsettings <- umap::umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

jmarks <- c("K27m3", "K9m3")
names(jmarks) <- jmarks

data(H3K27me3_ProjectionOutput)  # out.proj.K27m3
data(H3K9me3_ProjectionOutput)  # out.proj.K9m3

out.proj.lst <- list("K27m3" = out.proj.K27m3,
                     "K9m3" = out.proj.K9m3)

dat.umap.merge <- lapply(jmarks, function(jmark){
  out.proj.mark <- out.proj.lst[[jmark]]
  tm.result <- topicmodels::posterior(out.proj.mark$out.lda)

  topics.mat.orig <- tm.result$topics
  topics.mat.proj <- out.proj.mark$out.lda.predict$topics

  umap.out <- umap::umap(topics.mat.orig, config = jsettings)
  rownames(umap.out$layout) <- rownames(topics.mat.orig)
  dat.umap <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE)
  umap.proj <- predict(umap.out, topics.mat.proj)
  dat.umap.proj <- data.frame(cell = rownames(umap.proj), umap1 = umap.proj[, 1], umap2 = umap.proj[, 2], stringsAsFactors = FALSE)

  dat.umap.merge <- rbind(dat.umap %>% mutate(cond = "single"), dat.umap.proj %>% mutate(cond = "double")) %>%
    rowwise() %>%
    mutate(plate = ClipLast(cell, jsep = "_", jsep.out = "_")) %>%
    left_join(., annots.dat)
  dat.umap.merge$mark <- jmark
  return(dat.umap.merge)
})
#> Joining, by = "cell"
#> Joining, by = "cell"


# merge the single umaps into one, seaprate by umap2
dat.merge.singles <- dat.umap.merge %>%
  bind_rows() %>%
  group_by(mark) %>%
  mutate(umap1 = scale(umap1, center = TRUE, scale = TRUE),
         umap2 = scale(umap2, center = TRUE, scale = TRUE))

yrange <- diff(range(dat.merge.singles$umap2)) / 1.5
xrange <- diff(range(dat.merge.singles$umap1)) / 1.5

dat.merge.singles <- dat.merge.singles %>%
  rowwise() %>%
  mutate(umap2.shift = umap2,
         umap1.shift = ifelse(mark == jmarks[[1]], umap1 + xrange, umap1 - xrange))
```

``` r
m <- ggplot(dat.merge.singles, aes(x = -1 * umap1.shift, y = umap2.shift, color = celltype, shape = mark, group = cell)) +
  geom_point() +
  geom_path(size = 0.1, alpha = 0.1) +
  theme_bw() +
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = cbPalette) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  ggtitle(paste( "Linked map of H3K27me3 (left) and H3K9me3 (right)")) +
  xlab("UMAP1 (shifted)") + ylab("UMAP2")
print(m)
```

<img src="scChIX-vignette_files/figure-markdown_github/plot_umap2-1.png" width="80%" />
