# Jake Yeung
# Date of Creation: 2021-08-12
# File: ~/projects/scChIX/analysis_scripts/pseudotime/3-fit_celltypes_K9m3.multicore.discrete.R
# Fit K9me3 data but use K36 celltype annotations



rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scChIX)

ncores <- 32

# Functions ---------------------------------------------------------------



# Constants ---------------------------------------------------------------



hubprefix <- "/home/jyeung/hub_oudenaarden"

jmark <- "K9m3"
jname <- "manual2nocenter_K36_K9m3_K36-K9m3"

outdir <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis/glm_fits_outputs/by_clusters", jname)
dir.create(outdir)
outf <- file.path(outdir, paste0("glm_poisson_fits_output.clusters.", jname, ".", jmark, ".RData"))
print(outf)
assertthat::assert_that(!file.exists(outf))


# Load rawcounts  --------------------------------------------------------------

inf.obj <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/LDA_scchix_outputs/from_pipeline_unmixed_singles_LDA_together/var_filtered_", jname, "/lda_outputs.scchix_inputs_clstr_by_celltype_K36-K9m3.removeNA_FALSE-merged_mat.", jmark, ".K-30.binarize.FALSE/ldaOut.scchix_inputs_clstr_by_celltype_K36-K9m3.removeNA_FALSE-merged_mat.", jmark, ".K-30.Robj"))
assertthat::assert_that(file.exists(inf.obj))

load(inf.obj, v=T)

# count.mat[1:5, 1:5]

# # marks sometimes contain other marks? doublecheck
# inf.mat <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/count_tables_dedup_from_merged/counts_tables_10000/gastru_merged_", jmark, ".rowsdeduped.tagged.sorted.countTable.binsize_10000.csv"))
# mat <- ReadMatSlideWinFormat(inf.mat)

# Load celltype annotations -----------------------------------------------

inf.ctype <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis/metadata/from_scchix_downstream_plots/metadata_K36-K9me3.", jmark, ".2021-08-09.txt"))
dat.ctype <- fread(inf.ctype)

cells.keep <- dat.ctype$cell

ggplot(dat.ctype, aes(x = umap1, y = umap2, color = cluster)) +
  geom_point() +
  facet_wrap(~type) +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# make Epithelial the reference celltype

dat.ctype <- dat.ctype %>%
  rowwise() %>%
  mutate(cluster = ifelse(cluster == "Epithelial", "aEpithelial", cluster),
         cluster = ifelse(startsWith(cluster, "NeuralTubeNeuralProgs"), "NeuralTubeNeuralProgs", cluster))

m.umap <- ggplot(dat.ctype, aes(x = umap1, y = umap2, color = cluster)) +
  geom_point() +
  facet_wrap(~type) +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# rename NeuralTubeNeuralProgs2 and 3 to no numbers


# Get norm constants ------------------------------------------------------


# Fit each gene  ---------------------------------------------------------------


cnames <- colnames(count.mat)
dat.annots.filt.mark <- dat.ctype
ncuts.cells.mark <- data.frame(cell = colnames(count.mat), ncuts.total = colSums(count.mat), stringsAsFactors = FALSE)

# # test on one
# set.seed(0)
# jrow.i <- sample(seq_len(nrow(count.mat)), size = 1)
# jbin <- rownames(count.mat)[jrow.i]
# jrow <- count.mat[jrow.i, ]
# jfit.out <- FitGlmRowClusters.withse(jrow, cnames, dat.annots.filt.mark, ncuts.cells.mark, jbin = jbin, returnobj = FALSE, with.se = TRUE)

# fit all


jrow.names <- rownames(count.mat)
names(jrow.names) <- jrow.names

print("fitting genes")
system.time(
  jfits.lst <- parallel::mclapply(jrow.names, function(jrow.name){
    jrow <- count.mat[jrow.name, ]
    # jout <- scChIX::FitGlmRowPtime.withse(jrow = jrow, cnames = cells.keep, dat.annots.filt.mark = dat.annots.filt, ncuts.cells.mark = ncuts.cells, jbin = jrow.name , returnobj = FALSE, with.se = TRUE)
    jout <- scChIX::FitGlmRowClusters.withse(jrow, cnames, dat.annots.filt.mark, ncuts.cells.mark, jbin = jrow.name, returnobj = FALSE, with.se = TRUE)
    return(jout)
  }, mc.cores = ncores)
)

save(jfits.lst, dat.annots.filt.mark, ncuts.cells.mark, count.mat, dat.ctype, file = outf)



