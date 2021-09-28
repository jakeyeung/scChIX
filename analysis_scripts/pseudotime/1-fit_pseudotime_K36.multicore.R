# Jake Yeung
# Date of Creation: 2021-08-07
# File: ~/projects/scChIX/analysis_scripts/pseudotime/6-fit_pseudotime_K36.multicore.R
#

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(topicmodels)

library(hash)
library(igraph)
library(umap)

library(scchicFuncs)
library(scChIX)
library(JFuncs)

# Functions ---------------------------------------------------------------






# Constants ---------------------------------------------------------------

ncores <- 16

jmark <- "K36"
outmain <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/glm_fits_outputs"
jprefix <- "var_filtered"
jsuffix <- "manual2nocenter_K36_K9m3_K36-K9m3"
jname <- paste(jprefix, jsuffix, sep = "_")
outdir <- file.path(outmain, jsuffix)
dir.create(outdir)

outf <- file.path(outdir, paste0("glm_poisson_fits_output.", jsuffix, ".", jmark, ".", Sys.Date(), ".RData"))
outpdf <- file.path(outdir, paste0("pseudotime_plots.", jsuffix, ".", jmark, ".", Sys.Date(), ".pdf"))

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123


# Check whether run on projection or on mixed -----------------------------

# projection
# /home/jyeung/hub_oudenaarden

inf.meta.proj <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/scchix_downstream_plots/celltyping_after_scchix/var_filtered_manual2nocenter_K36_K9m3_K36-K9m3/celltyping_", jmark, "_first_try.2021-08-02.txt")
assertthat::assert_that(file.exists(inf.meta.proj))

dat.meta.proj <- fread(inf.meta.proj)

m1 <- ggplot(dat.meta.proj, aes(x = umap1, y = umap2, color = cluster)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# mixed
# inf.mixed <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/LDA_scchix_outputs/from_pipeline_unmixed_singles_LDA_together/var_filtered_manual2nocenter_K36_K9m3_K36-K9m3/lda_outputs.scchix_inputs_clstr_by_celltype_K36-K9m3.removeNA_FALSE-merged_mat.K9m3.K-30.binarize.FALSE/ldaOut.scchix_inputs_clstr_by_celltype_K36-K9m3.removeNA_FALSE-merged_mat.", jmark, ".K-30.Robj")
inf.mixed <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/LDA_scchix_outputs/from_pipeline_unmixed_singles_LDA_together/var_filtered_manual2nocenter_K36_K9m3_K36-K9m3/lda_outputs.scchix_inputs_clstr_by_celltype_K36-K9m3.removeNA_FALSE-merged_mat.", jmark, ".K-30.binarize.FALSE/ldaOut.scchix_inputs_clstr_by_celltype_K36-K9m3.removeNA_FALSE-merged_mat.", jmark, ".K-30.Robj")
assertthat::assert_that(file.exists(inf.mixed))

load(inf.mixed, v=T)

tm.result <- posterior(out.lda)

dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings = jsettings)

m2 <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

cell2ctype <- hash::hash(dat.meta.proj$cell, dat.meta.proj$cluster)

# compare with dat.umap
dat.umap$ctype <- sapply(dat.umap$cell, function(x) cell2ctype[[x]])

m3 <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = ctype)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Fit on the unmixed data  -----------------------------------------------

# load ptime from K9 fit
inf.ptime <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/metadata/from_pseudotime/metadata_pseudotime.2021-08-05.txt"

dat.ptime <- fread(inf.ptime)

cell2ptime <- hash::hash(dat.ptime$cell, dat.ptime$ptime)

dat.umap$ptime <- sapply(dat.umap$cell, function(x) AssignHash(x, cell2ptime, NA))

m5 <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = ptime)) +
  geom_point() +
  scale_color_viridis_c() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Fit genes that follow ptime  --------------------------------------------------------------

# load raw counts

cells.keep <- subset(dat.umap, !is.na(ptime))$cell
dat.umap.filt <- subset(dat.umap, cell %in% cells.keep)
count.mat.filt <- count.mat[, cells.keep]

# fit poisson regression
dat.annots.filt <- subset(dat.umap.filt, cell %in% cells.keep, select = c(cell, ptime))
ncuts.cells <- data.frame(cell = colnames(count.mat.filt), ncuts.total = colSums(count.mat.filt), stringsAsFactors = FALSE)

# jrow <- count.mat.filt[1, ]
# jbin <- rownames(count.mat.filt)[1]
# dat.fit.glm <- FitGlmRowPtime.withse(jrow = jrow, cnames = cells.keep, dat.annots.filt.mark = dat.annots.filt, ncuts.cells.mark = ncuts.cells, jbin =jbin , returnobj = FALSE, with.se = TRUE)
#
# fit all genes

jrow.names <- rownames(count.mat.filt)
names(jrow.names) <- jrow.names

print("fitting genes")
system.time(
  jfits.lst <- parallel::mclapply(jrow.names, function(jrow.name){
    jrow <- count.mat.filt[jrow.name, ]
    jout <- scChIX::FitGlmRowPtime.withse(jrow = jrow, cnames = cells.keep, dat.annots.filt.mark = dat.annots.filt, ncuts.cells.mark = ncuts.cells, jbin = jrow.name , returnobj = FALSE, with.se = TRUE)
    return(jout)
  }, mc.cores = ncores)
)

# Ssave outputs -----------------------------------------------------------

save(jfits.lst, dat.annots.filt, ncuts.cells, count.mat.filt, dat.umap.filt, file = outf)

print("Making plots")

pdf(outpdf, useDingbats = FALSE)
print(m1)
print(m2)
print(m3)
print(m5)
dev.off()
