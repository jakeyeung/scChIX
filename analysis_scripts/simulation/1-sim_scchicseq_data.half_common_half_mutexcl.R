# Jake Yeung
# Date of Creation: 2021-09-20
# File: ~/projects/scChIX/analysis_scripts/simulation/1-sim_scchicseq_data.half_common_half_mutexcl.R
#



rm(list=ls())


# library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(irlba)
library(hash)
library(igraph)
library(umap)
library(scchicFuncs)

# library(devtools)
# install_github("bowang-lab/simATAC")

library(simATAC)

# Constants ---------------------------------------------------------------


hubprefix <- "/home/jyeung/hub_oudenaarden"
outdir <- file.path(hubprefix, "jyeung/data/dblchic/simulation_data/frac_mutexcl_0.5/snakemake_inputs/countmats")
dir.create(outdir, recursive = TRUE)

jspec <- "mm10"
jseed <- 0
jncells.per.rep <- 250
# jnbins <- 5000
jnbins <- 10000
jlibmean <- 12
jlibsd <- 1
jlibp <- 0.5
jzp <- 1
# frac.mutual.excl <- 0.5
shuffle.rows.seed <- jseed

sim.params <- list(jspec = jspec,
                   jseed = jseed,
                   jncells.per.rep = jncells.per.rep,
                   jnbins = jnbins,
                   jlibmean = jlibmean,
                   jlibsd = jlibsd,
                   jlibp = jlibp,
                   jzp = jzp)


# Run simulation ----------------------------------------------------------

ctypes <- c("A", "B", "C")
names(ctypes) <- ctypes

ctype.params <- list(ctype.name = c("A", "B", "C"),
                     jseed = c(0, 123, 999),
                     shuffle.rows.seed = c(0, 123, 999),
                     frac.mutual.excl = c(0.5, 0.5, 0.5))

ctype.params.lst <- lapply(ctypes, function(ctype){
  i <- which(ctype.params$ctype.name == ctype)
  return(list(ctype = ctype,
              jseed = ctype.params$jseed[[i]],
              shuffle.rows.seed = ctype.params$shuffle.rows.seed[[i]],
              frac.mutual.excl = ctype.params$frac.mutual.excl[[i]]))
})

ctype.sim.counts.lst <- lapply(ctypes, function(jctype){
  SimulateChICseq(ctype.name = jctype,
                  jseed = ctype.params.lst[[jctype]]$jseed,
                  shuffle.rows.seed = ctype.params.lst[[jctype]]$shuffle.rows.seed,
                  frac.mutual.excl = ctype.params.lst[[jctype]]$frac.mutual.excl)
})


# Create simulated matrix -------------------------------------------------


jmark1 <- "mark1"
jmark2 <- "mark2"
jmarkdbl <- paste(jmark1, jmark2, sep = "-")

mat.mark1 <- ConcatMat(ctype.sim.counts.lst, ctypes, jmark = "mark1")
mat.mark2 <- ConcatMat(ctype.sim.counts.lst, ctypes, jmark = "mark2")
mat.markdbl <- ConcatMat(ctype.sim.counts.lst, ctypes, jmark = "markdbl")

mat.mark.lst <- list(mat.mark1, mat.mark2, mat.markdbl)
names(mat.mark.lst) <- c(jmark1, jmark2, jmarkdbl)

jmarks <- names(mat.mark.lst)
names(jmarks) <- jmarks

jsettings <- umap.defaults
jsettings$n_neighbors <- 150
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

dat.umap.lst <- lapply(mat.mark.lst, function(jmat){
  print(head(jmat[1:5, 1:5]))
  lsi.out <- scchicFuncs::RunLSI(as.matrix(jmat))
  dat.umap.mark <- scchicFuncs::DoUmapAndLouvain(topics.mat = lsi.out$u, jsettings = jsettings)
  return(dat.umap.mark)
})

m.lst <- lapply(dat.umap.lst, function(jdat){
  m <- ggplot(jdat, aes(x = umap1, y = umap2, color = louvain)) +
    geom_point() +
    theme_bw() +
    scale_color_manual(values = cbPalette) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(m)
})

print(m.lst)

# Write to output ---------------------------------------------------------

sim.files <- file.path(outdir, paste0("ATAC_simulator_params_and_outputs.RData"))

save(ctype.params.lst, ctype.sim.counts.lst, file = sim.files)

# save each mark
for (jmark in jmarks){
  print(jmark)
  outf <- file.path(outdir, paste0("countmat_var_filt.", jmark, ".rds"))
  saveRDS(mat.mark.lst[[jmark]], file = outf)
}



# Checks  -----------------------------------------------------------------


dat.meta.lst <- lapply(jmarks, function(jmark){
  dat.meta.tmp <- data.frame(cell = colnames(mat.mark.lst[[jmark]]), stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(ctype = substr(cell, start = nchar(cell), nchar(cell)))
})

cnames.keep.lst.lst <- lapply(jmarks, function(jmark){
  jlist <- split(x = dat.meta.lst[[jmark]], f = dat.meta.lst[[jmark]]$ctype)
  lapply(jlist, function(x) x$cell)
})

mats.pbulk <- lapply(jmarks, function(jmark){
  vec.lst <- SumAcrossClusters(mat.mark.lst[[jmark]], cnames.keep.lst = cnames.keep.lst.lst[[jmark]])
  dat <- as.data.frame(vec.lst)
})


plot(mats.pbulk$mark1$A, mats.pbulk$mark2$A)

plot(mats.pbulk$mark1$A, mats.pbulk$mark1$C)

plot(mats.pbulk$mark2$A, mats.pbulk$mark2$B)
plot(mats.pbulk$mark2$A, mats.pbulk$mark2$C)

plot(mats.pbulk$`mark1-mark2`$A, mats.pbulk$`mark1-mark2`$B)
plot(mats.pbulk$`mark1-mark2`$A, mats.pbulk$`mark1-mark2`$C)



plot(mats.pbulk$mark1$A + mats.pbulk$mark2$A, mats.pbulk$`mark1-mark2`$A)

plot(mats.pbulk$mark1$B + mats.pbulk$mark2$B, mats.pbulk$`mark1-mark2`$B)

plot(mats.pbulk$mark1$C + mats.pbulk$mark2$C, mats.pbulk$`mark1-mark2`$C)

plot(mats.pbulk$mark1$C + mats.pbulk$mark2$C, mats.pbulk$`mark1-mark2`$C, log = "xy")

