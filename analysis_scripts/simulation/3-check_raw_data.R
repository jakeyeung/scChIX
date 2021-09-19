# Jake Yeung
# Date of Creation: 2021-09-19
# File: ~/projects/scChIX/analysis_scripts/simulation/3-check_raw_data.R
#

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)





# Read raw data  ----------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

inf.sim <- file.path(hubprefix, "jyeung/data/dblchic/simulation_data/snakemake_inputs/countmats/ATAC_simulator_params_and_outputs.RData")
load(inf.sim, v=T)


jmarks <- c("mark1", "mark2", "mark1-mark2"); names(jmarks) <- jmarks
inf.mats <- lapply(jmarks, function(jmark){
  inf <- file.path(hubprefix, paste0("jyeung/data/dblchic/simulation_data/snakemake_inputs/countmats/countmat_var_filt.", jmark, ".rds"))
})

mats <- lapply(inf.mats, readRDS)


# Load meta  --------------------------------------------------------------

inf.meta.lst <- lapply(jmarks, function(jmark){
  inf.meta <- file.path(hubprefix, paste0("jyeung/data/dblchic/simulation_data/snakemake_outputs/objs_from_LDA/celltyping_output_filt.", jmark, ".rds"))
  return(inf.meta)
})

dat.meta.lst <- lapply(jmarks, function(jmark){
  jinf <- inf.meta.lst[[jmark]]
  dat.meta.tmp <- readRDS(jinf) %>%
    rowwise() %>%
    mutate(ctype = substr(cell, start = nchar(cell), nchar(cell)),
           mark = jmark)
})


cnames.keep.lst.lst <- lapply(jmarks, function(jmark){
  jlist <- split(x = dat.meta.lst[[jmark]], f = dat.meta.lst[[jmark]]$ctype)
  lapply(jlist, function(x) x$cell)
})


mats.pbulk <- lapply(jmarks, function(jmark){
  vec.lst <- SumAcrossClusters(mats[[jmark]], cnames.keep.lst = cnames.keep.lst.lst[[jmark]])
  dat <- as.data.frame(vec.lst)
})


plot(mats.pbulk$mark1$A, mats.pbulk$mark2$A)

plot(mats.pbulk$mark1$A, mats.pbulk$mark1$C)

plot(mats.pbulk$mark2$A, mats.pbulk$mark2$B)
plot(mats.pbulk$mark2$A, mats.pbulk$mark2$C)

plot(mats.pbulk$`mark1-mark2`$A, mats.pbulk$`mark1-mark2`$B)
plot(mats.pbulk$`mark1-mark2`$A, mats.pbulk$`mark1-mark2`$C)



# Add up two marks together, do they resemble?  ---------------------------


plot(mats.pbulk$mark1$A + mats.pbulk$mark2$A, mats.pbulk$`mark1-mark2`$A)

plot(mats.pbulk$mark1$B + mats.pbulk$mark2$B, mats.pbulk$`mark1-mark2`$B)

plot(mats.pbulk$mark1$C + mats.pbulk$mark2$C, mats.pbulk$`mark1-mark2`$C)

plot(mats.pbulk$mark1$C + mats.pbulk$mark2$C, mats.pbulk$`mark1-mark2`$C, log = "xy")


x.sum <- log10(mats.pbulk$mark1$C + mats.pbulk$mark2$C)
x.dbl <- log10(mats.pbulk$`mark1-mark2`$C)

# plot(log10(x.sum), log10(x.dbl))
plot(x.sum, x.dbl)
abline(a = 0.2, b = 1)

bins.check <- x.sum * 1 + 0.2 < x.dbl
col.vec <- ifelse(bins.check, "blue", "red")
plot(x.sum, x.dbl, col = col.vec)

length(which(bins.check)) / nrow(mats.pbulk$`mark1-mark2`)

bin.names <- rownames(mats.pbulk$`mark1-mark2`)[bins.check]

# Check discrepancy bins  -------------------------------------------------

bin.dat <- ctype.sim.counts.lst$C$bin.data %>%
  rowwise() %>%
  mutate(bin.check = Bin %in% bin.names)

ggplot(bin.dat, aes(x = BinMean, fill = bin.check)) +
  geom_density(alpha = 0.25) +
  scale_x_log10() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
