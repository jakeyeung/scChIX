# Jake Yeung
# Date of Creation: 2021-06-29
# File: ~/projects/scChIX/analysis_scripts/1-explore_raw_count_data_K9me3_less_stringent_filter.R
# Use less stringnet var filter

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

# Load raw counts  --------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"


jmarks <- c("K36", "K9m3", "K36-K9m3")
names(jmarks) <- jmarks
jmark1 <- jmarks[[1]]
jmark2 <- jmarks[[2]]
jmarkdbl <- jmarks[[3]]

jsuffix <- "50000"
filtdir <- "filtered_counts_varcutoffmin_0.03"

# load mine
infs.lst <- lapply(jmarks, function(jmark) file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/count_tables/counts_tables_", jsuffix, "/", jmark, "/cbind_out/", filtdir, "/count_mat.", jmark, ".countcutoffmin_1000.TAcutoff_0.5.rds")))
lapply(infs.lst, function(x) file.exists(x))
mats.lst <- lapply(infs.lst, function(inf) readRDS(inf))

print(lapply(mats.lst, dim))

# Get common rows  --------------------------------------------------------

rows.lst <- lapply(mats.lst, function(jmat) rownames(jmat))
rows.common <- Reduce(intersect, x = rows.lst)


mats.lst.common <- lapply(mats.lst, function(jmat) jmat[rows.common, ])

mats.cbind.common <- do.call(cbind, mats.lst.common)

# Plot LSI all three together  --------------------------------------------

library(JFuncs)
library(irlba)

library(igraph)
library(umap)
library(hash)

jcheck <- scchicFuncs::RunLSI(count.mat = as.matrix(mats.cbind.common), n.components = 50)
jcheck2 <- scchicFuncs::RunLSI(count.mat = as.matrix(mats.lst$K36), n.components = 50)

lapply(jcheck, dim)
lapply(jcheck2, dim)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

dat.umap <- DoUmapAndLouvain(jcheck$u, jsettings = jsettings)

dat.umap2 <- DoUmapAndLouvain(jcheck2$u, jsettings = jsettings)

ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(dat.umap2, aes(x = umap1, y = umap2, color = louvain)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Assign stages -----------------------------------------------------------


# Filter same timepoints  -------------------------------------------------

# get timepoints?
source("/home/jyeung/projects/gastru_scchic/scripts/Rfunctions/QCFunctionsGastru.R")
dat.umap$stage <- sapply(dat.umap$cell, function(cell) StageToNumeric(GetStage(PreprocessSamp(cell))))
stages.keep <- unique(subset(dat.umap, mark == "K36-K9m3")$stage)
lapply(jmarks, function(jmark){
  unique(subset(dat.umap, mark == jmark)$stage)
})

# keep all stages



# Remove bad cells  -------------------------------------------------------

dat.umap$mark <- sapply(dat.umap$cell, function(x){
  if (grepl(jmarkdbl, x)){
    jmarktmp <- jmarkdbl
  } else if (grepl(jmark1, x)){
    jmarktmp <- jmark1
  } else if (grepl(jmark2, x)){
    jmarktmp <- jmark2
  }
  return(jmarktmp)
})

# dat.umap <- dat.umap %>%
#   rowwise()

umap1min <- 1
umap2max <- -1
# badclusts <- c("5", "8")
badclusts <- c("5")


ggplot(dat.umap, aes(x = umap1, y = umap2, color = mark)) +
  geom_point(alpha = 0.25) + theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap, aes(x = umap1, y = umap2, color = mark)) +
  geom_point() + theme_bw() + facet_wrap(~mark) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) +
  geom_point(alpha = 0.25) + theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap %>% mutate(louvain = louvain %in% badclusts), aes(x = umap1, y = umap2, color = louvain)) +
  geom_point(alpha = 0.25) + theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap %>% filter(mark == "K9m3") %>% mutate(louvain = louvain %in% badclusts), aes(x = umap1, y = umap2, color = louvain)) +
  geom_point(alpha = 0.25) + theme_bw() +
  facet_wrap(~stage) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap %>% filter(mark == "K36-K9m3") %>% mutate(louvain = louvain %in% badclusts), aes(x = umap1, y = umap2, color = louvain)) +
  geom_point(alpha = 0.25) + theme_bw() +
  facet_wrap(~stage) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap %>% filter(louvain != "3"), aes(x = umap1, y = umap2, color = mark)) +
  geom_point(alpha = 0.25) + theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap %>% filter(louvain != "3"), aes(x = umap1, y = umap2, color = mark)) +
  geom_point(alpha = 0.25) + theme_bw() + facet_wrap(~mark) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dat.umap.sum <- dat.umap %>%
  group_by(louvain, mark) %>%
  summarise(ncell = length(cell))

# find clustesr with lots of overlap
dat.umap.sum2 <- subset(dat.umap.sum, mark != jmarkdbl) %>%
  group_by(louvain) %>%
  mutate(nfrac = ncell / sum(ncell))



# Remove K9m3s that cluster with K36  -------------------------------------

# umap1min <- 0
# umap2max <- 1
# badclusts = c("")

dat.umap.midfilt <- subset(dat.umap, !louvain %in% badclusts)
dat.umap.midfilt2 <- dat.umap.midfilt %>% filter( ! ( grepl("K9m3", mark) & umap1 > umap1min & umap2 < umap2max) )

ggplot(dat.umap.midfilt2, aes(x = umap1, y = umap2, color = mark)) +
  geom_point(alpha = 0.25) + theme_bw() + facet_wrap(~mark) +
  coord_cartesian(xlim = c(-7, 7), ylim = c(-10, 10)) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# check k36me3 still OK ?

ggplot(dat.umap2 %>% filter(cell %in% dat.umap.midfilt2$cell), aes(x = umap1, y = umap2)) +
  geom_point(alpha = 0.25) + theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Write to output  --------------------------------------------------------

outdir <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis/filtered_count_tables_", filtdir), jmarks.str)
dir.create(outdir)

for (jmark in jmarks){
  print(jmark)
  outrds <- file.path(outdir, paste0("count_tables.", jsuffix, ".", jmark, ".", Sys.Date(), ".rds"))
  outmeta <- file.path(outdir, paste0("meta_data.", jsuffix, ".", jmark, ".", Sys.Date(), ".txt"))
  jsub <- subset(dat.umap.midfilt2, mark == jmark)
  cells.keep <- jsub$cell
  jmat.sub <- mats.lst.common[[jmark]][rows.common, cells.keep]
  print(dim(jmat.sub))
  saveRDS(jmat.sub, file = outrds)
  fwrite(jsub, file = outmeta)
}



