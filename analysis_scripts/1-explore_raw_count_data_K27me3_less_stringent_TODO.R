# Jake Yeung
# Date of Creation: 2021-06-30
# File: ~/projects/scChIX/analysis_scripts/1-explore_raw_count_data_K27me3_less_stringent_TODO.R
#


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

# Load raw counts  --------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"


# jmarks <- c("K36", "K9m3", "K36-K9m3")
jmarks <- c("K36", "K27", "K36-K27")
jmarks.str <- paste(jmarks, collapse = "_")
names(jmarks) <- jmarks
jmark1 <- jmarks[[1]]
jmark2 <- jmarks[[2]]
jmarkdbl <- jmarks[[3]]

jsuffix <- "50000"
filtdir <- "filtered_counts_varcutoffmin_0.03"

# load mine
infs.lst <- lapply(jmarks, function(jmark) file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/count_tables/counts_tables_", jsuffix, "/", jmark, "/cbind_out/", filtdir, "/count_mat.", jmark, ".countcutoffmin_1000.TAcutoff_0.5.rds")))

assertthat::assert_that(all(sapply(infs.lst, function(x) file.exists(x))))

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
jcheck2 <- scchicFuncs::RunLSI(count.mat = as.matrix(mats.lst[[jmark1]]), n.components = 50)

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

# Filter same timepoints  -------------------------------------------------


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


# get timepoints?
source("/home/jyeung/projects/gastru_scchic/scripts/Rfunctions/QCFunctionsGastru.R")
dat.umap$stage <- sapply(dat.umap$cell, function(cell) StageToNumeric(GetStage(PreprocessSamp(cell))))
stages.keep <- unique(subset(dat.umap, mark == "K36-K9m3")$stage)
lapply(jmarks, function(jmark){
  unique(subset(dat.umap, mark == jmark)$stage)
})

# keep all stages


# Remove bad cells  -------------------------------------------------------


ggplot(dat.umap, aes(x = umap1, y = umap2, color = mark)) +
  geom_point(alpha = 0.25) + theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap, aes(x = umap1, y = umap2, color = mark)) +
  geom_point() + theme_bw() + facet_wrap(~mark) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) +
  geom_point(alpha = 0.25) + theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) +
  geom_point(alpha = 0.8) + theme_bw() + facet_wrap(~mark) +
  scale_color_manual(values = cbPalette) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# badclust <- "5"
badclusts <- c("8")
ggplot(dat.umap %>% filter(!louvain %in% badclusts), aes(x = umap1, y = umap2, color = mark)) +
  geom_point(alpha = 0.25) + theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap %>% filter(!louvain %in% badclusts), aes(x = umap1, y = umap2, color = mark)) +
  geom_point(alpha = 0.25) + theme_bw() + facet_wrap(~mark) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap %>% filter(!louvain %in% badclusts), aes(x = umap1, y = umap2, color = louvain)) +
  geom_point(alpha = 0.85) + theme_bw() + facet_wrap(~mark) +
  scale_color_manual(values = cbPalette) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Remove K27m3s that cluster with K36  -------------------------------------

umap1max <- 1
umap2max <- 0

ggplot(dat.umap %>% filter(!louvain %in% badclusts), aes(x = umap1, y = umap2, color = mark)) +
  geom_point(alpha = 0.25) + theme_bw() + facet_wrap(~mark) +
  geom_vline(xintercept = umap1max) + geom_hline(yintercept = umap2min) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dat.umap.midfilt <- subset(dat.umap, !louvain %in% badclusts)
dat.umap.midfilt2 <- dat.umap.midfilt %>% filter( ! ( mark == jmarkdbl & umap1 < umap1max & umap2 < umap2max) )

ggplot(dat.umap.midfilt2, aes(x = umap1, y = umap2, color = mark)) +
  geom_point(alpha = 0.25) + theme_bw() + facet_wrap(~mark) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# check k36me3 still OK ?

ggplot(dat.umap2 %>% filter(cell %in% dat.umap.midfilt2$cell), aes(x = umap1, y = umap2)) +
  geom_point(alpha = 0.25) + theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Filter stages that are matching  ----------------------------------------



# Write to output  --------------------------------------------------------

outdir <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis/filtered_count_tables_", filtdir), jmarks.str)
dir.create(outdir)

for (jmark in jmarks){
  print(jmark)
  outrds <- file.path(outdir, paste0(jmarks.str, ".count_tables.", jsuffix, ".", jmark, ".", Sys.Date(), ".rds"))
  outmeta <- file.path(outdir, paste0(jmarks.str, ".meta_data.", jsuffix, ".", jmark, ".", Sys.Date(), ".txt"))
  if (file.exists(outrds)){
    print(paste("file exists", outrds))
    next
  }
  if (file.exists(outmeta)){
    print(paste("file exists", outmeta))
    next
  }
  jsub <- subset(dat.umap.midfilt2, mark == jmark)
  cells.keep <- jsub$cell
  jmat.sub <- mats.lst.common[[jmark]][rows.common, cells.keep]
  print(dim(jmat.sub))
  saveRDS(jmat.sub, file = outrds)
  fwrite(jsub, file = outmeta)
}



