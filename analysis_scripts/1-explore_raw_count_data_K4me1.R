# Jake Yeung
# Date of Creation: 2021-06-29
# File: ~/projects/scChIX/analysis_scripts/1-explore_raw_count_data_K4me1.R
#

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

hubprefix <- "/home/jyeung/hub_oudenaarden"

# # check
# inf1 <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/count_tables/counts_tables_50000/K36-K4m1/E9p5-CB6-K36-K4m1-190409-1.sorted.tagged.countTable.binsize_50000.csv"
# dat1 <- ReadMatSlideWinFormat(inf1)
rds1 <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/count_tables/counts_tables_50000/K36-K4m1/cbind_out/K36-K4m1.countTable.binsize_50000.rds")
datrds1 <- readRDS(rds1)

# Load raw counts  --------------------------------------------------------


# jmarks <- c("K36", "K9m3", "K36-K9m3")
jmarks <- c("K36", "K4m1", "K36-K4m1", "K9m3")
jmarks.str <- paste(jmarks[1:3], collapse = "_")
names(jmarks) <- jmarks
jmark1 <- jmarks[[1]]
jmark2 <- jmarks[[2]]
jmarkdbl <- jmarks[[3]]
jmarkneg <- jmarks[[4]]

jmarks.sub <- jmarks[jmarks != jmarkneg]

jsuffix <- "50000"

# load mine
infs.lst <- lapply(jmarks, function(jmark) file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/count_tables/counts_tables_", jsuffix, "/", jmark, "/cbind_out/filtered_counts/count_mat.", jmark, ".countcutoffmin_1000.TAcutoff_0.5.rds")))

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


# Remove bad cells  -------------------------------------------------------

dat.umap$mark <- sapply(dat.umap$cell, function(x){
  if (grepl(jmarkdbl, x)){
    jmarktmp <- jmarkdbl
  } else if (grepl(jmark1, x)){
    jmarktmp <- jmark1
  } else if (grepl(jmark2, x)){
    jmarktmp <- jmark2
  } else if (grepl(jmarkneg, x)){
    jmarktmp <- jmarkneg
  }
  return(jmarktmp)
})

ggplot(dat.umap, aes(x = umap1, y = umap2, color = mark)) +
  geom_point(alpha = 0.25) + theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap, aes(x = umap1, y = umap2, color = mark)) +
  geom_point() + theme_bw() + facet_wrap(~mark) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) +
  geom_point(alpha = 0.9) + theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

badclust <- "10"
ggplot(dat.umap %>% filter(louvain != badclust), aes(x = umap1, y = umap2, color = mark)) +
  geom_point(alpha = 0.25) + theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap %>% filter(louvain != badclust), aes(x = umap1, y = umap2, color = mark)) +
  geom_point(alpha = 0.25) + theme_bw() + facet_wrap(~mark) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Filter good cells for K36 -----------------------------------------------

inf.k36.meta <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis/filtered_count_tables/K36_K9_K36-K9/meta_data.50000.K36.2021-06-28.txt")
dat.k36.meta <- fread(inf.k36.meta)
cells.k36 <- dat.k36.meta$cell

# fillter out some marks
dat.umap.k36 <- subset(dat.umap, cell %in% cells.k36 & mark != jmarkneg)
dat.umap.k4m1 <- subset(dat.umap, mark == jmark2)

# get timepoints?
source("/home/jyeung/projects/gastru_scchic/scripts/Rfunctions/QCFunctionsGastru.R")
dat.umap$stage <- sapply(dat.umap$cell, function(cell) StageToNumeric(GetStage(PreprocessSamp(cell))))
umap2min <- 1
stages.keep <- unique(subset(dat.umap, mark == jmarkdbl)$stage)
dat.umap.filt <- subset(dat.umap, stage %in% stages.keep & louvain != badclust)
dat.umap.k36 <- subset(dat.umap.filt, cell %in% cells.k36 & louvain != badclust)
dat.umap.k4m1 <- subset(dat.umap.filt, mark == jmark2 & louvain != badclust & umap2 > umap2min)
dat.umap.dbl <- subset(dat.umap.filt, mark == jmarkdbl & umap2 > umap2min)

dat.umap.midfilt2 <- bind_rows(list(dat.umap.k36, dat.umap.k4m1, dat.umap.dbl))


# Write to output  --------------------------------------------------------

outdir <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis/filtered_count_tables"), jmarks.str)
dir.create(outdir)

for (jmark in jmarks.sub){
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



