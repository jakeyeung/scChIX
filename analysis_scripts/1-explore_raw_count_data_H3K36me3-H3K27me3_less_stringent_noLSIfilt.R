# Jake Yeung
# Date of Creation: 2021-07-12
# File: ~/projects/scChIX/analysis_scripts/1-explore_raw_count_data_H3K36me3-H3K27me3_less_stringent_noLSIfilt.R
#


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

# Load raw counts  --------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"


jmarks <- c("K36", "K27", "K36-K27")
names(jmarks) <- jmarks
jmark1 <- jmarks[[1]]
jmark2 <- jmarks[[2]]
jmarkdbl <- jmarks[[3]]
jmarks.str <- paste(jmarks, collapse = "_")

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


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

# LSI each mat

dat.umap.lst <- lapply(mats.lst.common, function(jmat){
  jcheck <- scchicFuncs::RunLSI(count.mat = as.matrix(jmat), n.components = 50)
  dat.umap <- DoUmapAndLouvain(jcheck$u, jsettings = jsettings) %>%
    rowwise() %>%
    mutate(plate = ClipLast(cell, jsep = "_"),
           experi = ClipLast(plate, jsep = "-"),
           stage = strsplit(cell, split = "-")[[1]][[1]])
  return(dat.umap)
})

print(head(dat.umap.lst[[1]]))

outdir <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis/filtered_count_tables_", filtdir, ".noLSIfilt"), jmarks.str)
dir.create(outdir)

for (jmark in jmarks){
  print(jmark)
  outrds <- file.path(outdir, paste0("count_tables.", jsuffix, ".", jmark, ".", Sys.Date(), ".rds"))
  outmeta <- file.path(outdir, paste0("meta_data.", jsuffix, ".", jmark, ".", Sys.Date(), ".txt"))
  jsub <- dat.umap.lst[[jmark]]
  jmat.sub <- mats.lst.common[[jmark]][rows.common, ]
  print(dim(jmat.sub))
  saveRDS(jmat.sub, file = outrds)
  fwrite(jsub, file = outmeta)
}



