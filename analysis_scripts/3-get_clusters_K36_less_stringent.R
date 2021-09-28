# Jake Yeung
# Date of Creation: 2021-07-12
# File: ~/projects/scChIX/analysis_scripts/3-get_clusters_K36_less_stringent.R
# description


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)

library(scchicFuncs)

library(hash)
library(igraph)
library(umap)

library(mixtools)

library(ggrepel)
source("/home/jyeung/projects/gastru_scchic/scripts/Rfunctions/QCFunctionsGastru.R")

nn <- 30
jsettings <- umap.defaults
jsettings$n_neighbors <- nn
jsettings$min_dist <- 0.2
jsettings$random_state <- 123
jsettings$spread <- 5

outdir <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_less_stringent_dbl_cleaned")
dir.create(outdir)
# outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/celltyping_FromTopics"

# jsettings2 <- jsettings
# jsettings2$n_neighbors <- 15
# jsettings2$min_dist <- 0.2
# jsettings2$spread <- 5

# Load  -------------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

jsuffix <- "50000"

jmarks <- c("K36", "K9m3", "K36-K9m3")
names(jmarks) <- jmarks

# jmark <- jmarks[[1]]

jdate <- "2021-06-29"
infs <- lapply(jmarks, function(jmark){
  print(jmark)
  inf <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/LDA_outputs/ldaAnalysis_", jsuffix, "/lda_outputs.count_tables.", jsuffix, ".", jmark, ".", jdate, ".K-30.binarize.FALSE/ldaOut.count_tables.", jsuffix, ".", jmark, ".", jdate, ".K-30.Robj"))
  assertthat::assert_that(file.exists(inf))
  return(inf)
})

out.lst <- lapply(infs, function(inf){
  load(inf, v=T)
  return(list(out.lda = out.lda, count.mat = count.mat))
})


count.mat.lst <- lapply(out.lst, function(jout){
  # load(inf, v=T)  # out.lda, count.mat
  # return(count.mat)
  # return(jout$out.lda)
  return(jout$count.mat)
})

out.lda.lst <- lapply(out.lst, function(jout){
  # load(inf, v=T)  # out.lda, count.mat
  # return(out.lda)
  return(jout$out.lda)
})

tm.result.lst <- lapply(jmarks, function(jmark){
  # load(inf, v=T)  # out.lda
  out.lda <- out.lda.lst[[jmark]]
  tm.result <- posterior(out.lda)
  tm.result <- AddTopicToTmResult(tm.result)
  return(tm.result)
})

dat.umap.lst <- lapply(tm.result.lst, function(tm.result){
  dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings = jsettings) %>%
    rowwise() %>%
    mutate(plate = ClipLast(cell, jsep = "_"),
           experi = ClipLast(plate, jsep = "-"))
  return(dat.umap)
})



cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

m.lst <- lapply(jmarks, function(jmark){
  jdat <- dat.umap.lst[[jmark]] %>%
    rowwise() %>%
    mutate(stage = strsplit(cell, split = "-")[[1]][[1]])
  m <- ggplot(jdat, aes(x = umap1, y = umap2, color = louvain)) +
    geom_point() +
    theme_bw() +
    facet_wrap(~stage) +
    ggtitle(jmark) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})
print(m.lst[[2]])


inf.cells.stringent <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis/filtered_count_tables/K36_K9_K36-K9/meta_data.50000.K9m3.2021-06-28.txt"))
dat.cells.stringent <- fread(inf.cells.stringent)

cells.stringent <- dat.cells.stringent$cell

jdat.check <- dat.umap.lst$K9m3 %>%
  rowwise() %>%
  mutate(stage = strsplit(cell, split = "-")[[1]][[1]],
         in.stringent = cell %in% cells.stringent)

m <- ggplot(jdat.check, aes(x = umap1, y = umap2, color = in.stringent)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~stage) +
  ggtitle("K9m3") +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

m <- ggplot(jdat.check, aes(x = umap1, y = umap2, color = louvain)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~stage) +
  ggtitle("K9m3") +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

m <- ggplot(jdat.check %>% filter(stage == "E9p5"), aes(x = umap1, y = umap2, color = louvain)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~plate) +
  ggtitle("K9m3") +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


table(subset(jdat.check, louvain == 8)$in.stringent)

for (jmarktmp in jmarks){
  tm.result <- tm.result.lst[[jmarktmp]]
  # dat.umap.long <- dat.umap.lst[[jmarktmp]]

  # add topic info
  dat.topics <- data.frame(cell = rownames(tm.result$topics), tm.result$topics, stringsAsFactors = FALSE)
  dat.umap.long <- left_join(dat.umap.lst[[jmarktmp]], dat.topics)

  pdf(file.path(outdir, paste0("celltyping_alltopics_from_", jmarktmp, ".", Sys.Date(), ".pdf")))


  # Double check by louvain clustering  -------------------------------------

  m.louv <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = as.character(louvain))) +
    geom_point() +
    scale_color_manual(values = cbPalette) +
    ggtitle(jmarktmp) +
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m.louv)

  m.louv.experi <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = as.character(louvain))) +
    geom_point() +
    scale_color_manual(values = cbPalette) +
    ggtitle(jmarktmp) +
    facet_wrap(~experi) +
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  m.louv.plate <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = as.character(louvain))) +
    geom_point() +
    scale_color_manual(values = cbPalette) +
    ggtitle(jmarktmp) +
    facet_wrap(~plate) +
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  print(m.louv)
  print(m.louv.plate)
  print(m.louv.experi)

  dev.off()



  # Write outputs:   ----------------------------------------------------------

  dat.merge <- dat.umap.long %>%
    dplyr::rename(cluster = louvain)

  outrds.meta <- file.path(outdir, paste0("celltyping_output_LessStringent.", jmarktmp, ".", Sys.Date(), ".rds"))
  saveRDS(dat.merge, file = outrds.meta)

  outrds.countmat <- file.path(outdir, paste0("countmat_output_LessStringent.", jmarktmp, ".", Sys.Date(), ".rds"))
  saveRDS(count.mat.lst[[jmarktmp]], file = outrds.countmat)

  outrds.outlda <- file.path(outdir, paste0("lda_output_LessStringent.", jmarktmp, ".", Sys.Date(), ".rds"))
  saveRDS(out.lda.lst[[jmarktmp]], file = outrds.outlda)

  # outrdata <- file.path(outdir, paste0("celltyping_output.", jmarktmp, ".", Sys.Date(), ".RData"))

  # save RData with dat.merge and out.lda


}




