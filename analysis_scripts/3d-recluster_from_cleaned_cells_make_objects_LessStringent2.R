# Jake Yeung
# Date of Creation: 2021-07-14
# File: ~/projects/scChIX/analysis_scripts/3d-recluster_from_cleaned_cells_make_objects_LessStringent_BadPlatesFilt.R
#


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(scChIX)

library(topicmodels)
library(hash)
library(igraph)
library(umap)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123


# Load new LDA after cleaning  --------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"


jsuffix <- "filtered_count_tables_filtered_counts_varcutoffmin_0.03.noLSIfilt"
# jdate <- "2021-07-12"
# jmarks <- c("K36", "K9m3", "K36-K9m3")
jdate <- "2021-07-12"
jmarks <- c("K36", "K27", "K36-K27")
names(jmarks) <- jmarks
jstr <- paste(jmarks, collapse = "_")

jmain <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_from_less_stringent2"))
# assertthat::assert_that(!dir.exists(jmain))
dir.create(jmain, recursive = FALSE)
outmain <- file.path(jmain, jmarks[[2]])
assertthat::assert_that(!dir.exists(outmain))
dir.create(outmain, recursive = FALSE)

jname <- "count_tables.50000"

for (jmark in jmarks){
  infrobj.check <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/LDA_outputs/ldaAnalysis_50000_", jsuffix, "_", jstr, "/lda_outputs.", jname, ".", jmark, ".", jdate, ".K-30.binarize.FALSE/ldaOut.", jname, ".", jmark, ".", jdate, ".K-30.Robj"))
  assertthat::assert_that(file.exists(infrobj.check))
  print(infrobj.check)
}



for (jmark in jmarks){
  infrobj <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/LDA_outputs/ldaAnalysis_50000_", jsuffix, "_", jstr, "/lda_outputs.", jname, ".", jmark, ".", jdate, ".K-30.binarize.FALSE/ldaOut.", jname, ".", jmark, ".", jdate, ".K-30.Robj"))
  load(infrobj, v=T)

  tm.result <- posterior(out.lda)

  dat.umap.long <- DoUmapAndLouvain(tm.result$topics, jsettings = jsettings)

  # annotate
  dat.umap.long <- dat.umap.long %>%
    rowwise() %>%
    mutate(plate = ClipLast(cell, jsep = "_"),
           experi = ClipLast(plate, jsep = "-"),
           stage =  strsplit(cell, "-")[[1]][[1]],
           cluster = paste("cluster", louvain, sep = "")) %>%
    dplyr::select(-louvain)


  # Write tables  -----------------------------------------------------------

  # jdate <- Sys.Date()
  # write count tables, metadata, and LDA objects
  outdir <- file.path(outmain, jsuffix)
  dir.create(outdir)

  m <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = cluster)) +
    geom_point() +
    theme_bw() +
    ggtitle(paste(jsuffix, jmark)) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  m1 <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = cluster)) +
    geom_point() +
    theme_bw() +
    facet_wrap(~plate) +
    ggtitle(paste(jsuffix, jmark)) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  m2 <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = cluster)) +
    geom_point() +
    theme_bw() +
    facet_wrap(~experi) +
    ggtitle(paste(jsuffix, jmark)) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  m3 <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = cluster)) +
    geom_point() +
    theme_bw() +
    facet_wrap(~stage) +
    ggtitle(paste(jsuffix, jmark)) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  outmat <- file.path(outdir, paste0("countmat_output_filt.", jmark, ".", Sys.Date(), ".rds"))
  outmeta <- file.path(outdir, paste0("celltyping_output_filt.", jmark, ".", Sys.Date(), ".rds"))
  outlda <- file.path(outdir, paste0("lda_output_filt.", jmark, ".", Sys.Date(), ".rds"))
  outpdf <- file.path(outdir, paste0("plots.", jmark, ".", Sys.Date(), ".pdf"))

  saveRDS(count.mat, file = outmat)
  saveRDS(dat.umap.long, file = outmeta)
  saveRDS(out.lda, file = outlda)

  pdf(outpdf, useDingbats = FALSE)
    print(m)
    print(m1)
    print(m2)
    print(m3)
  dev.off()
}

