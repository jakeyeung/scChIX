# Jake Yeung
# Date of Creation: 2021-07-13
# File: ~/projects/scChIX/analysis_scripts/3d-recluster_from_cleaned_cells_make_objects_LessStringent.R
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


jdate <- "2021-07-12"
jsuffix <- "filtered_count_tables_filtered_counts_varcutoffmin_0.03.noLSIfilt"
jmarks <- c("K36", "K9m3", "K36-K9m3")
# jmarks <- c("K36", "K27", "K36-K27")
names(jmarks) <- jmarks

outmain <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_from_less_stringent.", jmarks[[3]])
dir.create(outmain)

jname <- "count_tables.50000"

for (jmark in jmarks){
  infrobj.check <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/LDA_outputs/ldaAnalysis_50000_", jsuffix, "/lda_outputs.", jname, ".", jmark, ".", jdate, ".K-30.binarize.FALSE/ldaOut.", jname, ".", jmark, ".", jdate, ".K-30.Robj"))
  assertthat::assert_that(file.exists(infrobj.check))
}


for (jmark in jmarks){
  infrobj <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/LDA_outputs/ldaAnalysis_50000_", jsuffix, "/lda_outputs.", jname, ".", jmark, ".", jdate, ".K-30.binarize.FALSE/ldaOut.", jname, ".", jmark, ".", jdate, ".K-30.Robj"))
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

  outmat <- file.path(outdir, paste0("countmat_output_filt.", jmark, ".", Sys.Date(), ".rds"))
  outmeta <- file.path(outdir, paste0("celltyping_output_filt.", jmark, ".", Sys.Date(), ".rds"))
  outlda <- file.path(outdir, paste0("lda_output_filt.", jmark, ".", Sys.Date(), ".rds"))
  outpdf <- file.path(outdir, paste0("plots.", jmark, ".", Sys.Date(), ".pdf"))

  saveRDS(count.mat, file = outmat)
  saveRDS(dat.umap.long, file = outmeta)
  saveRDS(out.lda, file = outlda)

  pdf(outpdf, useDingbats = FALSE)
    print(m)
  dev.off()
}


# # jcheck <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_from_less_stringent/filtered_count_tables_filtered_counts_varcutoffmin_0.03.noLSIfilt/celltyping_output_filt.K36-K9m3.2021-07-13.rds")
# jcheck <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_from_less_stringent/filtered_count_tables_filtered_counts_varcutoffmin_0.03.noLSIfilt/celltyping_output_filt.K36.2021-07-13.rds")
# dat.check <- readRDS(jcheck)
#
# ggplot(dat.check, aes(x = umap1, y = umap2, color = cluster)) +
#   geom_point() +
#   facet_wrap(~plate) +
#   theme_bw() +
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#

