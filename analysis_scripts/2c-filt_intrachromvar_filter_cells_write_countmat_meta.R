# Jake Yeung
# Date of Creation: 2021-07-15
# File: ~/projects/scChIX/analysis_scripts/2c-filt_intrachromvar_write_countmat_featuresfilt.R
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

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")

# jquantile <- 0.2
jquantile <- 0.3

# Load  -------------------------------------------------------------------

jdate.output <- "2021-07-15"
jdate <- "2021-07-14"
# jmarks <- c("K36", "K27", "K36-K27")
jmarks <- c("K36", "K9m3", "K36-K9m3")
names(jmarks) <- jmarks
jstr <- paste(jmarks, collapse = "_")

hubprefix <- "/home/jyeung/hub_oudenaarden"
inmain <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_from_less_stringent2/", jmarks[[2]], "/filtered_count_tables_filtered_counts_varcutoffmin_0.03.noLSIfilt")
dir.exists(inmain)


# Load metas  -------------------------------------------------------------

infs.lda <- lapply(jmarks, function(jmark){
  inf.lda <- file.path(inmain, paste0("lda_output_filt.", jmark, ".2021-07-14.rds"))
  assertthat::assert_that(file.exists(inf.lda))
  return(inf.lda)
})

infs.meta <- lapply(jmarks, function(jmark){
  inf.meta <- file.path(inmain, paste0("celltyping_output_filt.", jmark, ".2021-07-14.rds"))
  assertthat::assert_that(file.exists(inf.meta))
  return(inf.meta)

})

infs.countmat <- lapply(jmarks, function(jmark){
  inf.countmat <- file.path(inmain, paste0("countmat_output_filt.", jmark, ".2021-07-14.rds"))
  assertthat::assert_that(file.exists(inf.countmat))
  return(inf.countmat)
})

ldas.lst <- lapply(infs.lda, readRDS)
dat.umap.lst <- lapply(infs.meta, readRDS)
countmats.lst <- lapply(infs.countmat, readRDS)



# Load LDAs  --------------------------------------------------------------

tm.result.lst <- lapply(ldas.lst, function(out.lda){
  tm.result <- posterior(out.lda)
})

dat.impute.log.lst <- lapply(tm.result.lst, function(tm.result){
  t(log2(tm.result$topics %*% tm.result$terms))
})

dat.var.lst <- lapply(dat.impute.log.lst, function(dat.impute.log){
  CalculateVarAll(dat.impute.log, jchromos)
})

dat.umap.merge.lst <- lapply(jmarks, function(jmark){
  dat.umap.merge <- left_join(dat.umap.lst[[jmark]], dat.var.lst[[jmark]])
})

m.lst <- lapply(jmarks, function(jmark){
  m <- ggplot(dat.umap.merge.lst[[jmark]], aes(x = umap1, y = umap2, color = log(cell.var.within.sum.norm))) +
    geom_point() +
    theme_bw() +
    scale_color_viridis_c() +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})



# ggplot(dat.umap.merge, aes(x = umap1, y = umap2, color = log(cell.var.within.sum.norm))) +
#   geom_point() +
#   theme_bw() +
#   scale_color_viridis_c() +
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

jcutoff.lst <- lapply(jmarks, function(jmark){
  jcutoff <- quantile(dat.umap.merge.lst[[jmark]]$cell.var.within.sum.norm, jquantile)
})


dat.umap.merge.annot.lst <- lapply(jmarks, function(jmark){
  dat.umap.merge.lst[[jmark]] %>%
    rowwise() %>%
    mutate(is.bad = cell.var.within.sum.norm <= jcutoff.lst[[jmark]])
})


m.cutoff.lst <- lapply(jmarks, function(jmark){
  m <- ggplot(dat.umap.merge.annot.lst[[jmark]], aes(x = log(cell.var.within.sum.norm))) +
    geom_density() +
    geom_vline(xintercept = log(jcutoff.lst[[jmark]])) +
    ggtitle(jmark) +
    theme_bw() +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})



m.cutoff.check.lst <- lapply(jmarks, function(jmark){
  m <- ggplot(dat.umap.merge.annot.lst[[jmark]], aes(x = umap1, y = umap2, color = is.bad)) +
    geom_point() +
    theme_bw() +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})


m.cutoff.check.plate.lst <- lapply(jmarks, function(jmark){
  m <- ggplot(dat.umap.merge.annot.lst[[jmark]], aes(x = umap1, y = umap2, color = is.bad)) +
    geom_point() +
    theme_bw() +
    facet_wrap(~plate) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})



# Write new counts, metas -------------------------------------------------

outmain <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_from_less_stringent2/var_filtered_", jquantile)
dir.create(outmain)
outdir <- file.path(outmain, jstr)
dir.create(outdir)

outpdf <- file.path(outdir, paste0("plots_check_intrachrom_var.", jstr, ".", jdate.output, ".pdf"))
pdf(outpdf, useDingbats = FALSE)
  print(m.lst)
  print(m.cutoff.lst)
  print(m.cutoff.check.lst)
  print(m.cutoff.check.plate.lst)
dev.off()

for (jmark in jmarks){
  print(jmark)
  # write dat impute vars
  dat.var.tmp <- dat.umap.merge.annot.lst[[jmark]]
  outf.var.tmp <- file.path(outdir, paste0("intrachrom_var_outputs.", jmark, ".", jdate.output, ".rds"))


  # write filtered meta data
  dat.meta.tmp <- subset(dat.var.tmp, !is.bad)
  outf.meta.tmp <- file.path(outdir, paste0("celltyping_var_filt.", jmark, ".", jdate.output, ".rds"))

  # write filtered count mat
  cells.keep.tmp <- dat.meta.tmp$cell
  dat.countmat.tmp <- countmats.lst[[jmark]][, cells.keep.tmp]
  outf.countmat.tmp <- file.path(outdir, paste0("countmat_var_filt.", jmark, ".", jdate.output, ".rds"))

  print(dim(countmats.lst[[jmark]]))
  print(dim(dat.countmat.tmp))

  saveRDS(dat.var.tmp, outf.var.tmp)
  saveRDS(dat.meta.tmp, outf.meta.tmp)
  saveRDS(dat.countmat.tmp, outf.countmat.tmp)
}

