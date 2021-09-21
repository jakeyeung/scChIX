# Jake Yeung
# Date of Creation: 2021-09-19
# File: ~/projects/scChIX/analysis_scripts/simulation/2-scchix_downstream_checks.R
#


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(dplyr)

library(topicmodels)


# Functions ---------------------------------------------------------------



# Constants ---------------------------------------------------------------


hubprefix <- "/home/jyeung/hub_oudenaarden"

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
jmarks <- c("mark1", "mark2", "mark1-mark2"); names(jmarks) <- jmarks

# frac.mutexcl.str <- "0.5"
# frac.mutexcl.str <- "1.0"
frac.mutexcl.str.vec <- c("0.01", "0.5", "1.0")

for (frac.mutexcl.str in frac.mutexcl.str.vec){


  print(frac.mutexcl.str)


  dmain <- file.path(hubprefix, paste0("jyeung/data/dblchic/simulation_data/frac_mutexcl_", frac.mutexcl.str))
  assertthat::assert_that(dir.exists(dmain))

  outdir <- file.path(hubprefix, "jyeung/data/dblchic/from_cluster/simulation_data_downstream/plots")
  assertthat::assert_that(dir.exists(outdir))

  outpdf <- file.path(outdir, paste0("plot_downstream_overlaps.frac_mutexcl_", frac.mutexcl.str, ".pdf"))
  outrds1 <- file.path(outdir, paste0("dat_umap_joined.downstream_overlaps.frac_mutexcl_", frac.mutexcl.str, ".rds"))
  outrds2 <- file.path(outdir, paste0("bin_means.downstream_overlaps.frac_mutexcl_", frac.mutexcl.str, ".rds"))
  outrds3 <- file.path(outdir, paste0("LL_mat_avg_byctype.frac_mutexcl_", frac.mutexcl.str, ".rds"))

  pdf(file = outpdf, useDingbats = FALSE)

  # Load sim params ---------------------------------------------------------

  inf.sim <- file.path(dmain, "snakemake_inputs/countmats/ATAC_simulator_params_and_outputs.RData")
  load(inf.sim, v=T)


  # Load unmixing mat -------------------------------------------------------

  inf.prob <- file.path(dmain, "snakemake_outputs/scchix_unmixing_downstream/scchix_inputs_clstr_by_celltype-prob_mat.mark1-mark2_to_mark1.txt")
  dat.prob <- fread(inf.prob)
  mat.prob <- as.matrix(subset(dat.prob, select = -V1))
  rownames(mat.prob) <- dat.prob$V1

  plot(density(mat.prob[, 1]))



  # Load unmixing  ----------------------------------------------------------



  # Check outputs -----------------------------------------------------------

  inf.input <- file.path(dmain, "snakemake_outputs/scchix_inputs_objs/scchix_inputs_clstr_by_celltype_mark1-mark2.removeNA_FALSE.RData")
  load(inf.input, v=T)

  inf.output <- file.path(dmain, "snakemake_outputs/scchix_outputs_objs/scchix_inputs_clstr_by_celltype_mark1-mark2.RData")
  load(inf.output, v=T)

  fits.out <- act.repress.coord.lst
  w.lst <- sapply(fits.out, function(x){
    return(x$w)
  })
  plot(hist(w.lst))  # more active than repressive? Why is that?

  plot(density(w.lst))

  out.dat <- data.frame(cell = names(w.lst), w = w.lst, stringsAsFactors = FALSE)

  # bimodal???
  ggplot(out.dat, aes(x = w)) + geom_histogram() +
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


  # get active, repress index for each cell
  cell.vec <- names(fits.out)
  names(cell.vec) <- cell.vec
  coords.dbl <- lapply(cell.vec, function(jcell){
    jfit <- fits.out[[jcell]]
    jweight <- fits.out[[jcell]]$w
    p.mat <- SoftMax(jfit$ll.mat)
    jcoord <- which(jfit$ll.mat == max(jfit$ll.mat), arr.ind = TRUE)
    jmax <- max(p.mat)
    mark1.clstname <- rownames(jfit$ll.mat)[[jcoord[[1]]]]
    mark2.clstname <- colnames(jfit$ll.mat)[[jcoord[[2]]]]
    out.dat <- data.frame(cell = jcell, louv.act = mark1.clstname, louv.repress = mark2.clstname, lnprob = jmax, w = jweight, stringsAsFactors = FALSE)
    return(out.dat)
  }) %>%
    bind_rows()

  library(ggforce)
  ggplot(coords.dbl, aes(x = as.character(louv.act), y = as.character(louv.repress))) +
    geom_point(position = position_jitternormal(sd_x = 0.1, sd_y = 0.1), alpha = 0.25) +
    theme_bw() +
    scale_color_manual(values = cbPalette, na.value = "grey85") +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    xlab("Active Clusters") + ylab("Repressed Clusters") +
    ggtitle("Each dot is a doubel stained cell,\nX-Y shows the cluster pair it is assigned")



  # Plot UMAP  --------------------------------------------------------------

  inf.meta.dbl <- file.path(dmain, "snakemake_outputs/objs_from_LDA/celltyping_output_filt.mark1-mark2.rds")
  dat.meta.dbl <- readRDS(inf.meta.dbl) %>%
    rowwise() %>%
    mutate(ctype = strsplit(cell, "_")[[1]][[3]])


  ggplot(dat.meta.dbl, aes(x = umap1, y = umap2, color = cluster)) +
    geom_point() +
    facet_wrap(~ctype) +
    theme_bw() +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  coords.dbl.annot <- left_join(coords.dbl, subset(dat.meta.dbl, select = c(cell, umap1, umap2, ctype)))

  ggplot(coords.dbl.annot, aes(x = as.character(louv.act), y = as.character(louv.repress), color = ctype)) +
    geom_point(position = position_jitternormal(sd_x = 0.1, sd_y = 0.1), alpha = 0.25) +
    theme_bw() +
    scale_color_manual(values = cbPalette, na.value = "grey85") +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    xlab("Active Clusters") + ylab("Repressed Clusters") +
    ggtitle("Each dot is a doubel stained cell,\nX-Y shows the cluster pair it is assigned")

  # plot likelihood matrix

  coords.dbl.annot.byctype <- coords.dbl.annot %>%
    ungroup() %>%
    group_by(ctype, louv.act, louv.repress) %>%
    summarise(lnprob = mean(lnprob),
              prob = exp(lnprob))
  exp(coords.dbl$lnprob)




  # Check other  ------------------------------------------------------------


  inf.meta.lst <- lapply(jmarks, function(jmark){
    inf.meta <- file.path(dmain, paste0("snakemake_outputs/objs_from_LDA/celltyping_output_filt.", jmark, ".rds"))
    return(inf.meta)
  })

  dat.meta.lst <- lapply(jmarks, function(jmark){
    jinf <- inf.meta.lst[[jmark]]
    dat.meta.tmp <- readRDS(jinf) %>%
      rowwise() %>%
      mutate(ctype = substr(cell, start = nchar(cell), nchar(cell)),
             mark = jmark)
  })

  m.meta.ctype <- lapply(jmarks, function(jmark){
    jdat <- dat.meta.lst[[jmark]]
    ggplot(jdat, aes(x = umap1, y = umap2, color = ctype)) +
      geom_point() +
      theme_bw() +
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  })

  m.meta.cluster <- lapply(jmarks, function(jmark){
    jdat <- dat.meta.lst[[jmark]]
    ggplot(jdat, aes(x = umap1, y = umap2, color = cluster)) +
      geom_point() +
      theme_bw() +
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  })

  JFuncs::multiplot(m.meta.ctype[[1]], m.meta.ctype[[2]], m.meta.ctype[[3]], cols = 3)
  JFuncs::multiplot(m.meta.cluster[[1]], m.meta.cluster[[2]], m.meta.cluster[[3]], cols = 3)


  # Connect the UMAPs  ------------------------------------------------------


  inf1 <- file.path(dmain, "snakemake_outputs/projection_output.mark1.RData")
  inf2 <- file.path(dmain, "snakemake_outputs/projection_output.mark2.RData")

  load(inf1, v=T)
  out.lda.predict1 <- out.lda.predict
  count.mat.proj1 <- count.mat.proj
  out.objs1 <- out.objs
  load(inf2, v=T)
  out.lda.predict2 <- out.lda.predict
  count.mat.proj2 <- count.mat.proj
  out.objs2 <- out.objs


  # Plot outputs ------------------------------------------------------------

  jsettings1 <- umap.defaults
  jsettings1$n_neighbors <- 30
  jsettings1$min_dist <- 0.1
  jsettings1$random_state <- 123
  jsettings1$spread <- 1

  jsettings2 <- umap.defaults
  jsettings2$n_neighbors <- 30
  jsettings2$min_dist <- 0.1
  jsettings2$random_state <- 123
  jsettings2$spread <- 1

  umap.out1 <- umap(posterior(out.objs1$out.lda)$topics, config = jsettings1)
  umap.out2 <- umap(posterior(out.objs2$out.lda)$topics, config = jsettings2)


  umap.out1.pred <- predict(umap.out1, out.lda.predict1$topics)
  umap.out2.pred <- predict(umap.out2, out.lda.predict2$topics)

  dat.umap.out1.pred <- data.frame(cell = rownames(umap.out1.pred), umap1 = umap.out1.pred[, 1], umap2 = umap.out1.pred[, 2], experi = "double", stringsAsFactors = FALSE)
  dat.umap.out2.pred <- data.frame(cell = rownames(umap.out2.pred), umap1 = umap.out2.pred[, 1], umap2 = umap.out2.pred[, 2], experi = "double", stringsAsFactors = FALSE)

  dat.umap.out1 <- data.frame(cell = rownames(umap.out1$layout), umap1 = umap.out1$layout[, 1], umap2 = umap.out1$layout[, 2], experi = "single")
  dat.umap.out2 <- data.frame(cell = rownames(umap.out2$layout), umap1 = umap.out2$layout[, 1], umap2 = umap.out2$layout[, 2], experi = "single")

  dat.umap.out1.merged <- rbind(dat.umap.out1, dat.umap.out1.pred) %>%
    mutate(mark = "mark1")
  dat.umap.out2.merged <- rbind(dat.umap.out2, dat.umap.out2.pred) %>%
    mutate(mark = "mark2")

  ggplot(dat.umap.out1.merged, aes(x = umap1, y = umap2)) +
    geom_point() +
    theme_bw() +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  ggplot(dat.umap.out2.merged, aes(x = umap1, y = umap2)) +
    geom_point() +
    theme_bw() +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


  # Join umaps --------------------------------------------------------------

  dat.umap.merged <- rbind(dat.umap.out1.merged, dat.umap.out2.merged)

  dat.final.annots <- dat.meta.lst %>%
    bind_rows()

  dat.final.annots.sub <- subset(dat.final.annots, select = c(cell, cluster, ctype))

  dat.umap.joined <- left_join(dat.umap.merged, dat.final.annots.sub) %>%
    group_by(mark) %>%
    mutate(umap1.scale = scale(umap1, center = TRUE, scale = TRUE),
           umap2.scale = scale(umap2, center = TRUE, scale = TRUE),
           umap1.scale.shift = ifelse(mark == "mark1", umap1.scale - 3, umap1.scale + 3))


  ggplot(dat.umap.joined, aes(x = umap1, y = umap2, color = ctype)) +
    geom_point() +
    facet_wrap(~mark, nrow = 1) +
    theme_bw() +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")


  ggplot(dat.umap.joined, aes(x = umap1.scale.shift, y = umap2.scale, color = ctype, group = cell)) +
    geom_point() +
    geom_path(alpha = 0.02) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    scale_color_manual(values = cbPalette, na.value = "grey85") +
    theme_bw() +
    theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

  ggplot(dat.umap.joined, aes(x = umap1.scale.shift, y = umap2.scale, color = ctype, group = cell)) +
    geom_point() +
    geom_path(alpha = 0.02) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    scale_color_manual(values = cbPalette, na.value = "grey85") +
    facet_wrap(~cluster) +
    theme_bw() +
    theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")




  # Verify we can infer mutual and overlapping bins  ------------------------

  inf.mats <- lapply(jmarks, function(jmark){
    inf <- file.path(dmain, paste0("snakemake_inputs/countmats/countmat_var_filt.", jmark, ".rds"))
    assertthat::assert_that(file.exists(inf))
    return(inf)
  })
  mats <- lapply(inf.mats, readRDS)



  jctypes <- c("A", "B", "C"); names(jctypes) <- jctypes
  jctype <- "A"

  bin.means.merged.lst <- lapply(jctypes, function(jctype){
    grepsuffix <- paste0("_", jctype, "$")
    bin.means.merged <- MergeMarksEstimateOverlap(bin.dat = ctype.sim.counts.lst[[jctype]]$bin.data)
    bin.means.merged$ctype <- jctype
    return(bin.means.merged)
  })

  # add in real estimates
  cell.annot.vec <- sapply(colnames(mat.prob), function(x){
    substr(x, start = nchar(x), stop = nchar(x))
  })
  cell.annot <- data.frame(ctype = cell.annot.vec, cell = names(cell.annot.vec), stringsAsFactors = FALSE)
  cell.annot.lst <- lapply(split(cell.annot, f = cell.annot$ctype), function(x) x$cell)

  vec.prob.ctype <- MeanAcrossClusters(mat.prob, cell.annot.lst)

  dat.overlap.est.lst <- lapply(jctypes, function(jctype){
    dat.overlap.est <- data.frame(Bin = names(vec.prob.ctype[[jctype]]),
                                  overlap.estimate = vec.prob.ctype[[jctype]],
                                  ctype = jctype,
                                  stringsAsFactors = FALSE)
    return(dat.overlap.est)
  })

  bin.means.merged.est.lst <- lapply(jctypes, function(jctype){
    left_join(bin.means.merged.lst[[jctype]], dat.overlap.est.lst[[jctype]], by = c("Bin", "ctype"))
  })

  m.check <- lapply(jctypes, function(jctype){
    m <- ggplot(bin.means.merged.est.lst[[jctype]], aes(x = BinMean.dbl, y = overlap.estimate)) +
      geom_point(alpha = 0.1) +
      theme_bw() +
      ggtitle(paste("Bin overlap inference.\nCelltype:", jctype)) +
      xlab("Ground Truth Overlap") +
      ylab("Estimated Overlap from scChIX") +
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    return(m)
  })

  JFuncs::multiplot(m.check[[1]], m.check[[2]], m.check[[3]], cols = 3)



  # get average heatmap across celltypes

  LL.mat.avg.byctype <- lapply(jctypes, function(jctype){
    jmat.lst <- lapply(fits.out[cell.annot.lst[[jctype]]], function(x){
      jmat.out <- x$ll.mat
    })
    jmat.avg <- Reduce("+", jmat.lst) / length(jmat.lst)
  })

  LL.long.byctype <- lapply(LL.mat.avg.byctype, function(jmat){
    jlst <- lapply(as.data.frame(jmat) %>% as.list(), function(jvec){
      names(jvec) <- paste("mark1", letters[1:3], sep = "_")
      return(jvec)
    })
    names(jlst) <- paste("mark2", letters[4:6], sep = "_")

    mark2.cnames <- names(jlst); names(mark2.cnames) <- mark2.cnames
    jdat.long <- lapply(mark2.cnames, function(mark2.cname){
      x <- jlst[[mark2.cname]]
      jdat.tmp <- data.frame(Mark1_ClusterName = names(x), Mark2_ClusterName = mark2.cname, LL = x, stringsAsFactors = FALSE)
    }) %>%
      bind_rows()
    return(jdat.long)
  })

  m.LL.heatmap.byctype <- lapply(jctypes, function(jctype){
    jdat.long.tmp <- LL.long.byctype[[jctype]]
    m.tmp <- ggplot(jdat.long.tmp, aes(x = Mark1_ClusterName, y = Mark2_ClusterName, fill = LL)) +
      geom_tile() +
      scale_fill_viridis_c() +
      theme_bw() +
      xlab("") +
      ylab("") +
      ggtitle(jctype) +
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom")
    return(m.tmp)
  })

  print(m.LL.heatmap.byctype)
  JFuncs::multiplot(m.LL.heatmap.byctype[[1]], m.LL.heatmap.byctype[[2]], m.LL.heatmap.byctype[[3]], cols = 3)

  dev.off()

  saveRDS(dat.umap.joined, file = outrds1)
  saveRDS(bin.means.merged.est.lst, file = outrds2)
  saveRDS(LL.mat.avg.byctype, file = outrds3)
}



