# Jake Yeung
# Date of Creation: 2021-07-16
# File: ~/projects/scChIX/analysis_scripts/4-analyze_scChIX_outputs_projections_FromPipelineGenomewide.R
# description

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)

library(ggforce)

library(JFuncs)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

hubprefix <- "/home/jyeung/hub_oudenaarden"


# Load matrix  ------------------------------------------------------------

# jquants <- c("0.15", "0.2")
# jquants <- c("0.2")
# jquants <- c("0.15")
# jquants <- c("0.15")

# jdate <- "2021-07-16"
# jquant <- "0.15"

# jdate <- "2021-07-19"
# jquant <- "manual"

# jdate <- "2021-07-20"
# jquant <- "manual2"

# jdate <- "2021-07-22"
# jquant <- "manual2noblood"

# jdate <- "2021-07-23"
# jquant <- "manual2nocenter"

# jdate <- "2021-07-24"
# jquant <- "manual2nocenternoE8"

jdate <- "2021-07-29"
jquant <- "manual2nocenternoE8unifyK36"


jmark1 <- "K36"; jmark2 <- "K27"; jmarks <- c(jmark1, jmark2); jmarkdbl <- paste(jmark1, jmark2, sep = "-")
# jmark1 <- "K36"; jmark2 <- "K9m3"; jmarks <- c(jmark1, jmark2); jmarkdbl <- paste(jmark1, jmark2, sep = "-")

# for (jquant in jquants){
  names(jmarks) <- jmarks
  jmarkdbl <- paste(c(jmark1, jmark2), collapse = "-")

  jstr <- paste(c(jmarks, jmarkdbl), collapse = "_")

  jprefix <- "var_filtered"
  jname <- paste(jprefix, jquant, jstr, sep = "_")

  infrdata <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/scchix_pipeline_from_LDA/scchix_outputs_objs/", jname, "/unmix_scchix_inputs_clstr_by_celltype_", jmarkdbl, ".removeNA_FALSE.RData"))
  assertthat::assert_that(file.exists(infrdata))

  outmain <- "/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/scchix_downstream_plots"
  outdir <- file.path(outmain, jname)
  dir.create(outdir)


  outpdf <- file.path(outdir, paste0("scchix_downstream_", jstr, ".", Sys.Date(), ".pdf"))
  pdf(outpdf, useDingbats = FALSE)

  # jsuffix <- "dbl_k36_k9m3_cleaned"
  infs <- lapply(jmarks, function(jmark){
    inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/LDA_scchix_outputs/from_pipeline/", jname, "/", jstr, "/Gastru_Unmixed_DblMark.", jname, ".", jmark, ".RData")
    assertthat::assert_that(file.exists(inf))
    return(inf)
  })


  # load mats
  # metamain <- file.path(paste0("/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_", jprefix))
  metamain <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/scchix_pipeline_from_LDA/objs_from_LDA/", jname))
  assertthat::assert_that(dir.exists(metamain))

  dat.meta.lst <- lapply(jmarks, function(jmark){
    infmeta <- file.path(metamain, paste0("celltyping_output_filt.", jmark, ".", jdate, ".rds"))
    print(infmeta)
    dat.meta <- readRDS(infmeta)
    return(dat.meta)
  })


  dat.umap.merge.lst <- lapply(jmarks, function(jmarktmp){

    # jmarktmp <- jmarks[[1]]
    load(infs[[jmarktmp]], v=T)

    tm.orig <- posterior(out.objs$out.lda)
    umap.out <- umap(tm.orig$topics, config = jsettings)
    dat.umap.orig <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE)
    dat.umap.orig <- DoLouvain(topics.mat = tm.orig$topics, custom.settings.louv = jsettings, dat.umap.long = dat.umap.orig)
    dat.umap.orig.annot <- left_join(dat.umap.orig %>% mutate(type = "single") %>% dplyr::select(-louvain),
                                     dat.meta.lst[[jmarktmp]] %>% dplyr::select(c(cell, cluster)), by = "cell")


    # add projections
    umap.out.pred.layout <- predict(umap.out, data = out.lda.predict$topics)
    dat.umap.pred.annot <- data.frame(cell = rownames(umap.out.pred.layout), umap1 = umap.out.pred.layout[, 1], umap2 = umap.out.pred.layout[, 2], stringsAsFactors = FALSE) %>%
      mutate(type = "dbl",
             cluster = "na")

    dat.umap.merge <- rbind(dat.umap.orig.annot, dat.umap.pred.annot)
    dat.umap.merge$mark <- jmarktmp
    dat.umap.merge <- dat.umap.merge %>%
      rowwise() %>%
      mutate(stage = as.character(strsplit(cell, split = "-")[[1]][[1]]))
    return(dat.umap.merge)
  })


  cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
  m.lst <- lapply(jmarks, function(jmarktmp){
    dat.umap.merge <- dat.umap.merge.lst[[jmarktmp]]
    m <- ggplot(dat.umap.merge, aes(x = umap1, y = umap2, color = cluster)) +
      geom_point() +
      facet_wrap(~type) +
      theme_bw() +
      scale_color_manual(values = cbPalette) +
      ggtitle(jmarktmp) +
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  })

  multiplot(m.lst[[1]], m.lst[[2]], cols = 1)


  m.stage.lst <- lapply(jmarks, function(jmarktmp){
    dat.umap.merge <- dat.umap.merge.lst[[jmarktmp]]
    m <- ggplot(dat.umap.merge, aes(x = umap1, y = umap2, color = stage)) +
      geom_point() +
      facet_wrap(~type) +
      theme_bw() +
      scale_color_manual(values = cbPalette) +
      ggtitle(jmarktmp) +
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  })

  multiplot(m.stage.lst[[1]], m.stage.lst[[2]], cols = 1)


  m.lst2 <- lapply(jmarks, function(jmarktmp){
    dat.umap.merge <- dat.umap.merge.lst[[jmarktmp]]
    m <- ggplot(dat.umap.merge, aes(x = umap1, y = umap2, color = cluster)) +
      geom_point() +
      # facet_wrap(~type) +
      theme_bw() +
      scale_color_manual(values = cbPalette) +
      ggtitle(jmarktmp) +
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  })
  multiplot(m.lst2[[1]], m.lst2[[2]], cols = 1)



  # Split things up  --------------------------------------------------------


  dat.merge.rbind <- bind_rows(dat.umap.merge.lst) %>%
    group_by(mark) %>%
    mutate(umap1.scale = scale(umap1, center = TRUE, scale = TRUE),
           umap2.scale = scale(umap2, center = TRUE, scale = TRUE),
           umap1.shift = ifelse(mark == "K36", umap1.scale - 5, umap1.scale + 5)) %>%
    rowwise() %>%
    mutate(stage = strsplit(cell, split = "-")[[1]][[1]])


  ggplot(dat.merge.rbind, aes(x = umap1.shift, y = -1 * umap2.scale, group = cell)) +
    geom_point() +
    ggtitle("Double + single cells") +
    geom_path(alpha = 0.05) +
    theme_bw() +
    theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  ggplot(dat.merge.rbind, aes(x = umap1.shift, y = 1 * umap2.scale, group = cell)) +
    geom_point() +
    ggtitle("Double + single cells") +
    geom_path(alpha = 0.05) +
    theme_bw() +
    theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  ggplot(dat.merge.rbind %>% filter(type == "dbl"), aes(x = umap1.shift, y = umap2.scale, group = cell)) +
    geom_point() +
    ggtitle("Double-incubated cells only") +
    geom_path(alpha = 0.05) +
    theme_bw() +
    theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  ggplot(dat.merge.rbind %>% filter(type == "dbl"), aes(x = umap1.shift, y = umap2.scale, group = cell, color = stage)) +
    geom_point() +
    ggtitle("Double-incubated cells only") +
    geom_path(alpha = 0.05) +
    facet_wrap(~stage) +
    theme_bw() +
    theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  ggplot(dat.merge.rbind, aes(x = umap1.shift, y = umap2.scale, group = cell, color = stage)) +
    geom_point() +
    ggtitle("Double + single cells") +
    geom_path(alpha = 0.05) +
    facet_wrap(~stage) +
    theme_bw() +
    theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#   cells.keep <- subset(dat.merge.rbind, umap2.scale < -1)$cell
#   cells.keep <- subset(dat.merge.rbind, mark == "K9m3" & cluster == "cluster5")$cell
#   cells.keep <- subset(dat.merge.rbind, umap2.scale > 2 & umap1.shift > 0)$cell
#
#   ggplot(dat.merge.rbind  %>% mutate(highlight = cell %in% cells.keep) %>% arrange(highlight), aes(x = umap1.shift, y = umap2.scale, group = cell, color = highlight)) +
#     geom_point() +
#     geom_path(alpha = 0.05) +
#     theme_bw() +
#     theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#
#   ggplot(dat.merge.rbind  %>% filter(cell %in% cells.keep), aes(x = umap1.shift, y = umap2.scale, group = cell)) +
#     geom_point() +
#     geom_path(alpha = 0.05) +
#     theme_bw() +
#     theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#

  # Plot 2D -----------------------------------------------------------------

  load(infrdata, v=T)

  fits.out <- act.repress.coord.lst
  w.lst <- sapply(fits.out, function(x) x$w)


  cell.vec <- names(fits.out)
  names(cell.vec) <- cell.vec

  # Check louvins  ----------------------------------------------------------


  # if louvains are now from clusters need eto rethink jcoord
  coords.dbl <- lapply(cell.vec, function(jcell){
    jfit <- fits.out[[jcell]]
    jweight <- fits.out[[jcell]]$w
    p.mat <- SoftMax(jfit$ll.mat)
    jcoord <- which(jfit$ll.mat == max(jfit$ll.mat), arr.ind = TRUE)
    jmax <- max(p.mat)

    # rows are active, columns are repress I THINK?
    # TODO: assumes underscores be careful!
    jlouv.act <- rownames(p.mat)[[jcoord[[1]]]]
    jlouv.repress <- colnames(p.mat)[[jcoord[[2]]]]

    if (grepl("_", jlouv.act)){
      jlouv.act <- strsplit(jlouv.act, split = "_")[[1]][[2]]
    }
    if (grepl("_", jlouv.repress)){
      jlouv.repress <- strsplit(jlouv.repress, split = "_")[[1]][[2]]
    }
    out.dat <- data.frame(cell = jcell, louv.act = jlouv.act, louv.repress = jlouv.repress, lnprob = jmax, w = jweight, stringsAsFactors = FALSE)
    return(out.dat)
  }) %>%
    bind_rows()


  m.grid <- ggplot(coords.dbl, aes(x = louv.act, y = louv.repress, color = w)) +
    geom_point(alpha = 0.25, position = ggforce::position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
    theme_bw() +
    scale_color_viridis_c() +
    theme(aspect.ratio=0.6) +
    ggtitle("Each dot is a double stained cell,\nX-Y shows the cluster pair it is assigned")
  print(m.grid)


  dev.off()


# }



