# Jake Yeung
# Date of Creation: 2021-07-03
# File: ~/projects/scChIX/analysis_scripts/3-get_clusters_K36.R
# Get clusters K36
# use 50kb?

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

nn <- 15
jsettings <- umap.defaults
jsettings$n_neighbors <- nn
jsettings$min_dist <- 0.2
jsettings$random_state <- 123
jsettings$spread <- 5

outdir <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_NN_", nn, "_check_plates")
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

infs <- lapply(jmarks, function(jmark){
  print(jmark)
  inf <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/LDA_outputs/ldaAnalysis_", jsuffix, "/lda_outputs.count_tables.", jsuffix, ".", jmark, ".2021-06-28.K-30.binarize.FALSE/ldaOut.count_tables.", jsuffix, ".", jmark, ".2021-06-28.K-30.Robj"))
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
# jmark.test <- "K9m3"
# dat.umap.test <- DoUmapAndLouvain(tm.result.lst[[jmark.test]]$topics, jsettings = jsettings2)
# ggplot(dat.umap.test, aes(x = umap1, y = umap2, color = louvain)) +
#   geom_point() +
#   theme_bw() + ggtitle(jmark.test) +
#   scale_color_manual(values = cbPalette) +
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Assign clusters using topics --------------------------------------------
# K36

# 1) Erythroids (topic10)
# 2) Neuron progenitors (topic11)
# 3) Epithelial cells? (topic21)
# 4) Endothelial cells (topic17)
# 5) Cardiac muscle progenitors? (topic5)
# 6) Neural progenitors, neural tube (topic1)
# 7) Chondrocyte progenitors (topic6)
# 8) White blood cells (topic13)

topics.keep <- c("topic10", "topic11", "topic21", "topic17", "topic5", "topic1", "topic6", "topic13")
topics.keep.name <- c("Erythroid", "NeuronProgs", "Epithelial", "Endothelial", "CardiacMuscProg", "NeuralTube", "ChondrocyteProgs", "WhiteBlood")
names(topics.keep.name) <- topics.keep
assertthat::assert_that(length(topics.keep) == length(topics.keep.name))



# Assign cells using Gaussian mixture models  -----------------------------


jthres <- 0.5
# label cells by topic
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

jmarktmp <- "K36"


for (jmarktmp in jmarks){


  tm.result <- tm.result.lst[[jmarktmp]]
  # dat.umap.long <- dat.umap.lst[[jmarktmp]]

  # add topic info
  dat.topics <- data.frame(cell = rownames(tm.result$topics), tm.result$topics, stringsAsFactors = FALSE)
  dat.umap.long <- left_join(dat.umap.lst[[jmarktmp]], dat.topics)


  pdf(file.path(outdir, paste0("celltyping_alltopics_from_", jmarktmp, ".", Sys.Date(), ".pdf")))

  topics.all <- colnames(tm.result$topics)
  # mm.celltype.lst <- lapply(topics.keep, function(jtopic){
  mm.celltype.lst <- lapply(topics.all, function(jtopic){
    print(jtopic)
    # plot with topic loadings
    m.umap <- PlotXYWithColor(dat.umap.long, xvar = "umap1", yvar = "umap2", cname = jtopic) + scale_color_viridis_c()
    print(m.umap)

    tvec.raw <- sort(tm.result$topics[, jtopic])
    # transform
    tvec <- log(tvec.raw / (1 - tvec.raw))
    # xline <- quantile(tvec, )
    mm <- normalmixEM(x = tvec, lambda = c(0.9, 0.1), mu = c(-5, -1), sigma = c(2, 1), k = 2)
    (xline <- min(mm$x[which(mm$posterior[, 1] < jthres)]))
    xcells <- names(tvec)[which(tvec > xline)]
    print(paste(length(xcells), "/", length(tvec), "assigned to", jtopic))
    # get topic value for assigned cells
    tvec.raw.filt <- tvec.raw[xcells]

    (xline <- min(mm$x[which(mm$posterior[, 1] < jthres)]))
    # xline <- max(mm$x[indx.btwn][which(post.filt[, 2] > 0.5)])
    plot(density(tvec), main = paste(jtopic, jthres), xlab = "Log Odds [log(p / (1 - p))]")
    abline(v = xline, col = 'blue')
    plot.mixEM(mm, whichplots = 2, xlab2 = "Log Odds [log(p / (1 - p))]", main2 = paste(jtopic, jthres))
    abline(v = xline, col = 'blue')

    cells.keep <- xcells
    m.check <- PlotXYWithColor(dat.umap.long %>% mutate(is.celltype = cell %in% cells.keep), xvar = "umap1", yvar = "umap2", cname = "is.celltype", jtitle = paste(jtopic, jthres), cont.color = FALSE, col.palette = cbPalette)
    print(m.check)

    return(list(topic = jtopic, topic.weight = tvec.raw.filt, celltype = xcells, mm = mm, threshold = xline))
  })


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

  outrds.meta <- file.path(outdir, paste0("celltyping_output.", jmarktmp, ".", Sys.Date(), ".rds"))
  saveRDS(dat.merge, file = outrds.meta)

  outrds.countmat <- file.path(outdir, paste0("countmat_output.", jmarktmp, ".", Sys.Date(), ".rds"))
  saveRDS(count.mat.lst[[jmarktmp]], file = outrds.countmat)

  outrds.outlda <- file.path(outdir, paste0("lda_output.", jmarktmp, ".", Sys.Date(), ".rds"))
  saveRDS(out.lda.lst[[jmarktmp]], file = outrds.outlda)

  # outrdata <- file.path(outdir, paste0("celltyping_output.", jmarktmp, ".", Sys.Date(), ".RData"))

  # save RData with dat.merge and out.lda


}






# # plot all cells assigned to at least a toipc
# dat.ctypes <- lapply(mm.celltype.lst, function(jdat){
#   jtop <- jdat$topic
#   print(jtop)
#   cells.keep <- jdat$celltype
#   if (length(cells.keep) == 0){
#     jdat.out <- data.frame(NULL)
#   } else {
#     jdat.out <- data.frame(topic = jtop, cell = cells.keep, celltype = topics.keep.name[[jtop]], stringsAsFactors = FALSE)
#   }
#   return(jdat.out)
# }) %>%
#   bind_rows()
#
# dat.ctypes.merge <- left_join(dat.umap.long, dat.ctypes)
#
# m <- ggplot(dat.ctypes.merge, aes(x = umap1, y = umap2, color = celltype)) +
#   geom_point() +
#   theme_bw() +
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   scale_color_manual(values = cbPalette, na.value = "grey50")
# print(m)
#
# m <- ggplot(dat.ctypes.merge, aes(x = umap1, y = umap2, color = topic)) +
#   geom_point() +
#   theme_bw() +
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   scale_color_manual(values = cbPalette, na.value = "grey50")
# print(m)


