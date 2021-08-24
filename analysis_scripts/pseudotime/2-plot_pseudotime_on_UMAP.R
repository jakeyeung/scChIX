# Jake Yeung
# Date of Creation: 2021-08-10
# File: ~/projects/scChIX/analysis_scripts/pseudotime/2-plot_pseudotime_on_UMAP.R
#

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

library(topicmodels)

library(scchicFuncs)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123


hubprefix <- "/home/jyeung/hub_oudenaarden"


cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

# Load matrix  ------------------------------------------------------------

jdate <- "2021-07-23"
jquant <- "manual2nocenter"

jmark1 <- "K36"; jmark2 <- "K9m3"; jmarks <- c(jmark1, jmark2); jmarkdbl <- paste(jmark1, jmark2, sep = "-")

names(jmarks) <- jmarks
jmarkdbl <- paste(c(jmark1, jmark2), collapse = "-")

jstr <- paste(c(jmarks, jmarkdbl), collapse = "_")

jprefix <- "var_filtered"
jname <- paste(jprefix, jquant, jstr, sep = "_")

infrdata <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/scchix_pipeline_from_LDA/scchix_outputs_objs/", jname, "/unmix_scchix_inputs_clstr_by_celltype_", jmarkdbl, ".removeNA_FALSE.RData"))
assertthat::assert_that(file.exists(infrdata))



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




# Load data  --------------------------------------------------------------


out.lst.lst <- lapply(jmarks, function(jmarktmp){
  load(infs[[jmarktmp]], v=T)  # out.objs, out.lda.predict
  out.lst <- list(out.objs = out.objs, out.lda.predict = out.lda.predict)
})

tm.lst <- lapply(out.lst.lst, function(x) posterior(x$out.objs$out.lda))

umap.objs.lst <- lapply(jmarks, function(jmarktmp){
  tm.orig <- tm.lst[[jmarktmp]]
  umap.out <- umap(tm.orig$topics, config = jsettings)
  return(umap.out)
})

dat.umap.merge.lst <- lapply(jmarks, function(jmarktmp){
  umap.out <- umap.objs.lst[[jmarktmp]]
  tm.orig <- tm.lst[[jmarktmp]]
  out.lda.predict <- out.lst.lst[[jmarktmp]]$out.lda.predict
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



# Calculate imputed values  -----------------------------------------------


terms.mat <- tm.lst$K9m3$terms
topics.mat1 <- tm.lst$K9m3$topics
topics.mat2 <- out.lst.lst$K9m3$out.lda.predict$topics
topics.mat.merged <- rbind(topics.mat1, topics.mat2)

dat.impute.log2 <- t(log2(topics.mat.merged %*% terms.mat))

# plot top hit
jbin <- "chr3:106150000-106200000"
jbin <- "chr2:177500000-177550000"
jbin <- "chr9:3000000-3050000"

jvec <- data.frame(cell = colnames(dat.impute.log2), value = dat.impute.log2[jbin, ], stringsAsFactors = FALSE)

dat.check <- left_join(dat.umap.merge.lst$K9m3, jvec)

ggplot(dat.check, aes(x = umap1, y = umap2, color = value)) +
  geom_point() +
  theme_bw() +
  scale_color_viridis_c() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

inf.ptime <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis/metadata/from_scchix_downstream_plots/metadata_pseudotime_K36-K9me3.K9m3.2021-08-10.txt")
dat.ptime <- fread(inf.ptime)

cells.keep <- dat.ptime$cell

ggplot(dat.check %>% filter(cell %in% cells.keep & umap2 > -2 & type == "dbl"), aes(x = umap1, y = umap2, color = value)) +
  geom_point() +
  theme_bw() +
  scale_color_viridis_c() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~stage)

cells.single <- rownames(topics.mat1)

stages.vec <- sapply(cells.single, function(x) strsplit(x, split = "-")[[1]][[1]])

unique(stages.vec)

ggplot(dat.check %>% filter(cell %in% cells.keep & umap2 > -2 & type == "single"), aes(x = umap1, y = umap2, color = value)) +
  geom_point() +
  theme_bw() +
  scale_color_viridis_c() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~stage)


ptime.hash <- hash::hash(dat.ptime$cell, dat.ptime$ptime)

# Try the prob mat --------------------------------------------------------

inf.prob <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/tagged_bams/merged/split_bams/K36-K9m3_mapq_40/probabilityMatrix_linearInterpolated.csv")
dat.prob <- fread(inf.prob)


jbin.i <- which(dat.prob$V1 == jbin)

jbin.range <- seq(jbin.i - 5, jbin.i + 5)

mat.prob <- as.matrix(dat.prob %>% dplyr::select(-V1))
rownames(mat.prob) <- dat.prob$V1
mat.prob <- t(mat.prob)

# sort cells by pseudotime
cells.keep.i <- which(rownames(mat.prob) %in% cells.keep)
mat.prob.filt <- mat.prob[cells.keep.i, ]

heatmap(mat.prob.filt[, jbin.range], Rowv = NA, Colv = NA)


# plot log2 in heatmap



jbin <- "chr3:106150000-106200000"
jbin <- "chr9:3000000-3050000"
jbin <- "chr2:177500000-177550000"
jbins <- c("chr2:177500000-177550000", "chr3:106150000-106200000", "chr9:3000000-3050000")
jrange <- 15

for (jbin in jbins){

  jbin.i2 <- which(rownames(dat.impute.log2) == jbin)
  jbin.range2 <- seq(jbin.i2 - jrange, jbin.i2 + jrange)

  cells.keep.i2 <- which(colnames(dat.impute.log2) %in% cells.keep)
  mat.impute.filt <- t(dat.impute.log2[jbin.range2, cells.keep.i2])
  order.i <- order(sapply(rownames(mat.impute.filt), function(x) ptime.hash[[x]]), decreasing = TRUE)
  mat.impute.filt <- mat.impute.filt[order.i, ]

  rowlabs <- sapply(rownames(mat.impute.filt), function(x) strsplit(x, split = "-")[[1]][[1]])
  rowlabs.i <- as.numeric(factor(rowlabs, levels = c("E9p5", "E10p5", "E11p5")))
  rowcols <- sapply(rowlabs.i, function(i) cbPalette[[i]])

  # get ptime from colnames

  rowsptime <- sapply(rownames(mat.impute.filt), function(x) ptime.hash[[x]])

  plot(rowsptime)
  plot(dat.ptime$ptime)

  mat.impute.filt <- scale(mat.impute.filt, center = TRUE, scale = FALSE)


  jbinstr <- gsub(":", "-", jbin)
  outpdf <- file.path(paste0("/home/jyeung/hub_oudenaarden/jyeung/data/tmp/heatmap_K9me3_pseudotime.", jbinstr, ".pdf"))
  pdf(outpdf, useDingbats = FALSE)
  heatmap3::heatmap3(mat.impute.filt, Rowv = NA, Colv = NA, RowSideColors = rowcols, main = jbin, scale = "none")
  dev.off()
}


# Show imputed values on UMAP  --------------------------------------------




