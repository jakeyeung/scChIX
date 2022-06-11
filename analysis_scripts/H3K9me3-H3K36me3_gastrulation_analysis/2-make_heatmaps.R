# Jake Yeung
# Date of Creation: 2021-09-28
# File: ~/projects/scChIX/analysis_scripts/H3K9me3-H3K36me3_gastrulation_analysis/2-make_heatmaps.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(heatmap3)

hubprefix <- "/Users/yeung/hub_oudenaarden"
jmarks <- c("K36", "K9m3"); names(jmarks) <- jmarks

# Load meta ---------------------------------------------------------------

outpdf <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis_macbook/heatmaps_downstream/K36_K9me3_RNAseq_heatmaps.", Sys.Date(), ".pdf"))
pdf(outpdf, useDingbats = FALSE)

inf.meta <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis_macbook/objs_from_macbook/metadata_flipped.rds"))
dat.meta.all.flipped <- readRDS(inf.meta)

# Load tm results ---------------------------------------------------------

inf.tm <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis_macbook/objs_from_macbook/tm_results_merged_lst.rds")
tm.lst <- readRDS(inf.tm)


# Get imputed -------------------------------------------------------------

dat.imputed.lst <- lapply(tm.lst, function(tm){
  t(log2(tm$topics %*% tm$terms))
})


# Load genes and  bins -------------------------------------------------------------

inf.genes <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis_macbook/objs_from_macbook/top_bins_genes.keeptop_250.rds")
annot.genes <- readRDS(inf.genes)

inf.annot <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis_macbook/objs_from_macbook/bin_annotations_from_biomart.rds"))
bins.annot <- readRDS(inf.annot)

bins2genes <- hash::hash(bins.annot$out2.df.closest$region_coord, bins.annot$out2.df.closest$gene)

inf.topicsannot <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis_macbook/objs_from_macbook/topics_celltypes_annotation.rds")
topicsannot <- readRDS(inf.topicsannot)

topics.hash <- hash::hash(topicsannot)
topics.hash.inv <- hash::invert(topics.hash)

# Load RNA-seq data  ------------------------------------------------------

inf.rnaseq <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis_macbook/objs_from_macbook/cao_rnaseq_reference_pbulk_zscore_ctypes_merged.rds")
dat.ref <- readRDS(inf.rnaseq)



# Make heatmap  -----------------------------------------------------------


clstrs.order <- c("Erythroid", "WhiteBloodCells", "Endothelial", "NeuralTubeNeuralProgs", "Neurons", "SchwannCellPrecursor", "Epithelial", "Stromal", "ConnectiveTissueProg")
names(clstrs.order) <- clstrs.order
topics.order <- sapply(clstrs.order, function(x) topics.hash.inv[[x]])

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
cbPalette.sub <- cbPalette[1:length(clstrs.order)]
cols.hash <- hash::hash(clstrs.order, cbPalette.sub)

bins.vec <- unlist(lapply(annot.genes[topics.order], function(jlst) jlst$bins.keep))

bins.ctype.colvec <- lapply(clstrs.order, function(clstr){
  bin.vec <- annot.genes[[topics.order[[clstr]]]]$bins.keep
  assertthat::assert_that(length(bin.vec) >= 1)
  jcol <- cols.hash[[clstr]]
  colvec <- rep(jcol, length(bin.vec))
  names(colvec) <- bin.vec
  return(colvec)
}) %>%
  unlist() 

genes.ctype.colvec <- lapply(clstrs.order, function(clstr){
  gene.vec <- annot.genes[[topics.order[[clstr]]]]$genes.keep
  assertthat::assert_that(length(gene.vec) >= 1)
  jcol <- cols.hash[[clstr]]
  colvec <- rep(jcol, length(gene.vec))
  names(colvec) <- gene.vec
  return(colvec)
}) %>%
  unlist()
gene2col <- hash::hash(sapply(names(genes.ctype.colvec), function(x) strsplit(x, split = "\\.")[[1]][[2]]), genes.ctype.colvec)


# order cells by celltypes
dat.ordered.lst <- lapply(jmarks, function(jmark){
  dat.ordered <- dat.meta.all.flipped %>%
    ungroup() %>%
    filter(mark == jmark) %>%
    mutate(cluster = gsub("Precusor", "Precursor", cluster)) %>%
    mutate(cluster = factor(cluster, levels = clstrs.order)) %>%
    arrange(cluster) %>%
    rowwise() 
  dat.ordered$colorcode <- sapply(as.character(dat.ordered$cluster), function(x) cols.hash[[x]])
  return(dat.ordered)
})

# dat.ordered <- dat.meta.all.flipped %>%
#   ungroup() %>%
#   filter(mark == "K36") %>%
#   mutate(cluster = gsub("Precusor", "Precursor", cluster)) %>%
#   mutate(cluster = factor(cluster, levels = clstrs.order)) %>%
#   arrange(cluster) %>%
#   rowwise() 
# dat.ordered$colorcode <- sapply(as.character(dat.ordered$cluster), function(x) cols.hash[[x]])
# cells.ordered <- dat.ordered$cell

cols.color.ordered.lst <- lapply(dat.ordered.lst, function(jdat){
  jdat$colorcode
})

# cols.color.ordered <- dat.ordered$colorcode
rows.color.ordered <- rev(bins.ctype.colvec)

mat.arranged.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  dat.imputed.lst[[jmark]][rev(bins.vec), dat.ordered.lst[[jmark]]$cell]
})

mat.arranged.merged <- do.call(cbind, mat.arranged.lst)
jfloor <- -16
mat.arranged.merged[which(mat.arranged.merged <= jfloor)] <- jfloor

# heatmap3::heatmap3()


# label cell types and gene lists
jmarktmp <- "K9m3"
heatmap3::heatmap3(mat.arranged.lst[[jmarktmp]], Rowv = NA, Colv = NA, labRow = NA, labCol = NA, scale = "row", 
                   ColSideColors = cols.color.ordered.lst[[jmarktmp]], RowSideColors = rows.color.ordered, 
                   main = paste(jmarktmp, "celltyping"))

# plot both
heatmap3::heatmap3(mat.arranged.merged, Rowv = NA, Colv = NA, labRow = NA, labCol = NA, scale = "row", 
                   ColSideColors = unlist(cols.color.ordered.lst), RowSideColors = rows.color.ordered, 
                   main = paste("Celltyping, both celltyping"))


# Scale using H3K36me3, then apply that scaling to H3K9me3  ---------------

jfloor2 <- -16
mat.arranged.floor.lst <- lapply(mat.arranged.lst, function(jmat){
  jmat[which(jmat < jfloor2)] <- jfloor2
  return(jmat)
})

mat.arranged.k36.scale <- t(scale(t(mat.arranged.floor.lst$K36), center = TRUE, scale = TRUE))

jcenter <- attr(mat.arranged.k36.scale, "scaled:center")
jscale <- attr(mat.arranged.k36.scale, "scaled:scale")

# apply to k9
mat.arranged.k9.scale <- sweep(mat.arranged.floor.lst$K9m3, MARGIN = 1, STATS = jcenter, FUN = "-")
mat.arranged.k9.scale <- sweep(mat.arranged.k9.scale, MARGIN = 1, STATS = jscale, FUN = "*")

heatmap3::heatmap3(mat.arranged.k36.scale, Rowv = NA, Colv = NA, labRow = NA, labCol = NA, scale = "none", 
                   ColSideColors = cols.color.ordered.lst$K36, RowSideColors = rows.color.ordered, 
                   main = paste("K36", "celltyping"))

heatmap3::heatmap3(mat.arranged.k9.scale, Rowv = NA, Colv = NA, labRow = NA, labCol = NA, scale = "none", 
                   ColSideColors = cols.color.ordered.lst$K9, RowSideColors = rows.color.ordered, 
                   main = paste("K9m3", "celltyping"))

# Show RNA-seq data gets celltype properly  -------------------------------

# create pbulk RNAseq mat
mat.ref <- reshape2::dcast(data = dat.ref, formula = gene ~ celltype, value.var = "zscore")
rownames(mat.ref) <- mat.ref$gene
mat.ref$gene <- NULL


bins.vec.ordered <- rownames(mat.arranged.k36.scale)  # already reversed! So reverse the color scale
genes.vec.ordered <- sapply(bins.vec.ordered, function(b) AssignHash(x = b, jhash = bins2genes, null.fill = NA))
genes.vec.ordered.nona <- genes.vec.ordered[!is.na(genes.vec.ordered)]
cols.vec.ordered <- bins.ctype.colvec[!is.na(genes.vec.ordered)]
genes.vec.ordered.nona.filt <- genes.vec.ordered.nona[genes.vec.ordered.nona %in% rownames(mat.ref)]
cols.vec.ordered <- cols.vec.ordered[genes.vec.ordered.nona %in% rownames(mat.ref)]

rows.keep <- genes.vec.ordered.nona.filt

mat.ref.ordered <- mat.ref[rows.keep, ]
print(colnames(mat.ref.ordered))
pbulk.order <- c("Erythroid", 
                 "X31.White.blood.cells", 
                 "X20.Endothelial.cells", 
                 "NeuralTubeAndProgs", 
                 "NeuronalProgenitors", 
                 "X23.Schwann.cell.precursor", 
                 "X6.Epithelial.cells", 
                 "MesenchymalProgenitors", 
                 "X34.Cardiac.muscle.lineages")

pbulk.others <- colnames(mat.ref.ordered)[!colnames(mat.ref.ordered) %in% pbulk.order]
cnames.ordered <- c(pbulk.order, pbulk.others)

mat.ref.ordered <- mat.ref[rows.keep, cnames.ordered]

# get cols for rows and cols 
rows.color.ordered.genes <- rev(cols.vec.ordered)  # need to be reversed because rows already reversed
cols.color.ordered.genes <- cbPalette[1:length(cnames.ordered)]

heatmap3::heatmap3(mat.ref.ordered, Rowv = NA, Colv = NA, labRow = NA, labCol = NA, scale = "none", 
                   ColSideColors = cols.color.ordered.genes, RowSideColors = rows.color.ordered.genes, 
                   main = paste("RNAseq", "celltyping"))

# get colorcode 

dat.ordered.lst2 <- dat.ordered.lst
dat.ordered.lst2$K9m3$k9clsts <- sapply(dat.ordered.lst2$K9m3$cluster, function(x) ifelse(!x %in% c("Erythroid", "WhiteBloodCells"), "Other", x)) 

# shuffle rows, sort by clsts
dat.ordered.lst2$K9m3 <- dat.ordered.lst2$K9m3[sample(rownames(dat.ordered.lst2$K9m3), size = nrow(dat.ordered.lst2$K9m3), replace = FALSE), ]

outrds.meta.colored <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis_macbook/heatmaps_downstream/dat_meta_ordered_colorcoded.rds")
saveRDS(dat.ordered.lst2, file = outrds.meta.colored)

jlong <- dat.ordered.lst2 %>% bind_rows()
ggplot(jlong,  
       aes(x = umap1.shift2, y = umap2.flip, color = colorcode, group = cell)) + 
  geom_point() + 
  geom_path(alpha = 0.01) + 
  scale_color_identity() + 
  theme_bw() + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme(aspect.ratio=0.33, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.meta.all.flipped,  
       aes(x = umap1.shift2, y = umap2.flip, color = cluster, group = cell)) + 
  geom_point() + 
  geom_path(alpha = 0.01) + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme(aspect.ratio=0.33, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

dev.off()






