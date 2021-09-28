# Jake Yeung
# Date of Creation: 2021-08-25
# File: ~/projects/scChIX/analysis_scripts/7-check_total_counts_on_UMAP.R
# Eryths have different K36 to K9 ratio: is it due to K36 or due to K9?

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)



# Load raw counts ---------------------------------------------------------


jquant <- "manual2nocenterfilt2"
jmark1 <- "K36"; jmark2 <- "K9m3"; jmarks <- c(jmark1, jmark2); jmarkdbl <- paste(jmark1, jmark2, sep = "-")

names(jmarks) <- jmarks
jmarkdbl <- paste(c(jmark1, jmark2), collapse = "-")

jstr <- paste(c(jmarks, jmarkdbl), collapse = "_")

jprefix <- "var_filtered"
jname <- paste(jprefix, jquant, jstr, sep = "_")



inf.lda.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.lda.tmp <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/dblchic/gastrulation/LDA_scchix_outputs/from_pipeline_unmixed_singles_LDA_together/var_filtered_manual2nocenterfilt2_K36_K9m3_K36-K9m3/lda_outputs.scchix_inputs_clstr_by_celltype_K36-K9m3.removeNA_FALSE-merged_mat.", jmark, ".K-30.binarize.FALSE/ldaOut.scchix_inputs_clstr_by_celltype_K36-K9m3.removeNA_FALSE-merged_mat.", jmark, ".K-30.Robj")
  assertthat::assert_that(file.exists(inf.lda.tmp))
  return(inf.lda.tmp)
})

out.objs <- lapply(inf.lda.lst, function(inf.lda){
  load(inf.lda, v=T)
  return(list(out.lda = out.lda, count.mat = count.mat))
})

count.mat.lst <- lapply(out.objs, function(jout){
  jout$count.mat
})

ncuts.long <- lapply(count.mat.lst, function(jmat){
  data.frame(cell = colnames(jmat), ncuts = colSums(jmat), stringsAsFactors = FALSE)
}) %>%
  bind_rows()


# Load metas --------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"
inf.meta <- file.path(hubprefix, "jyeung/data/dblchic/gastrulation/from_analysis/metadata/from_demux_cleaned_var_filtered_manual2nocenterfilt2_K36_K9m3_K36-K9m3/demux_cleaned_filtered_var_filtered_manual2nocenterfilt2_K36_K9m3_K36-K9m3.2021-08-24.filt2.spread_7.single_and_dbl.txt")

dat.meta <- fread(inf.meta) %>%
  left_join(., ncuts.long) %>%
  rowwise() %>%
  mutate(plate = ClipLast(cell, jsep = "_"),
         cluster = ifelse(cluster == "Erythroid", "aErythroid", cluster))

ggplot(dat.meta, aes(x = umap1.shift, y = umap2.scale, color = cluster)) +
  geom_point() +
  theme_bw() + theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.meta, aes(x = umap1.shift, y = umap2.scale, color = log10(ncuts))) +
  geom_point() +
  scale_color_viridis_c() +
  theme_bw() + theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

m1 <- ggplot(dat.meta %>% filter(mark == "K36" & type == "single"), aes(x = cluster, y = ncuts)) +
  geom_point() +
  geom_boxplot() +
  scale_y_log10 () +
  facet_wrap(~plate) +
  theme_bw() +
  ggtitle("K36") +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

plates.eryth <- c("E10p5-B6C-K9m3-190409-2", "E9p5-CB6-K9m3-190409-1", "E9p5-CB6-K9m3-190409-2")
m2 <- ggplot(dat.meta %>% filter(mark == "K9m3" & type == "single" & !plate %in% plates.eryth), aes(x = cluster, y = ncuts)) +
  geom_point() +
  geom_boxplot() +
  scale_y_log10 () +
  facet_wrap(~plate) +
  theme_bw() +
  ggtitle("K9m3") +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

JFuncs::multiplot(m1, m2, cols = 2)


