# Jake Yeung
# Date of Creation: 2022-07-15
# File: ~/projects/scChIX/analysis_scripts/H3K9me3-H3K36me3_gastrulation_analysis/17-load_processed_data_write_outputs.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scChIX)

jmarks <- c("K36", "K9m3", "K36-K9m3"); names(jmarks) <- jmarks
jmarks.singles <- c("K36", "K9m3"); names(jmarks.singles) <- jmarks.singles


# Load processed files ----------------------------------------------------

# load cleaned countmats

indir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChIX_gastrulation/from_cluster/snakemake_runs/K36_K9m3_K36-K9m3/snakemake_inputs/countmats"

count.mats <- lapply(jmarks, function(jmark){
  readRDS(file.path(indir, paste0("countmat_var_filt.", jmark, ".rds")))
})




# Load outputs ------------------------------------------------------------

coords.dbl.annot <- readRDS("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChIX_gastrulation/from_analysis_macbook/H3K36me3_H3K9me3_downstream_objs/dbl_cell_assignment_to_cluster.2021-10-22.rds")

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

ctypes.k36 <- levels(coords.dbl.annot$cluster.k36)
ctypes.k9 <- levels(coords.dbl.annot$cluster.k9)
clst2color <- hash::hash(ctypes.k36, cbPalette[1:length(ctypes.k36)])

m.grid.reordered <- ggplot(coords.dbl.annot, aes(x = cluster.k36, y = cluster.k9, color = cluster)) +
  geom_point(alpha = 0.75, position = ggforce::position_jitternormal(sd_x = 0.08, sd_y = 0.08)) +
  theme_bw() +
  scale_color_manual(values = cbPalette) +
  theme(aspect.ratio=0.5, axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom") +
  ggtitle("Each dot is a double stained cell,\nX-Y shows the cluster pair it is assigned")
print(m.grid.reordered)


# inf.out <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChIX_gastrulation/from_cluster/snakemake_runs/K36_K9m3_K36-K9m3/snakemake_outputs/scchix_outputs_objs/scchix_inputs_clstr_by_celltype_K36-K9m3.RData"
# load(inf.out, v=T)
# fits.out <- act.repress.coord.lst
# w.lst <- sapply(fits.out, function(x) x$w)
# 
# # if louvains are now from clusters need eto rethink jcoord
# cell.vec <- names(fits.out)
# names(cell.vec) <- cell.vec
# coords.dbl <- lapply(cell.vec, function(jcell){
#   jfit <- fits.out[[jcell]]
#   jweight <- fits.out[[jcell]]$w
#   p.mat <- SoftMax(jfit$ll.mat)
#   jcoord <- which(jfit$ll.mat == max(jfit$ll.mat), arr.ind = TRUE)
#   jmax <- max(p.mat)
#   
#   # rows are active, columns are repress I THINK?
#   # TODO: assumes underscores be careful!
#   jlouv.act <- rownames(p.mat)[[jcoord[[1]]]]
#   jlouv.repress <- colnames(p.mat)[[jcoord[[2]]]]
#   
#   if (grepl("_", jlouv.act)){
#     jlouv.act <- strsplit(jlouv.act, split = "_")[[1]][[2]]
#   }
#   if (grepl("_", jlouv.repress)){
#     jlouv.repress <- strsplit(jlouv.repress, split = "_")[[1]][[2]]
#   }
#   out.dat <- data.frame(cell = jcell, louv.act = jlouv.act, louv.repress = jlouv.repress, lnprob = jmax, w = jweight, stringsAsFactors = FALSE)
#   return(out.dat)
# }) %>%
#   bind_rows()


# Annotate louv.act and louv.repress  -------------------------------------




# Load prob mat  ----------------------------------------------------------




# Load LDA projections ----------------------------------------------------

# out.lda.predict.K36 <- out.lda.predict
# count.mat.proj.K36 <- count.mat.proj
# out.lda.K36 <- out.objs$out.lda
# 
# load("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChIX_gastrulation/from_cluster/snakemake_runs/K36_K9m3_K36-K9m3/snakemake_outputs/projection_output.K9m3.RData", v=T)
# out.lda.predict.K9m3 <- out.lda.predict
# count.mat.proj.K9m3 <- count.mat.proj
# out.lda.K9m3 <- out.objs$out.lda


# Load gastru metadata annots  --------------------------------------------

# inf.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChIX_gastrulation/from_cluster/metadata/metadata_cleaned.2021-10-04.rds"
# dat.meta.annot.merged <- readRDS(inf.meta)

indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChIX_gastrulation/from_cluster/metadata"
dat.meta.annot.merged <- lapply(jmarks.singles, function(jmark){
  fread(file.path(indir.meta, paste0("metadata_cell_cluster_with_clustercol.", jmark, ".2021-12-03.txt")))
}) %>%
  bind_rows() %>%
  dplyr::select(c(-umap1, -umap2, -umap1.scale, -umap2.scale, -umap1.shift, -umap1.flip,))


# Write processed files  --------------------------------------------------


outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChIX_gastrulation/from_cluster/processed_data"
outpdf <- file.path(outdir, paste0("procesed_data_check_plots.", Sys.Date(), ".pdf"))

pdf(outpdf, useDingbats = FALSE)

# write dbl assignments
fwrite(coords.dbl.annot, file = file.path(outdir, paste0("scChIX_gastrulation_double_cell_assignment.txt")), sep = "\t")

for (jmark in jmarks){
  # write countmat
  outfcountmat <- file.path(outdir, paste0("scChIX_gastrulation_", jmark, "_countmat.txt"))
  write.table(as.matrix(count.mats[[jmark]]), file = outfcountmat, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
  cmd2 <- paste0("gzip ", outfcountmat)
  system(cmd2)
}

for (jmark in jmarks.singles){
  print(jmark)
  # write meta
  # write LDA projects
  infproj <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChIX_gastrulation/from_cluster/snakemake_runs/K36_K9m3_K36-K9m3/snakemake_outputs/projection_output.", jmark, ".RData")
  load(infproj, v=T)
  outfrdata <- file.path(outdir, paste0("scChIX_gastrulation_", jmark, "_LDA_projection_output.RData"))
  save(out.lda.predict, count.mat.proj, out.objs, file = outfrdata)
  # gzip the text files
}



outfmeta <- file.path(outdir, paste0("scChIX_gastrulation_H3K36me3_and_H3K9me3_metadata.txt"))
fwrite(dat.meta.annot.merged, file = outfmeta, sep = "\t", na = "NA")
cmd1 <- paste0("gzip ", outfmeta)
system(cmd1)

dev.off()