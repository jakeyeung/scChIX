# Jake Yeung
# Date of Creation: 2021-12-02
# File: ~/projects/scChIX/analysis_scripts/H3K9me3-H3K36me3_gastrulation_analysis/10-get_marker_genes.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(Seurat)
library(parallel)

jmarks <- c("K36", "K9m3"); names(jmarks) <- jmarks


# Load count tables -------------------------------------------------------

inf.counts <- lapply(jmarks, function(jmark){
  inf.tmp <- paste0("/Users/yeung/Dropbox/macbookpro_data/data/dblchic/gastrulation/snakemake_downstream_outputs/scchix_unmixing_downstream/scchix_inputs_clstr_by_celltype-merged_mat.", jmark, ".rds")
  assertthat::assert_that(file.exists(inf.tmp))
  return(inf.tmp)
})

count.mat.lst <- lapply(inf.counts, function(jinf){
  readRDS(jinf)
})

# load metadat
inf.meta <- "/Users/yeung/Dropbox/macbookpro_data/data/dblchic/gastrulation/snakemake_downstream_outputs/metadata_cleaned.2021-10-04.rds"
dat.meta.all <- readRDS(inf.meta) %>%
  mutate(cluster = as.character(cluster))

blood.ctypes <- c("WhiteBloodCells", "Erythroid")

# for K9me3 combine non-blood cells
dat.meta.k36 <- subset(dat.meta.all, mark == "K36") %>%
  as.data.frame() 
dat.meta.k9 <- subset(dat.meta.all, mark == "K9m3") %>%
  mutate(cluster = ifelse(cluster %in% c(blood.ctypes), cluster, "NonBlood")) %>%
  as.data.frame()

rownames(dat.meta.k36) <- dat.meta.k36$cell
rownames(dat.meta.k9) <- dat.meta.k9$cell

dat.meta.lst <- list(dat.meta.k36, dat.meta.k9); names(dat.meta.lst) <- jmarks

# Run marker genes --------------------------------------------------------


# Load metadata -----------------------------------------------------------

objs.lst <- lapply(jmarks, function(jmark){
  count.mat <- count.mat.lst[[jmark]]
  dat.meta <- dat.meta.lst[[jmark]]
  seuratobj <- CreateSeuratObject(counts = count.mat, meta.data = dat.meta.lst[[jmark]])
  Idents(object = seuratobj) <- seuratobj@meta.data$'cluster'
  return(seuratobj)
})


outdir <- "/Users/yeung/data/dblchic/gastrulation/H3K36me3_H3K9me3_celltyping"
dir.create(outdir)

# marker genes
maxclsts <- 10

for (jmark in jmarks){
  
  seuratobj <- objs.lst[[jmark]]
  
  clsts <- sort(as.character(unique(Idents(object = seuratobj)))); names(clsts) <- clsts
  
  print(clsts)
  
  print("Getting marker genes")
  print(jmark)
  outrds <- file.path(outdir, paste0("marker_genes_output.", jmark, ".", Sys.Date(), ".rds"))
  outtxtbase <- file.path(outdir, paste0("marker_genes_textoutput.", jmark, ".", Sys.Date()))
  
  print(length(clsts))
  ncores <- ifelse(length(clsts) > maxclsts, maxclsts, length(clsts))
  print(ncores)
  
  markers.lst <- mclapply(clsts, function(jclst){
    FindMarkers(object = seuratobj, ident.1 = jclst, ident.2 = NULL, test.use = "poisson")
  }, mc.cores = ncores)
  
  saveRDS(markers.lst, file = outrds)
  
  for (jclst in clsts){
    print(jclst)
    write.table(markers.lst[[jclst]], file = paste0(outtxtbase, "_", jclst, ".txt"), sep = "\t", row.names = TRUE, col.names = NA)
  }
}

print("Done")


