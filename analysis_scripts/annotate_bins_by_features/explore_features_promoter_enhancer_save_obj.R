# Jake Yeung
# Date of Creation: 2021-11-17
# File: ~/projects/scChIX/analysis_scripts/annotate_bins_by_features/explore_features_promoter_enhancer.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

CountFeaturesTotal <- function(inf.mat, dat.rz.long, label){
  
  dat.mat <- ReadMatTSSFormat(inf.mat)
  dat.mat.annot <- subset(fread(inf.mat, skip = 1), select = c("reference_name", "start", "end", "bname")) 
  
  dat.mat.annot$coord <- paste(dat.mat.annot$reference_name, paste(dat.mat.annot$start, dat.mat.annot$end, sep = "-"), sep = ":")
  
  rownames(dat.mat) <- dat.mat.annot$coord
  counts.pbulk <- data.frame(region = rownames(dat.mat), cuts = Matrix::rowSums(dat.mat), stringsAsFactors = FALSE)
  
  jsum <- table(dat.mat.annot$bname)
  jsum.filt <- jsum[which(jsum > 2)]
  print(jsum.filt)
  jsum.filt.comp <- jsum[which(jsum <= 2)]
  
  genes.all <- names(jsum)[which(jsum <= 2)]
  
  proms.i <- which(dat.mat.annot$bname %in% genes.all)
  enhs.i <- which(startsWith(x = dat.mat.annot$bname, prefix = "enhancer"))
  repet.i <- which(startsWith(x = dat.mat.annot$bname, prefix = "repetitive"))
  trans.i <- which(startsWith(x = dat.mat.annot$bname, prefix = "transposable"))
  
  # merge genes and call it "promoter" features
  # types <- c("proms", "enhs", "repet", "trans"); names(types) <- types
  regions.i.lst <- list("proms" = proms.i, 
                        "enhs" = enhs.i, 
                        "repet" = repet.i, 
                        "trans" = trans.i)
  
  cnames.keep.lst <- lapply(regions.i.lst, function(i.vec){
    counts.pbulk$region[i.vec]
  })
  
  dat.sum.by.regions <- SumAcrossClusters(t(dat.mat), cnames.keep.lst = cnames.keep.lst)
  dat.sum.by.regions <- as.data.frame(do.call(cbind, dat.sum.by.regions))
  dat.sum.by.regions$samp <- rownames(dat.sum.by.regions)
  
  
  dat.sum.by.regions.annot <- left_join(dat.sum.by.regions, dat.rz.long, by = "samp")
  
  dat.sum.by.regions.annot.long <- dat.sum.by.regions.annot %>%
    dplyr::select(proms, enhs, repet, trans, samp) %>%
    melt() %>%
    left_join(., subset(dat.rz.long, select = c(samp, total.count)))
  dat.sum.by.regions.annot.long$label <- label
  return(dat.sum.by.regions.annot.long)
}

# Load meta ---------------------------------------------------------------

inf.meta <- "/Users/yeung/Dropbox/macbookpro_data/data/dblchic/unfixed_revisions/metadata_original/BM_UnfixedTopics.FinalCellClusterTable.NAimputed.2020-04-28.RData"
load(inf.meta, v=T)  # dat.final.annots



# Get sum  ----------------------------------------------------------------

jmarks <- c("K4m1", "K27m3", "K4m1_K27m3"); names(jmarks) <- jmarks
indir <- "/Users/yeung/data/dblchic/qc_unfixed_2021/RZdat.all.deeper"

infs.rz <- lapply(jmarks, function(jmark){
  inf.rz.tmp <- file.path(indir, paste0("all_BM_", jmark, "_200119.LHcounts.csv.gz"))
  assertthat::assert_that(file.exists(inf.rz.tmp))
  return(inf.rz.tmp)
})


dat.rz.long <- lapply(jmarks, function(jmark){
  jinf <- infs.rz[[jmark]]
  print(jinf)
  dat.rz <- ReadLH.SummarizeTA(jinf, remove.nones = FALSE, na.to.zero = TRUE) %>%
    rowwise() %>%
    mutate(experi = ClipLast(samp),
           cellindx = paste("cell", strsplit(samp, "_")[[1]][[2]], sep = ""),
           mark = jmark)
  
}) %>%
  bind_rows()

# dat.sum.by.regions.annot <- left_join(dat.sum.by.regions, dat.rz.long, by = "samp")
# dat.sum.by.regions.annot.long <- dat.sum.by.regions.annot %>%
#   dplyr::select(proms, enhs, repet, trans, samp) %>%
#   melt() %>%
#   left_join(., subset(dat.rz.long, select = c(samp, total.count)))



# Load infs ---------------------------------------------------------------


indir <- "/Users/yeung/data/dblchic/annotate_bins_by_promoter_enh_repeats/unfixed"
fnames <- list.files(indir, pattern = "all.*.csv.gz", all.files = TRUE, full.names = FALSE)
ctypes <- sapply(fnames, function(x) strsplit(x, "\\.")[[1]][[2]], USE.NAMES = FALSE); names(ctypes) <- ctypes

jmarks <- c("K4m1", "K27m3"); names(jmarks) <- jmarks
jmarkctypes <- as.character(interaction(jmarks, ctypes, sep = "_")); names(jmarkctypes) <- jmarkctypes

infs <- lapply(jmarkctypes, function(jmarkctype){
  jmark <- strsplit(jmarkctype, split = "_")[[1]][[1]]
  jctype <- strsplit(jmarkctype, split = "_")[[1]][[2]]
  inf.tmp <- file.path(indir, paste0("all_BM_", jmark, "_200119.", jctype, ".sorted.countTable.features.mapq_40.TSSdist_1000.blfiltered.csv.gz"))
  assertthat::assert_that(file.exists(inf.tmp))
  return(inf.tmp)
})


dat.sum.by.regions.annot.long.lst <- lapply(jmarkctypes, function(jmarkctype){
  print(jmarkctype)
  CountFeaturesTotal(inf.mat = infs[[jmarkctype]], dat.rz.long = dat.rz.long, label = jmarkctype)
})


outrds <- "/Users/yeung/data/dblchic/annotate_bins_by_promoter_enh_repeats/dat_sum_by_regions.rds"
saveRDS(dat.sum.by.regions.annot.long.lst, file = outrds)




# 
# 
# 
# # Load mat ----------------------------------------------------------------
# 
# inf.mat <- "/Users/yeung/Dropbox/macbookpro_data/data/dblchic/annotate_bins_by_promoter_enh_repeats/unfixed/all_BM_K27m3_200119.DC.sorted.countTable.features.mapq_40.TSSdist_1000.blfiltered.csv.gz"
# dat.mat <- ReadMatTSSFormat(inf.mat)
# dat.mat.annot <- subset(fread(inf.mat, skip = 1), select = c("reference_name", "start", "end", "bname")) 
#  
# dat.mat.annot$coord <- paste(dat.mat.annot$reference_name, paste(dat.mat.annot$start, dat.mat.annot$end, sep = "-"), sep = ":")
# 
# rownames(dat.mat) <- dat.mat.annot$coord
# counts.pbulk <- data.frame(region = rownames(dat.mat), cuts = Matrix::rowSums(dat.mat), stringsAsFactors = FALSE)
# 
# jsum <- table(dat.mat.annot$bname)
# jsum.filt <- jsum[which(jsum > 2)]
# print(jsum.filt)
# jsum.filt.comp <- jsum[which(jsum <= 2)]
# 
# genes.all <- names(jsum)[which(jsum <= 2)]
# 
# proms.i <- which(dat.mat.annot$bname %in% genes.all)
# enhs.i <- which(startsWith(x = dat.mat.annot$bname, prefix = "enhancer"))
# repet.i <- which(startsWith(x = dat.mat.annot$bname, prefix = "repetitive"))
# trans.i <- which(startsWith(x = dat.mat.annot$bname, prefix = "transposable"))
# 
# # merge genes and call it "promoter" features
# # types <- c("proms", "enhs", "repet", "trans"); names(types) <- types
# regions.i.lst <- list("proms" = proms.i, 
#                       "enhs" = enhs.i, 
#                       "repet" = repet.i, 
#                       "trans" = trans.i)
# 
# cnames.keep.lst <- lapply(regions.i.lst, function(i.vec){
#   counts.pbulk$region[i.vec]
# })
# 
# dat.sum.by.regions <- SumAcrossClusters(t(dat.mat), cnames.keep.lst = cnames.keep.lst)
# dat.sum.by.regions <- as.data.frame(do.call(cbind, dat.sum.by.regions))
# dat.sum.by.regions$samp <- rownames(dat.sum.by.regions)
# 
# 
# 
# 
# 
# ggplot(dat.sum.by.regions.annot.long %>% filter(variable %in% c("proms", "enhs")), aes(x = variable, y = value / total.count)) + 
#   geom_boxplot() + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# 
# 
# 
# 
#  