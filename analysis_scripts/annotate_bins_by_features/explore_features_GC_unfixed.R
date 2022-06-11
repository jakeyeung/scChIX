# Jake Yeung
# Date of Creation: 2021-11-18
# File: ~/projects/scChIX/analysis_scripts/annotate_bins_by_features/explore_features_GC_unfixed.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(scChIX)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)

# Load prob matrix --------------------------------------------------------

# inf.mat <- "/Users/yeung/Dropbox/macbookpro_data/data/dblchic/annotate_bins_by_promoter_enh_repeats/prob_mats/MF_TopBins_autosomesOnly.FewerTopics2_KeepAllPlates_clstrbytopics_0_1_RemoveBadMixings.K27m3-K9m3-prob_mat.K27m3-K9m3_to_K27m3.txt"
inf.mat <- "/Users/yeung/Dropbox/macbookpro_data/data/dblchic/annotate_bins_by_promoter_enh_repeats/prob_mats/scchix_inputs_clstr_by_celltype-prob_mat.K4m1-K27m3_to_K4m1.txt.gz"
prob.dat <- fread(inf.mat)
rnames <- prob.dat$V1
prob.mat <- as.matrix(subset(prob.dat, select = -V1))
rownames(prob.mat) <- rnames

# write bins
# saveRDS(rnames, paste0("/Users/yeung/data/dblchic/unfixed_revisions/metadata_5kb/rownames_5kb_coords.", Sys.Date(), ".rds"))

# Load GC -----------------------------------------------------------------

# inf.gc <- "/Users/yeung/Dropbox/macbookpro_data/data/dblchic/annotate_bins_by_promoter_enh_repeats/gcs_genomewide.RData"
inf.gc <- "/Users/yeung/Dropbox/macbookpro_data/data/dblchic/unfixed_revisions/metadata_5kb/gcs_genomewide_5kb.2021-11-19.RData"
load(inf.gc, v=T)  # gr.gc.dat


# Load metadata  ----------------------------------------------------------

inf.meta <- "/Users/yeung/Dropbox/macbookpro_data/data/dblchic/unfixed_revisions/metadata_original/BM_UnfixedTopics.FinalCellClusterTable.NAimputed.2020-04-28.RData"
load(inf.meta, v=T)


cnames.orig <- c("BcellsNaive", "BcellsMemory", "Basophils", "Eosinophils", "InnateLymph", "BcellsPlasma")
cnames.new <- c("Bcells", "ProB", "Baso/Eosino", "Baso/Eosino", "Innate Lymph", "ProB")
cnames.hash <- hash::hash(cnames.orig, cnames.new)


inf.meta.withcolor <- "/Users/yeung/Dropbox/macbookpro_data/data/dblchic/unfixed_revisions/metadata_original/metadata_K4me1_K27me3_UMAPs_merged.remove_bad_cluster.rename_Bcells.2021-03-19.txt"
dat.meta.withcolor <- fread(inf.meta.withcolor) %>%
  rowwise() %>%
  mutate(cluster.renamed = AssignHash(x = cluster.knn, jhash = cnames.hash, null.fill = cluster.knn))

cluster2color <- hash::hash(dat.meta.withcolor$cluster.renamed, dat.meta.withcolor$color)


dat.final.annots <- dat.final.annots %>%
  rowwise() %>%
  dplyr::mutate(ctype  = cluster,
                ctype.new = AssignHash(x = ctype, jhash = cnames.hash, null.fill = ctype),
                clustercol = cluster2color[[ctype.new]])

dat.clst.toappend <- subset(dat.final.annots, select = c(ctype, ctype.new, clustercol))
dat.clst.toappend <- dat.clst.toappend[!duplicated(dat.clst.toappend), ]

# Correlate GC with prob --------------------------------------------------

# prob.vec <- rowMeans(prob.mat)

cnames.keep.lst <- lapply(split(dat.final.annots, f = dat.final.annots$cluster), function(x) x$cell)
prob.mat.byclst <- do.call(cbind, MeanAcrossClusters(prob.mat, cnames.keep.lst = cnames.keep.lst))

prob.mat.byclst.long <- melt(prob.mat.byclst)
colnames(prob.mat.byclst.long) <- c("bin", "ctype", "prob")

prob.mat.byclst.long <- prob.mat.byclst.long %>%
  left_join(., dat.clst.toappend)

# prob.mat.byclst.long <- prob.mat.byclst.long 
  # rowwise() %>%
  # mutate(clustercol = cluster2color[[as.character(ctype)]])



# correlate with GC content

ggplot(prob.mat.byclst.long, aes(x = prob, fill = clustercol)) + 
  geom_density() + 
  facet_wrap(~ctype) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# celltype-specific threhsold: 

clsts <- unique(dat.final.annots$cluster)
thres.vec <- sapply(clsts, function(x) ifelse(x == "Erythroblasts", 0.75, 0.5))
thres.hash <- hash::hash(clsts, thres.vec)

dat.merge <- prob.mat.byclst.long %>%
  left_join(., gr.gc.dat, by = c("bin" = "bname")) %>%
  # rowwise() %>%
  mutate(is.high = ifelse(ctype == "Erythroblasts", prob < 0.75, prob < 0.5))
  # mutate(is.high = prob > thres.hash[[as.character(ctype)]])

dat.cutoffs <- data.frame(ctype.new = unique(dat.merge$ctype.new), stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(cutoff = ifelse(ctype.new == "Erythroblasts", 0.75, 0.5))

library(ggrastr)

outpdf <- paste0("/Users/yeung/data/dblchic/annotate_bins_by_promoter_enh_repeats/plots/unfixed_GC_low_high_prob.", Sys.Date(), ".pdf")
pdf(outpdf, useDingbats = FALSE)

ggplot(dat.merge, aes(x = prob, fill = clustercol)) + 
  geom_density(alpha = 0.25) + 
  scale_fill_identity() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.merge, aes(x = prob, fill = clustercol)) + 
  geom_density(alpha = 1) + 
  scale_fill_identity() + 
  geom_vline(mapping = aes(xintercept = cutoff), data = dat.cutoffs, linetype = "dotted") + 
  facet_wrap(~ctype.new) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# ggplot(dat.merge, aes(x = gc, y = prob)) + 
#   geom_point_rast(alpha = 0.1) + 
#   facet_wrap(~ctype) + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# ggplot(dat.merge, aes(x = gc, y = log(prob / (1 - prob)))) + 
#   geom_point_rast(alpha = 0.1) + 
#   theme_bw() + 
#   facet_wrap(~ctype) + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(dat.merge, aes(x = gc, fill = is.high)) + 
  geom_density(alpha = 0.25) + 
  facet_wrap(~ctype.new) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.merge, aes(x = is.high, y = gc, fill = clustercol)) + 
  geom_boxplot() + 
  scale_fill_identity() + 
  facet_wrap(~ctype.new) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

jfit0 <- lm(formula = gc~ 1, data = dat.merge)
jfit1 <- lm(formula = gc~ is.high, data = dat.merge)
jfit2 <- lm(formula = gc ~ is.high + ctype, data = dat.merge)
jfit3 <- lm(formula = gc ~ is.high + ctype + is.high : ctype, data = dat.merge)

anova(jfit0, jfit1)
anova(jfit1, jfit2)
anova(jfit1, jfit3)

# Show distance between genes  --------------------------------------------


inf.tsspretty <- "/Users/yeung/Dropbox/macbookpro_data/data/databases/refseq/MmRefseqTss.chromorenamed.50000.again.cut.pretty.bed"
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
jbins.keep <- rnames
# get centers
dat.regions <- data.frame(bin = jbins.keep, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(chromo = GetChromo(x = bin), 
         start = as.numeric(GetStart(x = bin)), 
         end = as.numeric(GetEnd(x = bin)), 
         midpt = (start + end) / 2,
         binmid = paste(chromo, paste(midpt - 1, midpt + 1, sep = "-"), sep = ":"))

# binshash <- hash::hash(dat.regions$bin, dat.regions$binmid)

jbins.midpt <- paste(dat.regions$chromo, paste(dat.regions$midpt - 1, paste(dat.regions$midpt + 1), sep = "-"), sep = ":")

bins.annot.tmp <- AnnotateCoordsFromList(coords.vec = jbins.midpt, inf.tss = inf.tsspretty, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos)

bindists <- subset(bins.annot.tmp$regions.annotated, select = c(region_coord, distanceToTSS)) %>%
  left_join(., subset(dat.regions, select = c(bin, binmid)), by = c("region_coord" = "binmid"))

bindists$region_coord_orig <- bindists$bin


# Calculate distance to nearest bin, split by low or high -------------------------------------

dat.merge.withTSS <- left_join(dat.merge, bindists, by = c("bin" = "region_coord_orig"))

ggplot(dat.merge.withTSS, aes(y = abs(distanceToTSS) + 1, x = is.high, fill = clustercol)) + 
  geom_boxplot() + 
  facet_wrap(~ctype.new) + 
  scale_fill_identity() + 
  scale_y_log10() +
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dev.off()


outf <- paste0("/Users/yeung/Dropbox/macbookpro_data/data/dblchic/annotate_bins_by_promoter_enh_repeats/objs/dat_merge_withTSS.", Sys.Date(), ".rds")
saveRDS(dat.merge.withTSS, file = outf)


