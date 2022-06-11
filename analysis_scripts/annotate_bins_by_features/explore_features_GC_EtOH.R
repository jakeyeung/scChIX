# Jake Yeung
# Date of Creation: 2021-11-16
# File: ~/projects/scChIX/analysis_scripts/annotate_bins_by_features/explore_features.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(scChIX)

# Load prob matrix --------------------------------------------------------

inf.mat <- "/Users/yeung/Dropbox/macbookpro_data/data/dblchic/annotate_bins_by_promoter_enh_repeats/prob_mats/MF_TopBins_autosomesOnly.FewerTopics2_KeepAllPlates_clstrbytopics_0_1_RemoveBadMixings.K27m3-K9m3-prob_mat.K27m3-K9m3_to_K27m3.txt"
prob.dat <- fread(inf.mat)
rnames <- prob.dat$V1
prob.mat <- as.matrix(subset(prob.dat, select = -V1))
rownames(prob.mat) <- rnames

# Load GC -----------------------------------------------------------------

inf.gc <- "/Users/yeung/Dropbox/macbookpro_data/data/dblchic/annotate_bins_by_promoter_enh_repeats/gcs_genomewide.RData"
load(inf.gc, v=T)  # gr.gc.dat


# Load metadata  ----------------------------------------------------------

annots.dir <- "/Users/yeung/data/dblchic/annotations_EtOH"
jmarks.annot <- c("K27", "K9", "K27_K9"); names(jmarks.annot) <- jmarks.annot

ctypes <- c("Bcells", "Granulocytes", "NKcells"); names(ctypes) <- ctypes

annots.dat <- LoadCellAnnotsEtOH.MorePlates(annots.dir, jmarks.annot) %>%
  dplyr::select(cell.orig, celltype) %>%
  filter(celltype %in% ctypes)


# Correlate GC with prob --------------------------------------------------

prob.vec <- rowMeans(prob.mat)

cnames.keep.lst <- lapply(split(annots.dat, f = annots.dat$celltype), function(x) x$cell.orig)

prob.mat.byclst <- do.call(cbind, MeanAcrossClusters(prob.mat, cnames.keep.lst = cnames.keep.lst))

prob.mat.byclst.long <- melt(prob.mat.byclst)
colnames(prob.mat.byclst.long) <- c("bin", "ctype", "prob")

# correlate with GC content


dat.merge <- prob.mat.byclst.long %>%
  left_join(., gr.gc.dat, by = c("bin" = "bname")) %>%
  rowwise() %>%
  mutate(is.high = prob > 0.5)

cbPalette <- c("#696969", "#56B4E9", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")


outpdf <- paste0("/Users/yeung/data/dblchic/annotate_bins_by_promoter_enh_repeats/plots/EtOH_GC_low_high_prob.", Sys.Date(), ".pdf")
pdf(outpdf, useDingbats = FALSE)

ggplot(dat.merge, aes(x = prob, fill = ctype)) + 
  geom_density(alpha = 0.25) + 
  geom_vline(xintercept = 0.5, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.merge, aes(x = prob, fill = ctype)) + 
  geom_density(alpha = 1) + 
  facet_wrap(~ctype) + 
  geom_vline(xintercept = 0.5, linetype = "dotted") + 
  scale_fill_manual(values = cbPalette) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

ggplot(dat.merge, aes(x = gc, y = prob)) + 
  geom_point(alpha = 0.1) + 
  facet_wrap(~ctype) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.merge, aes(x = gc, y = log(prob / (1 - prob)))) + 
  geom_point(alpha = 0.1) + 
  theme_bw() + 
  facet_wrap(~ctype) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(dat.merge, aes(x = gc, fill = is.high, fill = ctype)) + 
  geom_density(alpha = 0.25) + 
  facet_wrap(~ctype) + 
  theme_bw() + 
  scale_fill_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.merge, aes(x = is.high, y = gc, fill = ctype)) + 
  geom_boxplot() + 
  scale_fill_manual(values = cbPalette) + 
  facet_wrap(~ctype) + 
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

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)

inf.tsspretty <- "/Users/yeung/Dropbox/macbookpro_data/data/databases/refseq/MmRefseqTss.chromorenamed.50000.again.cut.pretty.bed"
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
jbins.keep <- rnames
# get centers
dat.regions <- data.frame(bin = jbins.keep, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(chromo = GetChromo(x = bin), 
         start = as.numeric(GetStart(x = bin)), 
         end = as.numeric(GetEnd(x = bin)), 
         midpt = (start + end) / 2)
jbins.midpt <- paste(dat.regions$chromo, paste(dat.regions$midpt - 1, paste(dat.regions$midpt + 1), sep = "-"), sep = ":")

bins.annot.tmp <- AnnotateCoordsFromList(coords.vec = jbins.midpt, inf.tss = inf.tsspretty, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos)

bindists <- subset(bins.annot.tmp$regions.annotated, select = c(region_coord, distanceToTSS))
bindists$region_coord_orig <- jbins.keep


# Calculate distance to nearest bin, split by low or high -------------------------------------

dat.merge.withTSS <- left_join(dat.merge, bindists, by = c("bin" = "region_coord_orig"))

ggplot(dat.merge.withTSS, aes(y = abs(distanceToTSS) + 1, x = is.high, fill = ctype)) + 
  geom_boxplot() + 
  scale_fill_manual(values = cbPalette) + 
  facet_wrap(~ctype) + 
  scale_y_log10() +
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


dev.off()





