# Jake Yeung
# Date of Creation: 2021-10-21
# File: ~/projects/scChIX/analysis_scripts/H3K9me3-H3K36me3_gastrulation_analysis/6-make_QC_plots.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(topicmodels)

jmarks <- c("K36", "K9m3", "K36-K9m3"); names(jmarks) <- jmarks

indir <- "/Users/yeung/Dropbox/macbookpro_data/data/dblchic"
bsize <- "50000"

outdir <- "/Users/yeung/data/dblchic/gastrulation/primetime_plots"

outpdf <- file.path(outdir, paste0("qc_plots_total_cuts_ta_frac_intrachrom_var.", Sys.Date(), ".pdf"))

pdf(outpdf, useDingbats = FALSE)

# Load raw data -----------------------------------------------------------

infs.raw <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.tmp <- file.path(indir, paste0("from_hub/countmats/", jmark, ".countTable.binsize_", bsize, ".rds"))
  assertthat::assert_that(file.exists(inf.tmp))
  return(inf.tmp)
})

infs.rz <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.tmp <- file.path(indir, paste0("from_hub/RZ_counts_dedup/gastru_merged_", jmark, ".rowsdeduped.tagged.sorted.LH_counts.csv"))
  assertthat::assert_that(file.exists(inf.tmp))
  return(inf.tmp)
})

dats.raw <- lapply(infs.raw, function(jinf){
  readRDS(jinf)
})

dats.rz <- lapply(infs.rz, function(jinf){
  ReadLH.SummarizeTA(jinf)
})

dats.sums <- lapply(dats.raw, function(jdat){
  data.frame(cell = colnames(jdat), counts = Matrix::colSums(jdat), stringsAsFactors = FALSE)
})

infs.intra <- lapply(jmarks, function(jmark){
  inf.tmp <- file.path(indir, paste0("from_hub/intrachromvar_objs/K36_K9m3_K36-K9m3/intrachrom_var_outputs.", jmark, ".2021-07-15.rds"))
  assertthat::assert_that(file.exists(inf.tmp))
  return(inf.tmp)
})

dats.intrachrom <- lapply(jmarks, function(jmark){
  jinf <- infs.intra[[jmark]]
  readRDS(jinf) %>%
    left_join(., dats.rz[[jmark]], by = c("cell" = "samp"))
})

infs.lda <- lapply(jmarks, function(jmark){
  inf.tmp <- file.path(indir, paste0("from_hub/lda_objs/filtered_count_tables_filtered_counts_varcutoffmin_0.03.noLSIfilt/lda_output_filt.", jmark, ".2021-07-14.rds"))
  assertthat::assert_that(file.exists(inf.tmp))
  return(inf.tmp)
})

lda.lst <- lapply(infs.lda, function(jinf){
  readRDS(jinf)
})

tm.lst <- lapply(lda.lst, function(jlda){
  posterior(jlda)
})

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
dat.var.lst <- lapply(jmarks, function(jmark){
  tm <- tm.lst[[jmark]]
  dat.impute.log <- log2(t(tm$topics %*% tm$terms))
  dat.var <- CalculateVarAll(dat.impute.log, jchromos) %>%
    left_join(., dats.rz[[jmark]], by = c("cell" = "samp"))
})


# Load meta data ----------------------------------------------------------

inf.meta <- "/Users/yeung/Dropbox/macbookpro_data/data/dblchic/gastrulation/snakemake_downstream_outputs/metadata_cleaned.2021-10-04.rds"
dat.meta <- readRDS(inf.meta)
cells.good <- unique(dat.meta$cell)
plates.good <- unique(sapply(dat.meta$cell, function(x) ClipLast(x, jsep = "_")))


# Plot --------------------------------------------------------------------

dats.rz.annot <- lapply(dats.rz, function(jdat){
  jdat <- jdat %>%
    rowwise() %>%
    mutate(plate = ClipLast(samp, jsep = "_")) %>%
    left_join(., dat.meta, by = c("samp" = "cell")) %>%
    rowwise() %>%
    mutate(is.good = !is.na(cluster)) 
})

jmark <- jmarks[[1]]

lapply(jmarks, function(jmark){
  ggplot(dats.rz.annot[[jmark]] %>% filter(plate %in% plates.good), aes(x = log10(total.count), y = TA.frac, color = is.good)) + 
    geom_point(alpha = 0.5) + 
    theme_bw() + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})

ggplot(dats.rz.annot[[jmark]], aes(x = log10(total.count))) + 
  geom_histogram() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dats.rz.annot[[jmark]], aes(x = log10(total.count), fill = is.good)) + 
  geom_histogram(position = "dodge") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Do cutoffs histograms across colors -------------------------------------


dats.rz.annot.long <- lapply(jmarks, function(jmark){
  dats.rz.annot[[jmark]] %>%
    mutate(mark = jmark)
}) %>%
  bind_rows()

dats.rz.annot.long$mark <- factor(dats.rz.annot.long$mark, levels = jmarks)

ggplot(dats.rz.annot.long, aes(x = log10(total.count), fill = mark)) + 
  geom_density() + 
  facet_wrap(~mark) + 
  geom_vline(xintercept = 3, linetype = "dotted") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

ggplot(dats.rz.annot.long, aes(x = TA.frac, fill = mark)) + 
  geom_density() + 
  facet_wrap(~mark) + 
  geom_vline(xintercept = 0.5, linetype = "dotted") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

# Calculate bins across all cells ----------------------------------------------------------

bin.avgs.lst <- lapply(dats.raw, function(jdat){
  log10(rowMeans(jdat) * 1000 + 1)
})

# define peaks arbitrarily
bin.cutoff.lst <- lapply(bin.avgs.lst, function(bin.avgs){
  bin.cutoff <- mean(bin.avgs)
})

peaks.lst <- lapply(jmarks, function(jmark){
  bin.avgs <- bin.avgs.lst[[jmark]]
  bin.cutoff <- bin.cutoff.lst[[jmark]]
  plot(density(bin.avgs))
  abline(v = bin.cutoff)
  peaks <- names(bin.avgs)[bin.avgs > bin.cutoff]
  print(length(peaks) / length(bin.avgs))
  return(peaks)
})

# Calculate fraction of cells in or out of peaks  -------------------------

cell.sums.lst <- lapply(dats.raw, function(jmat){
  Matrix::colSums(jmat)
})

cell.sums.atpeaks.lst <- lapply(jmarks, function(jmark){
  jmat <- dats.raw[[jmark]]
  peaks <- peaks.lst[[jmark]]
  Matrix::colSums(jmat[peaks, ])
})

cell.data.lst <- lapply(jmarks, function(jmark){
  jmat <- dats.raw[[jmark]]
  data.frame(cell = colnames(jmat), total.count.from.mat = cell.sums.lst[[jmark]], count.at.peak = cell.sums.atpeaks.lst[[jmark]], 
             stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(frac.count.in.peak = count.at.peak / total.count.from.mat) %>%
    left_join(., dats.rz.annot[[jmark]], by = c("cell" = "samp"))
})


m.lst <- lapply(jmarks, function(jmark){
  # ggplot(cell.data.lst[[jmark]], aes(x = frac.count.in.peak, fill = is.good)) + 
  ggplot(cell.data.lst[[jmark]] %>% filter(plate %in% plates.good), aes(x = frac.count.in.peak, fill = is.good)) + 
    geom_density(alpha = 0.35) + 
    ggtitle(jmark) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})
print(m.lst)

m.nofill.lst <- lapply(jmarks, function(jmark){
  # ggplot(cell.data.lst[[jmark]], aes(x = frac.count.in.peak, fill = is.good)) + 
  ggplot(cell.data.lst[[jmark]] %>% filter(plate %in% plates.good), aes(x = frac.count.in.peak)) + 
    geom_density(alpha = 0.35) + 
    ggtitle(jmark) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})
print(m.nofill.lst)

m.points.lst <- lapply(jmarks, function(jmark){
  ggplot(cell.data.lst[[jmark]] %>% filter(plate %in% plates.good), aes(x = frac.count.in.peak, y = log10(total.count), color = is.good)) + 
    geom_point() + 
    ggtitle(jmark) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})
print(m.points.lst)

# 
# rz.annot.filt2 <- left_join(rz.annot.filt1, cell.data) %>%
#   filter(frac.count.in.peak >= frip.cutoff)
# 
# ggplot(cell.data, aes(x = frac.count.in.peak)) +
#   geom_density(alpha = 0.5) +
#   theme_bw() +
#   geom_vline(xintercept = frip.cutoff, linetype = "dotted") +
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 



# Load intrachromosomal variance ------------------------------------------

jquant <- 1

log.jcutoff.stringent.lst <- c(0.2, -0.5, -0.7) + 0.1
names(log.jcutoff.stringent.lst) <- jmarks

jcutoff.stringent.lst <- c(1.2, 0.7, 0.5)
names(jcutoff.stringent.lst) <- jmarks

dat.cutoff <- data.frame(mark = names(jcutoff.stringent.lst), cutoff = jcutoff.stringent.lst, stringsAsFactors = FALSE)
dat.cutoff$mark <- factor(dat.cutoff$mark, levels = jmarks)

dat.intravar.long <- lapply(jmarks, function(jmark){
  jsub1 <- dat.var.lst[[jmark]] %>%
    ungroup() %>%
    mutate(cell.var.within.sum.norm.capped = ifelse(cell.var.within.sum.norm > quantile(cell.var.within.sum.norm, jquant), quantile(cell.var.within.sum.norm, jquant), cell.var.within.sum.norm)) %>%
    left_join(., subset(dats.rz.annot[[jmark]], select = c(samp, is.good)), by = c("cell" = "samp")) %>%
    mutate(mark = jmark)
}) %>%
  bind_rows()
dat.intravar.long$mark <- factor(dat.intravar.long$mark, levels = jmarks)

m.intravar <- ggplot(dat.intravar.long, aes(x = cell.var.within.sum.norm, fill = mark)) + 
  geom_density() + 
  scale_x_log10() + 
  geom_vline(mapping = aes(xintercept = cutoff), data = dat.cutoff, linetype = "dotted") + 
  facet_wrap(~mark) + 
  # facet_grid(is.good ~ mark) + 
  xlab("Intrachromosomal Variance") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
print(m.intravar)

dev.off()