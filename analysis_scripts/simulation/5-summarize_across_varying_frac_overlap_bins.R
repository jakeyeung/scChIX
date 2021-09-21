# Jake Yeung
# Date of Creation: 2021-09-20
# File: ~/projects/scChIX/analysis_scripts/simulation/5-summarize_across_varying_frac_overlap_bins.R
#

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)



# Load downstream objs ----------------------------------------------------


hubprefix <- "/home/jyeung/hub_oudenaarden"
cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
cbPalette.ctype <- c("#FFB6C1", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
frac.mutexcl.str.vec <- c("0.01", "0.5", "1.0"); names(frac.mutexcl.str.vec) <- frac.mutexcl.str.vec
jmarks <- c("mark1", "mark2", "mark1-mark2"); names(jmarks) <- jmarks


outpdf <- file.path(hubprefix, paste0("jyeung/data/dblchic/from_cluster/simulation_data_downstream/plots/scchix_downstream_inference_checks_across_frac_overlapping_bins.", Sys.Date(), ".pdf"))

pdf(outpdf, useDingbats = FALSE)

dat.umap.lst <- lapply(frac.mutexcl.str.vec, function(jfrac){
  inf <- file.path(hubprefix, paste0("jyeung/data/dblchic/from_cluster/simulation_data_downstream/plots/dat_umap_joined.downstream_overlaps.frac_mutexcl_", jfrac, ".rds"))
  dat.umap.tmp <- readRDS(inf) %>%
    rowwise() %>%
    mutate(frac.mutexcl.str = jfrac,
           frac.mutexcl.recalc = as.numeric(ifelse(jfrac == "1.0", "0.99", jfrac)),
           frac.overlapping.bins = 1 - frac.mutexcl.recalc)
})

dat.binmeans.lst.lst <- lapply(frac.mutexcl.str.vec, function(jfrac){
  inf <- file.path(hubprefix, paste0("jyeung/data/dblchic/from_cluster/simulation_data_downstream/plots/bin_means.downstream_overlaps.frac_mutexcl_", jfrac, ".rds"))
  dat.binmeans.lst <- lapply(readRDS(inf), function(jdat){
    jdat$frac.mutexcl.str <- jfrac
    jdat$frac.mutexcl.recalc <- as.numeric(ifelse(jfrac == "1.0", "0.99", jfrac))
    jdat$frac.overlapping.bins <- 1 - jdat$frac.mutexcl.recalc
    return(jdat)
  })
})


# Calculate standard deviation vs mean -------------------------------------------

dat.binmeans.binned <- dat.binmeans.lst.lst$`1.0` %>%
  bind_rows() %>%
  rowwise() %>%
  mutate(xbin = round(BinMean.dbl, 2)) %>%
  group_by(xbin) %>%
  summarise(Mean = mean(overlap.estimate),
            StdDev = sd(overlap.estimate),
            CI.lower = quantile(overlap.estimate, 0.05),
            CI.upper = quantile(overlap.estimate, 0.95),
            CI.lower.centered = CI.lower - Mean,
            CI.upper.centered = CI.upper - Mean)

ggplot(dat.binmeans.binned, aes(x = xbin, y = StdDev)) +
  geom_point() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Fraction of Bin Signal Multiplexed to Mark 1") +
  ylab("Standard Deviation")

ggplot(dat.binmeans.binned, aes(x = xbin, y = Mean, ymin = CI.lower, ymax = CI.upper)) +
  geom_abline(slope = 1, color = 'blue', size = 1, alpha = 0.5) +
  geom_point() +
  geom_errorbar() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("True Fraction of Bin Signal Multiplexed to Mark 1") +
  ylab("Inferred Fraction of Bin Signal Multiplexed to Mark 1\n(ErrorBars: 95% CI)")

ggplot(dat.binmeans.binned, aes(x = xbin, y = 0, ymin = CI.lower.centered, ymax = CI.upper.centered)) +
  geom_point() +
  geom_errorbar() +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("True Fraction of Bin Signal Multiplexed to Mark 1") +
  ylab("95% Confidence Interval (+/- Fraction)") +
  coord_cartesian(ylim = c(-0.1, 0.1))

ggplot(dat.binmeans.binned, aes(x = xbin, y = 0, ymin = CI.lower.centered, ymax = CI.upper.centered)) +
  geom_point() +
  geom_errorbar() +
  theme_bw() +
  theme(aspect.ratio=0.33, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("True Fraction of Bin Signal Multiplexed to Mark 1") +
  ylab("95% Confidence Interval (+/- Fraction)") +
  coord_cartesian(ylim = c(-0.1, 0.1))


# Plot umaps  -------------------------------------------------------------



ggplot(dat.umap.lst %>% bind_rows(), aes(x = umap1.scale.shift, y = umap2.scale, color = ctype, group = cell)) +
  geom_point() +
  geom_path(alpha = 0.02) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_manual(values = cbPalette.ctype, na.value = "grey85") +
  facet_wrap(~frac.overlapping.bins) +
  theme_bw() +
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")



# Plot demux prob ---------------------------------------------------------

dat.binmeans.lst <- lapply(frac.mutexcl.str.vec, function(jfrac){
  dat.binmeans.lst.lst[[jfrac]] %>%
    bind_rows()
})

# plot true signal
# m.overlaps <- ggplot(bin.means.merged.lst %>% bind_rows(), aes(x = BinMean.dbl, fill = annot, y =..count../sum(..count..))) +
#     geom_histogram(bins = 50, position = "identity", alpha = 0.5) +
#     theme_bw() +
#     xlab("Fraction of signal belonging to mark1") +
#     ylab("Fraction of bins") +
#     theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

m.lst <- lapply(frac.mutexcl.str.vec, function(jfrac){
  jtitle <- ifelse(jfrac == "1.0", "0.99", jfrac)
  jtitle <- 1 - as.numeric(jtitle)
  m <- ggplot(dat.binmeans.lst[[jfrac]], aes(x = BinMean.dbl, fill = annot, y = ..count../sum(..count..))) +
    geom_histogram(bins = 50, position = "identity", alpha = 0.5) +
    ggtitle(paste("True Fraction of Signal Multiplexed to Mark 1. Frac overlap: ", jtitle)) +
    theme_bw() +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  return(m)
})

print(m.lst)
JFuncs::multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], cols = 3)


m.infer.lst <- lapply(frac.mutexcl.str.vec, function(jfrac){
  jtitle <- ifelse(jfrac == "1.0", "0.99", jfrac)
  jtitle <- 1 - as.numeric(jtitle)
  m <- ggplot(dat.binmeans.lst[[jfrac]], aes(x = overlap.estimate, fill = annot, y = ..count../sum(..count..))) +
    geom_histogram(bins = 50, position = "identity", alpha = 0.5) +
    ggtitle(paste("Inferred Fraction of Signal Multiplexed to Mark 1. Frac overlap: ", jtitle)) +
    theme_bw() +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  return(m)
})
print(m.infer.lst)
JFuncs::multiplot(m.infer.lst[[1]], m.infer.lst[[2]], m.infer.lst[[3]], cols = 3)

dev.off()



