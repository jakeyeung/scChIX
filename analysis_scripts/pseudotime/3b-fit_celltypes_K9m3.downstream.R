# Jake Yeung
# Date of Creation: 2021-08-12
# File: ~/projects/scChIX/analysis_scripts/pseudotime/3b-fit_celltypes_K9m3.downstream.R
#

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(JFuncs)
library(scchicFuncs)



# Load fits ---------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"
jmark <- "K9m3"
inf <- file.path(hubprefix, paste0("jyeung/data/dblchic/gastrulation/from_analysis/glm_fits_outputs/by_clusters/manual2nocenter_K36_K9m3_K36-K9m3/glm_poisson_fits_output.clusters.manual2nocenter_K36_K9m3_K36-K9m3.", jmark, ".RData"))
load(inf, v=T)


# Find top hits  ----------------------------------------------------------

pvals.vec <- lapply(jfits.lst, function(jfit){
  jfit$pval
}) %>%
  unlist()

plot(density(pvals.lst))


# Plot most significant ones ----------------------------------------------

qvals.vec <- p.adjust(pvals.vec)

qval.cutoff <- 10^-8

qvals.filt <- qvals.vec[which(qvals.vec <= qval.cutoff)]

# Cluster the outputs  ----------------------------------------------------

jbins.filt <- names(qvals.filt)
names(jbins.filt) <- jbins.filt

# plot a hit

set.seed(0)
jbin <- sample(jbins.filt, size = 1)

# get the parameters: maybe add SE too

bin2pval <- hash::hash(names(pvals.vec), pvals.vec)
bin2qval <- hash::hash(names(qvals.vec), qvals.vec)

merged.dat.lst <- lapply(jbins.filt, function(jbin){
  jpval <- bin2pval[[jbin]]
  jqval <- bin2qval[[jbin]]
  jparams.vec <- jfits.lst[[jbin]]
  estimate.i <- grepl(pattern = "^cluster.*Estimate$", names(jparams.vec))
  se.i <- grepl(pattern = "^cluster.*StdError$", names(jparams.vec))
  int.i <- grepl(pattern = "X.Intercept..Estimate", names(jparams.vec))

  estimate.dat <- data.frame(param = names(jparams.vec)[estimate.i], estimate = unlist(jparams.vec[estimate.i]) / log(2), stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(param = strsplit(param, split = "\\.")[[1]][[1]],
           param = gsub("^cluster", "", param))

  se.dat <- data.frame(param = names(jparams.vec)[se.i], se = unlist(jparams.vec[se.i]) / log(2), stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(param = strsplit(param, split = "\\.")[[1]][[1]],
           param = gsub("^cluster", "", param))

  merged.dat <- left_join(estimate.dat, se.dat, by = "param") %>%
    rowwise() %>%
    mutate(estimate = ifelse(abs(estimate) > 5, NA, estimate),
           se = ifelse(estimate == 0, NA, se))

  row.add <- data.frame("param" = "Epithelial", "estimate" = 0, "se" = 0, stringsAsFactors = FALSE)
  merged.dat <- rbind(merged.dat, row.add)

  merged.dat$bin <- jbin
  merged.dat$pval <- jpval
  merged.dat$qval <- jqval

  return(merged.dat)
})



merged.dat.filt.long <- lapply(jbins.filt, function(jbin){
  jdat <- merged.dat.lst[[jbin]]
  if (any(is.na(jdat$estimate))){
    return(NULL)
  } else{
    jdat$estimate.center <- jdat$estimate - mean(jdat$estimate)
    return(jdat)
  }
}) %>%
  bind_rows()


jbin <- names(which.min(pvals.vec))
jbin <- sample(jbins.filt, size = 1)

jbins.nonull <- unique(merged.dat.filt.long$bin)
jbin <- sample(jbins.nonull, size = 1)


ggplot(subset(merged.dat.filt.long, bin == jbin), aes(x = param, y = estimate, ymin = estimate - se, ymax = estimate + se)) +
  geom_point() +
  geom_errorbar(width = 0.25) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme_bw() +
  ggtitle(jbin) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# Do kmeans ---------------------------------------------------------------


ggplot(subset(merged.dat.filt.long, bin == jbin), aes(x = param, y = estimate.center, ymin = estimate.center - se, ymax = estimate.center + se)) +
  geom_point() +
  geom_errorbar(width = 0.25) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme_bw() +
  ggtitle(jbin) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


mat.estimates <- reshape2::dcast(data = merged.dat.filt.long, formula = bin ~ param, value.var = "estimate.center")
rownames(mat.estimates) <- mat.estimates$bin
mat.estimates$bin <- NULL

heatmap(as.matrix(mat.estimates))

outpdf <- "/home/jyeung/hub_oudenaarden/jyeung/data/tmp/heatmap_params_centered_K36.pdf"
pdf(outpdf, useDingbats = FALSE)
  heatmap3::heatmap3(as.matrix(mat.estimates))
dev.off()



