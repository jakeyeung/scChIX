# Jake Yeung
# Date of Creation: 2021-11-17
# File: ~/projects/scChIX/analysis_scripts/annotate_bins_by_features/explore_features_promoter_enhancer_load_obj.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

inrds <- "/Users/yeung/data/dblchic/annotate_bins_by_promoter_enh_repeats/dat_sum_by_regions.rds"
dat.lst <- readRDS(inrds)
dat.long <- dat.lst %>%
  bind_rows()

ggplot(dat.long %>% filter(variable %in% c("proms", "enhs")), aes(x = variable, y = value / total.count)) + 
  geom_boxplot() + 
  facet_wrap(~label) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Load prob mat and associate genes to prob mat?  ------------------------





