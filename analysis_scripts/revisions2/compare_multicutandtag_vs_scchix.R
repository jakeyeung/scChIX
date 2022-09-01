# Jake Yeung
# Date of Creation: 2022-08-30
# File: ~/projects/scChIX/analysis_scripts/revisions2/compare_multicutandtag_vs_scchix.R
# Compare fragments per cell between the two methods

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

# Load MCT frags per cell -------------------------------------------------

indir.mct <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChIX_revisions/Gopalan_et_al/raw_data/bams/outputs_mct"
indir.EtOH <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChIX_revisions/Gopalan_et_al/raw_data/bams/outputs_EtOH"
marks.mct <- c("H3K27ac", "H3K27me3"); names(marks.mct) <- marks.mct
marks.EtOH <- c("K27m3", "K9m3"); names(marks.EtOH) <- marks.EtOH

infs.mct <- lapply(marks.mct, function(mark){
  inf.tmp <- file.path(indir.mct, paste0("fragments_per_cell_", mark, "_filt.csv"))
  return(inf.tmp)
})

infs.EtOH <- lapply(marks.EtOH, function(mark){
  inf.tmp <- file.path(indir.EtOH, paste0("fragments_per_cell_EtOH_", mark, ".csv"))
  return(inf.tmp)
})

dats.mct <- lapply(marks.mct, function(jmark){
  jinf <- infs.mct[[jmark]]
  fread(jinf, header = FALSE) %>%
    mutate(dataset = "Multi-CUT&TAG", 
           mark = jmark) %>%
    dplyr::rename("barcode" = "V1",
                  "fragments_per_cell" = "V2") 
}) %>%
  bind_rows()

dats.EtOH <- lapply(marks.EtOH, function(jmark){
  jinf <- infs.EtOH[[jmark]]
  fread(jinf, header = FALSE) %>%
    mutate(dataset = "scChIX",
           mark = jmark) %>%
    dplyr::rename("barcode" = "V1",
                  "fragments_per_cell" = "V2")
}) %>%
  bind_rows()

# Make tables -------------------------------------------------------------

dats.mct.filt <- dats.mct %>% filter(fragments_per_cell > 100)
barcodes.keep <- intersect(subset(dats.mct.filt, dataset == "Multi-CUT&TAG" & mark == "H3K27ac")$barcode, 
                           subset(dats.mct.filt, dataset == "Multi-CUT&TAG" & mark == "H3K27me3")$barcode)
dats.mct.filt2 <- subset(dats.mct.filt, barcode %in% barcodes.keep)

dats.merged <- rbind(dats.mct.filt2, dats.EtOH) %>%
  mutate(mark = gsub("K27m3", "H3K27me3", mark), 
         mark = gsub("K9m3", "H3K9me3", mark))

options(scipen=100000)

dats.merged.summary <- dats.merged %>%
  group_by(dataset, mark) %>%
  summarise(ncells = length(barcode),
            mean_fragments_per_cell = mean(fragments_per_cell))

print(dats.merged.summary)

fwrite(dats.merged.summary, file = paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChIX_revisions/benchmarking/compare_with_scchix_reprocessed_summary.", Sys.Date(), ".txt"))
pdf(paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChIX_revisions/benchmarking/compare_with_scchix_reprocessed.", Sys.Date(), ".pdf"), useDingbats = FALSE)

ggplot(dats.merged, aes(x = dataset, fill = mark, y = fragments_per_cell)) + 
  scale_y_log10() + 
  geom_boxplot() + 
  xlab("") + 
  ylab("Fragments per cell") + 
  theme_bw(24) + 
  theme(aspect.ratio=2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

dev.off()
