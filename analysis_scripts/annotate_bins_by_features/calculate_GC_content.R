# Jake Yeung
# Date of Creation: 2021-11-18
# File: ~/projects/scChIX/analysis_scripts/annotate_bins_by_features/calculate_GC_content.R
# 

rm(list=ls())

library(scchicFuncs)
library(genoset)
library(BSgenome.Mmusculus.UCSC.mm10)

inf.rds <- "/Users/yeung/Dropbox/macbookpro_data/data/dblchic/unfixed_revisions/metadata_5kb/rownames_5kb_coords.2021-11-18.rds"
coords <- readRDS(inf.rds)

# remove chr?
# coords <- gsub("^chr", "", coords)

# coords <- params.long.merge$bin
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")

print(length(coords))
gr.dat <- data.frame(seqnames = sapply(coords, GetChromo), start = sapply(coords, GetStart), end = sapply(coords, GetEnd), bname = coords) %>%
  filter(seqnames %in% jchromos)

coords.filt <- as.character(gr.dat$bname)
print(head(gr.dat))
gr <- GenomicRanges::makeGRangesFromDataFrame(gr.dat)
names(gr) <- coords.filt

# print chromos
chromos.uniq <- unique(sapply(coords.filt, function(x) strsplit(x, ":")[[1]][[1]]))
print(chromos.uniq)

gr.gc <- calcGC(object = gr, bsgenome = BSgenome.Mmusculus.UCSC.mm10)
gr.gc.dat <- data.frame(bname = names(gr.gc), gc = gr.gc, stringsAsFactors = FALSE)

# split by chromosome? 

save(gr.gc.dat, file = paste0("/Users/yeung/data/dblchic/unfixed_revisions/metadata_5kb/gcs_genomewide_5kb.", Sys.Date(), ".RData"))
