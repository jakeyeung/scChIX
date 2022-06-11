# Jake Yeung
# Date of Creation: 2021-12-30
# File: ~/projects/scChIX/analysis_scripts/H3K9me3-H3K36me3_gastrulation_analysis/14-plot_gene_expression_Cao.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(forcats)

outpdf <- paste0("/Users/yeung/data/dblchic/gastrulation/primetime_plots/marker_genes_mRNA_Cao.", Sys.Date(), ".pdf")

# load meta for colors
inf.meta <- "/Users/yeung/Dropbox/macbookpro_data/data/dblchic/gastrulation/metadata/metadata_cell_cluster_with_clustercol.K36.2021-12-03.txt"
dat.meta <- fread(inf.meta)
  
clst2col <- hash::hash(dat.meta$cluster, dat.meta$clustercol)

inf.rnaseq <- "/Users/yeung/Dropbox/macbookpro_data/data/dblchic/gastrulation/from_cluster/cao_rnaseq_reference_pbulk_zscore_ctypes_merged.rds"
dat.ref <- readRDS(inf.rnaseq)

jgene <- "Pcdhga9"
jgene <- "Cdh6"

jgene <- "Sptb"

jgene <- "Diaph1"

jgene <- "Chil3"


jsub <- subset(dat.ref, grepl("Pcdhga", gene))

jgene <- "Tpm1"
jgene <- "Gata6"


jgene <- "Chil3"

jgene <- "Nell2"
clstrs.order.orig <- c("Erythroid", "WhiteBloodCells", "Endothelial", "NeuralTubeNeuralProgs", "Neurons", "SchwannCellPrecursor", "Epithelial", "Stromal", "ConnectiveTissueProg")
names(clstrs.order.orig) <- clstrs.order.orig

clstrs.order <- c("Erythroid", "WhiteBloodCells", "Endothelial", "NeuralTubeNeuralProgs", "Neurons", "SchwannCellPrecursor", "Epithelial", "MesenchymalProgs", "Cardiomyocytes")
pbulk.order <- c("Erythroid", 
                 "X31.White.blood.cells", 
                 "X20.Endothelial.cells", 
                 "NeuralTubeAndProgs", 
                 "NeuronalProgenitors", 
                 "X23.Schwann.cell.precursor", 
                 "X6.Epithelial.cells", 
                 "MesenchymalProgenitors", 
                 "X34.Cardiac.muscle.lineages")

dat.toappend <- data.frame(celltype.pretty = clstrs.order, clusters.orig = clstrs.order.orig, celltype = pbulk.order, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(clstcol = clst2col[[clusters.orig]])
  

# pbulk2clstr <- hash::hash(pbulk.order, clstrs.order)
# pbulk2clstr.orig <- hash::hash(pbulk.order, clstrs.order.orig)
jsub <- dat.ref %>% filter(gene == jgene & celltype %in% pbulk.order) %>%
  left_join(dat.toappend)
  # rowwise() %>%
  # mutate(celltype.pretty = AssignHash(x = celltype, jhash = pbulk2clstr, null.fill = celltype),
  #        celltype.pretty.orig = AssignHash(x = celltype, jhash = pbulk2clstr.orig, null.fill = celltype),)
  #        clstcol = clst2col[[celltype.pretty.orig]])

jaspect.ratio <- 2.5
pdf(outpdf, useDingbats = FALSE)


ggplot(jsub, aes(x = forcats::fct_reorder(.f = celltype.pretty, .x = counts, .fun = median, .desc = TRUE), y = counts, fill = clstcol)) + 
  geom_col() + 
  ggtitle(jgene) + 
  xlab("") + 
  scale_fill_identity() + 
  ylab("log2(CPM + 1)") + 
  theme_bw() + 
  geom_hline(yintercept = median(dat.ref$counts), linetype = "dotted") + 
  theme(aspect.ratio=jaspect.ratio, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


# Take all differentially expressed genes and plot  -----------------------

# "/Users/yeung/Dropbox/macbookpro_data/data/dblchic/gastrulation/H3K36me3_H3K9me3_celltyping/marker_genes_textoutput.K9m3.2021-12-02_Erythroid.txt"
ctypes <- c("Erythroid", "NonBlood", "WhiteBloodCells"); names(ctypes)<- ctypes

ctype.dir <- "/Users/yeung/data/dblchic/gastrulation/H3K36me3_H3K9me3_celltyping"

infs.lst <- lapply(ctypes, function(jctype){
  file.path(ctype.dir, paste0("marker_genes_textoutput.K9m3.2021-12-02_", jctype, ".txt"))
})

dat.markers.lst <- lapply(infs.lst, function(jinf){
  fread(jinf)
})

coords.lst <- lapply(dat.markers.lst, function(jdat){
  subset(jdat, avg_logFC > 0)$V1
})

# assign genes to nearest coords

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)

inf.tsspretty <- "/Users/yeung/Dropbox/macbookpro_data/data/databases/refseq/MmRefseqTss.chromorenamed.50000.again.cut.pretty.bed"
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")

bins.annot.lst <- lapply(coords.lst, function(coords.vec){
  print(head(coords.vec))
  bins.annot.tmp <- AnnotateCoordsFromList(coords.vec = coords.vec, inf.tss = inf.tsspretty, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos)
})

lapply(bins.annot.lst, function(jdat) dim(jdat$regions.annotated))

genes.lst <- lapply(bins.annot.lst, function(jdat) jdat$regions.annotated$SYMBOL)

genes.random <- sample(dat.ref$gene, size = 2000)

jctype <- "Erythroid"
# jsub <- subset(dat.ref, gene %in% genes.random)
jsub <- subset(dat.ref, gene %in% genes.lst[[jctype]] & celltype %in% pbulk.order) %>%
  left_join(., dat.toappend)

jsub.random <- subset(dat.ref, gene %in% genes.random & celltype %in% pbulk.order) %>%
  left_join(., dat.toappend)
  
m <- ggplot(jsub, aes(x = forcats::fct_reorder(.f = celltype.pretty, .x = zscore, .fun = median, .desc = TRUE), y = zscore, fill = clstcol)) + 
  geom_point() + 
  geom_boxplot() + 
  ggtitle(jctype) + 
  scale_fill_identity() + 
  xlab("") + 
  theme_bw() + 
  theme(aspect.ratio=jaspect.ratio, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
print(m)

m.random <- ggplot(jsub.random, aes(x = forcats::fct_reorder(.f = celltype.pretty, .x = zscore, .fun = median, .desc = TRUE), y = zscore, fill = clstcol)) + 
  geom_point() + 
  geom_boxplot() + 
  ggtitle("Random") + 
  scale_fill_identity() + 
  theme_bw() + 
  xlab("") + 
  theme(aspect.ratio=jaspect.ratio, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
print(m.random)

dev.off()


