# Jake Yeung
# cluster_make_lda_countmat_metadata_objects.R
# 2021-07-16
# DESCRIPTION
# 
#     Cluster and make objects for setting up scChIX
# 
# FOR HELP
# 
#     Rscript cluster_make_lda_countmat_metadata_objects.R --help
# 
# AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
# LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
# CREATED ON:  2021-07-16
# LAST CHANGE: see git log
# LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

suppressPackageStartupMessages(library("argparse"))

library(assertthat)

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
# library(scChIX)

library(topicmodels)
library(hash)
library(igraph)
library(umap)



# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument('-marks', metavar='Antibodies', nargs = "+", required=TRUE,
                                            help='Input LDA')
parser$add_argument('-infiles', metavar='INFILES', nargs = "+", required=TRUE,
                                            help='Input LDA files, one for each mark respectively')
parser$add_argument('-outdir', metavar='OUTDIR', required=TRUE,
                                            help='Output prefix to write countmat, LDA, metadata')
parser$add_argument('-n_neighbors', metavar='integer', type='integer', default = 30,
                                            help='Nearest neighbors for umap. Default 30')
parser$add_argument('-min_dist', metavar='Float', type='double', default = 0.1, 
                                            help='Floating point')
parser$add_argument('-random_state', metavar='integer', type='integer', default = 123, 
                                            help='Floating point')
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                        help="Print extra output [default]")

                                        
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

jsettings <- umap.defaults
jsettings$n_neighbors <- args$n_neighbors
jsettings$min_dist <- args$min_dist
jsettings$random_state <- args$random_state

# print some progress messages to stderr if "quietly" wasn't requested
if ( args$verbose ) { 
    print("Arguments:")
    print(args)
}

jmarks <- args$marks
names(jmarks) <- jmarks
infrobjs <- args$infiles
assertthat::assert_that(all(file.exists(infrobjs)))
assertthat::assert_that(length(infrobjs) == length(jmarks))

names(infrobjs) <- jmarks

outdir <- args$outdir

for (jmark in jmarks){

  # outmat <- file.path(outdir, paste0("countmat_output_filt.", jmark, ".", Sys.Date(), ".rds"))
  # outmeta <- file.path(outdir, paste0("celltyping_output_filt.", jmark, ".", Sys.Date(), ".rds"))
  # outlda <- file.path(outdir, paste0("lda_output_filt.", jmark, ".", Sys.Date(), ".rds"))
  # outpdf <- file.path(outdir, paste0("plots.", jmark, ".", Sys.Date(), ".pdf"))
  outmat <- file.path(outdir, paste0("countmat_output_filt.", jmark, ".rds"))
  # outmeta <- file.path(outdir, paste0("celltyping_output_filt.", jmark, "..rds"))
  # outlda <- file.path(outdir, paste0("lda_output_filt.", jmark, "..rds"))
  # outpdf <- file.path(outdir, paste0("plots.", jmark, "..pdf"))
  outmeta <- file.path(outdir, paste0("celltyping_output_filt.", jmark, ".rds"))
  outlda <- file.path(outdir, paste0("lda_output_filt.", jmark, ".rds"))
  outpdf <- file.path(outdir, paste0("plots.", jmark, ".pdf"))

  # skip if all exists
  outall <- c(outmat, outmeta, outlda, outpdf)
  if (all(file.exists(outall))){
    print(outall)
    print("All files exist, skipping")
    next
  }
  infrobj <- infrobjs[[jmark]]

  load(infrobj, v=T)

  tm.result <- posterior(out.lda)

  print("Running UMAP and Louvain...")

  dat.umap.long <- DoUmapAndLouvain(tm.result$topics, jsettings = jsettings)

  print("Annotating plate, experiment, stage, cluster ... ")
  # annotate
  dat.umap.long <- dat.umap.long %>%
    rowwise() %>%
    mutate(plate = ClipLast(cell, jsep = "_"),
           experi = ClipLast(plate, jsep = "-"),
           stage =  strsplit(cell, "-")[[1]][[1]],
           cluster = paste("cluster", louvain, sep = "")) %>%
    dplyr::select(-louvain)


  # Write tables  -----------------------------------------------------------

  print("Making plots...")
  m <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = cluster)) +
    geom_point() +
    theme_bw() +
    ggtitle(paste(jmark)) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  m1 <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = cluster)) +
    geom_point() +
    theme_bw() +
    facet_wrap(~plate) +
    ggtitle(paste(jmark)) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  m2 <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = cluster)) +
    geom_point() +
    theme_bw() +
    facet_wrap(~experi) +
    ggtitle(paste(jmark)) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  m3 <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = cluster)) +
    geom_point() +
    theme_bw() +
    facet_wrap(~stage) +
    ggtitle(paste( jmark)) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


  print("Saving rds objects...")
  saveRDS(count.mat, file = outmat)
  saveRDS(dat.umap.long, file = outmeta)
  saveRDS(out.lda, file = outlda)

  pdf(outpdf, useDingbats = FALSE)
    print(m)
    print(m1)
    print(m2)
    print(m3)
  dev.off() 
}

