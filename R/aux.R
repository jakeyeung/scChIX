# Jake Yeung
# Date of Creation: 2021-04-16
# File: ~/projects/scChIX/R/aux.R
#


AnnotateBins2.R4 <- function (terms.mat, top.thres = 0.995, inf.tss = "/Users/yeung/data/scchic/tables/gene_tss_winsize.50000.bed", 
                              txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", 
                              chromos.keep = c(paste("chr", seq(19), sep = ""), "chrX", 
                                               "chrY"), skip.split = FALSE) 
{
  assertthat::assert_that(file.exists(inf.tss))
  # assertthat::assert_that(class(terms.mat) == "matrix")
  regions <- data.frame(seqnames = sapply(colnames(terms.mat), 
                                          GetChromo), start = sapply(colnames(terms.mat), GetStart), 
                        end = sapply(colnames(terms.mat), GetEnd), stringsAsFactors = FALSE)
  rownames(regions) <- colnames(terms.mat)
  print("Chromos to keep")
  print(chromos.keep)
  regions <- subset(regions, seqnames %in% chromos.keep)
  regions.range <- makeGRangesFromDataFrame(as.data.frame(regions))
  regions.annotated <- as.data.frame(annotatePeak(regions.range, 
                                                  TxDb = txdb, annoDb = annodb))
  regions.annotated$region_coord <- names(regions.range)
  topic.regions <- lapply(seq(nrow(terms.mat)), function(clst) {
    return(SelectTopRegions(terms.mat[clst, ], colnames(terms.mat), 
                            method = "thres", method.val = top.thres))
  })
  print(paste("Using TSS definitions from:", inf.tss))
  terms.long <- data.frame(term = colnames(terms.mat), as.data.frame(t(terms.mat)), 
                           stringsAsFactors = FALSE) %>% gather(key = "topic", value = "weight", 
                                                                -term) %>% mutate(topic = gsub("X", "", topic)) %>% group_by(topic) %>% 
    arrange(desc(weight)) %>% mutate(rnk = seq(length(weight))) %>% 
    rowwise()
  terms.filt.top <- terms.long %>% rowwise()
  tss.dat <- fread(inf.tss, col.names = c("seqnames", "start", 
                                          "end", "tssname"))
  tss.dat$gene <- sapply(tss.dat$tssname, function(x) strsplit(x, 
                                                               ";")[[1]][[2]])
  annots.biomart <- regions.annotated %>% mutate(midpt = start + 
                                                   (end - start)/2) %>% filter(region_coord %in% terms.filt.top$term)
  annots.gr <- makeGRangesFromDataFrame(annots.biomart %>% 
                                          dplyr::select(seqnames, start, end, SYMBOL, region_coord, 
                                                        ENSEMBL), keep.extra.columns = TRUE)
  annots.tss.gr <- makeGRangesFromDataFrame(tss.dat, keep.extra.columns = TRUE)
  out <- findOverlaps(annots.tss.gr, annots.gr, type = "within")
  out2 <- findOverlaps(annots.gr, annots.tss.gr, type = "any")
  out2.df = data.frame(annots.gr[queryHits(out2), ], annots.tss.gr[subjectHits(out2), 
  ]) %>% mutate(midpt = start + round(width/2), midpt.1 = start.1 + 
                  round(width.1/2), dist.to.tss = midpt.1 - midpt)
  out2.df.closest <- out2.df %>% group_by(region_coord) %>% 
    filter(abs(dist.to.tss) == min(abs(dist.to.tss)))
  terms.new <- paste(out2.df.closest$region_coord, out2.df.closest$gene, 
                     sep = ";")
  terms.hash <- hash::hash(out2.df.closest$region_coord, terms.new)
  terms.annot <- terms.filt.top %>% mutate(termgene = ifelse(!is.null(terms.hash[[term]]), 
                                                             terms.hash[[term]], NA))
  terms.filt <- terms.filt.top %>% mutate(termgene = ifelse(!is.null(terms.hash[[term]]), 
                                                            terms.hash[[term]], NA)) %>% filter(!is.na(termgene))
  if (!skip.split) {
    terms.filt <- terms.filt %>% mutate(gene = sapply(termgene, 
                                                      function(x) strsplit(x, ";")[[1]][[2]])) %>% group_by(gene)
  }
  return(list(topic.regions = topic.regions, regions.annotated = regions.annotated, 
              terms.annot = terms.annot, out2.df.closest = out2.df.closest, 
              terms.filt = terms.filt))
}
ClipLast <- function(x, jsep = "-", jsep.out = NULL){
  # B6-13W1-BM-H3K4me3-1_269 -> B6-13W1-BM-H3K4me3
  if (is.null(jsep.out)){
    jsep.out <- jsep
  }
  jsplit <- strsplit(x, jsep)[[1]]
  # remove last one
  N <- length(jsplit) - 1
  return(paste(jsplit[1:N], collapse = jsep.out))
}



SoftMax <- function(x, return.log = TRUE, logfn = log){
  # numericallys table softmax by subtracting the maximum first
  # https://stackoverflow.com/questions/42599498/numercially-stable-softmax
  # x value are in log if .log = TRUE
  # numer <- log(exp(x - max(x)))
  # denom <- log(sum(exp(x - max(x))))
  numer <- logfn(exp(x - max(x)))
  denom <- logfn(sum(exp(x - max(x))))
  plog <- numer - denom
  if (return.log){
    return(plog)
  } else {
    return(exp(plog))
  }
}

