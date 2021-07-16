# Jake Yeung
# Date of Creation: 2019-06-18
# File: ~/projects/scchic_gastru/scripts_analysis/Rfunctions/QCFunctions.R
# QC functions

# Gastru samps: E7p5-EHF-CastBl6-H3K36me3-20190214-002_249

StageToNumeric <- function(stage, to.numeric = TRUE){
  # can be vector
  stage <- gsub("_EHF", "", stage)
  stage <- gsub("_EBLB", "", stage)
  stage <- gsub("p", "\\.", stage)
  stage <- sapply(stage, function(x) strsplit(x, "_")[[1]][[1]])
  if (to.numeric){
    stage.numer <- as.numeric(gsub("^E", "", stage))
    if (!all(!is.na(stage.numer))){
      print("Unxepcted NAs found")
      print(unique(stage))
    } 
    assertthat::assert_that(all(!is.na(stage.numer)))
    return(stage.numer)
  } else{
    return(stage)
  }
}


GetStage <- function(x, do.preprocess = FALSE){
  # preprocessed sample
  if (do.preprocess){
    x <- PreprocessSamp(x)
  } 
  strsplit(x, "-")[[1]][[1]]
}


GetPlateCoord <- function(cell, platecols = 24, is.zero.base = TRUE){
  # cell0 -> 1,1
  indx <- as.numeric(strsplit(cell, "cell")[[1]][[2]]) 
  if (is.zero.base){
    indx <- indx + 1
  }
  jcol <- indx %% platecols
  jcol <- ifelse(jcol == 0, platecols, jcol)
  jrow <- ceiling(indx / platecols)
  return(c(jrow, jcol))
}



PreprocessSamp <- function(x){
  # E8-13S-CastBl6-H3K27me3-20190315-001.filtered.RZcounts.csv -> combine E8-13S into E8_13S
  # E9p5-Bl6Cast-H3K36me3-H3K9me3-G1-20190404-002.filtered.RZcounts.csv -> combine into H3K36me3_H3K9me3
  x <- FixMarkName(x)
  if (GetLengthTime(x) == 2){
    timestr.orig <- paste(strsplit(x, "-")[[1]][1:2], collapse = "-")
    timestr.new <- gsub("-", "_", timestr.orig)
    x <- gsub(timestr.orig, timestr.new, x)
  } 
  if (GetNMarks(x) == 2){
    mark.indx <- which(startsWith(strsplit(x, "-")[[1]], prefix = "H3K"))
    markstr.orig <- paste(strsplit(x, "-")[[1]][mark.indx], collapse = "-")
    markstr.new <- gsub("-", "_", markstr.orig)
    x <- gsub(markstr.orig, markstr.new, x)
  } 
  return(x)
}

HasG1 <- function(x){
  # check for G1 string or not
  return(grepl("-G1-", x))
}

FixMarkName <- function(x){
  # E9p5-CastBl6-K36me3-K27me3-G1-20190404-001.filtered.RZcounts.csv -> E9p5-CastBl6-H3K36me3-H3K27me3-G1-20190404-001.filtered.RZcounts.csv
  return(gsub("-K", "-H3K", x))
}

FixMarkName2 <- function(x){
  # mm_K36me3-K27me3.rds -> mm_H3K36me3_H3K27me3.rds
  x <- gsub("-K", "-H3K", x)
  x <- gsub("_K", "_H3K", x)
  x <- gsub("-", "_", x)
  return(x)
}

GetNMarks <- function(x){
  # E9p5-Bl6Cast-H3K36me3-H3K9me3-G1-20190404-002.filtered.RZcounts.csv -> 2 marks
  # E8-13S-CastBl6-H3K27me3-20190315-001.filtered.RZcounts.csv -> 1 mark
  # assumes histone marks start with H3K
  return(length(which(startsWith(strsplit(x, "-")[[1]], prefix = "H3K"))))
}


GetLengthTime <- function(x){
  # E10-CastBl6-H3K36me3-G1-20192201-001.filtered.RZcounts.csv: length is 1
  # E8-9S-CastBl6-H3K9me3-20190315-002.filtered.RZcounts.csv: length is 2
  checkstr <- strsplit(x, "-")[[1]][[2]]
  strainstr <- c("CastBl6", "Bl6Cast")
  if (checkstr %in% strainstr){
    tlength <- 1
  } else {
    tlength <- 2
  }
  return(tlength)
}

GetMarkGastru <- function(x){
  # indx is 3 or 4 depending on sampname
  # E10-CastBl6-H3K36me3-G1-20192201-001.filtered.RZcounts.csv: length is 1
  # E8-9S-CastBl6-H3K9me3-20190315-002.filtered.RZcounts.csv: length is 2
  # if (GetLengthTime(x) == 1){
  #   indx <- 3
  # } else {
  #   indx <- 4
  # }
  indx <- 3
  return(strsplit(x, "-")[[1]][[indx]])
}

GetRepGastru <- function(x){
  # E10-CastBl6-H3K36me3-G1-20192201-001.filtered.RZcounts.csv: indx 6
  # E8-9S-CastBl6-H3K9me3-20190315-002.filtered.RZcounts.csv: length is also 6??
  indx <- length(strsplit(x, "-")[[1]])
  jrep <- strsplit(strsplit(x, "-")[[1]][[indx]], "_")[[1]][[1]]
  return(paste0("rep", jrep))
}

GetCellGastru <- function(x, shift = 0){
  indx <- length(strsplit(x, "_")[[1]])
  cell <- as.numeric(strsplit(x, "_")[[1]][[indx]]) + shift
  assertthat::assert_that(!is.na(cell))
  return(paste0("cell", as.character(cell)))
}

GetTime <- function(x){
  jtime <- strsplit(x, "-")[[1]][[1]]
  return(jtime)
}

GetStrain <- function(x){
  strain <- strsplit(x, "-")[[1]][[2]]
  return(strain)
}

ReadDinucGastru <- function(inf){
  dat <- fread(inf)[-1, ] %>%
    dplyr::rename(dinuc = sampleName) 
  dat.long <- gather(dat, key = "samp", value = "count", -dinuc) %>%
    rowwise() %>%
    mutate(samp = PreprocessSamp(samp)) %>%
    mutate(count = ifelse(is.na(count), 0, count),
           mark = GetMarkGastru(samp),
           repl = GetRepGastru(samp),
           cell = GetCellGastru(samp, shift = 0),
           time = GetTime(samp),
           strain = GetStrain(samp))
  return(dat.long)
}



# Get empty wells ---------------------------------------------------------


# #  1 index
GetEmptyWells <- function(indx = 0){
  if (indx == 1){
    indx.all <- seq(384)
    hascell.indx <- c(seq(1:356),seq(360:379)+360)
    empty.indx <- setdiff(indx.all, hascell.indx)
    empty.names <- paste("cell", empty.indx, sep = "")
  } else if (indx == 0){
    # 0 index
    indx.all <- seq(384) - 1
    hascell.indx <- c(seq(1:356),seq(360:379)+360) - 1
    empty.indx <- setdiff(indx.all, hascell.indx)
    empty.names <- paste("cell", empty.indx, sep = "")
  } else {
    warning("Index must be 0 or 1")
  }
  return(empty.names)
}


# Matrix functions --------------------------------------------------------

LoadMats <- function(infs, cells.keep, preprocess.names = TRUE){
  sparse.mats <- lapply(infs, function(inf) readRDS(inf))
  rnames <- lapply(sparse.mats, rownames)
  rnames.common <- Reduce(intersect, rnames)
  sparse.mats <- lapply(sparse.mats, function(x){
    return(x[rnames.common, ])
  })
  mat.merge <- do.call(cbind, sparse.mats)
  if (preprocess.names){
    colnames(mat.merge) <- sapply(colnames(mat.merge), PreprocessSamp)
  }
  cells.keep.i <- which(colnames(mat.merge) %in% cells.keep)
  mat.merge <- mat.merge[, cells.keep.i]
  assertthat::assert_that(ncol(mat.merge) > 0)
  return(mat.merge)
}


# Blacklist ---------------------------------------------------------------


LoadBlacklist <- function(inf = "data/blacklists/mm10.blacklist.bed.gz", asGR = TRUE){
  dat <- fread(inf, col.names = c("seqnames", "start", "end"))
  if (asGR){
    dat <- makeGRangesFromDataFrame(dat)
  }
  return(dat)
}

