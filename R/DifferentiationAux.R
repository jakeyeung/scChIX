# Jake Yeung
# Date of Creation: 2022-03-26
# File: ~/projects/scChIX/R/DifferentiationAux.R
# 

#' Get day given cell data frame
#' 
#' @param dat.meta Input dataframe containing column name "cell" indicating cell names for macrophage differentiation experiment
#' @return metadata annotated with day
#' @examples
#' dat.meta.annot <- GetDayFromCellDataFrame(dat.meta)
#' @export
GetDayFromCellDataFrame <- function(dat.meta){
  # assumes "cell" is in column name
  dat.meta.annot <- dat.meta %>%
    rowwise() %>%
    mutate(platename = scChIX::ClipLast(cell, jsep = "_"),
           experi = scChIX::ClipLast(platename, jsep = "-"),
           cellname = paste0("cell", strsplit(cell, split = "_")[[1]][[2]]),
           jcol = scChIX::GetPlateCoord(cellname, is.zero.base = FALSE)[[2]],
           jrow = scChIX::GetPlateCoord(cellname, is.zero.base = FALSE)[[1]],
           day = ceiling(jcol / 3))
  return(dat.meta.annot)
}

#' Get day given cell name
#' 
#' @param cell cell name for macrophage differentiation experiment
#' @return day from which cell was extracted
#' @examples
#' day <- GetDadyFromCellByRow(dat.meta$cell[[1]])
#' @export
GetDayFromCellByRow <- function(cell){
  # assumes "cell" is in column name
  platename <- scChIX:::ClipLast(cell, jsep = "_")
  experi <- scChIX:::ClipLast(platename, jsep = "-")
  cellname <- paste0("cell", strsplit(cell, split = "_")[[1]][[2]])
  jcol <- scChIX:::GetPlateCoord(cellname, is.zero.base = FALSE)[[2]]
  jrow <- scChIX:::GetPlateCoord(cellname, is.zero.base = FALSE)[[1]]
  day <- ceiling(jcol / 3)
  return(day)
}
