# Jake Yeung
# Date of Creation: 2022-03-26
# File: ~/projects/scChIX/R/DifferentiationAux.R
# 

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

GetDayFromCellByRow <- function(cell){
  # assumes "cell" is in column name
  platename <- scChIX::ClipLast(cell, jsep = "_")
  experi <- scChIX::ClipLast(platename, jsep = "-")
  cellname <- paste0("cell", strsplit(cell, split = "_")[[1]][[2]])
  jcol <- scChIX::GetPlateCoord(cellname, is.zero.base = FALSE)[[2]]
  jrow <- scChIX::GetPlateCoord(cellname, is.zero.base = FALSE)[[1]]
  day <- ceiling(jcol / 3)
  return(day)
}
