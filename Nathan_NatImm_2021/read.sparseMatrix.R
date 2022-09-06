
library(dplyr)
library(stringr)


#' Read sparseMatrix efficiently
#'
#' Reading a big matrix requires too much memory. If matrix is sparse, read it line-by-line and convert to sparseMatrix. Assume first row and first column are labels
#'
#' @param file file path
#' @param batchSize number of lines to process per batch 
read.sparseMatrix = function(file, batchSize=100){

  myCon = gzfile(file, open="r")

  lst = list()
  i = 1
  repeat{

    if( i == 1){
      # read 1 line
      pl = readLines(myCon, n = 1)
      # if first line, set colnames
      colNames = str_split(pl, "\t")[[1]][-1]
    }else{
      pl = readLines(myCon, n = batchSize)

      # If the line is empty, exit.
      if(identical(pl, character(0))){break} 

      # else set values and rownames   
      # using fixed() is faster, but assumes only ASCII    
      res = str_split(pl, fixed("\t"))
      rm(pl)

      mat = lapply(res, function(x){
        values = as(as.integer(x[-1]), "sparseMatrix")
        colnames(values) = x[1]
        values
        })  
      rm(res)
      mat = t(do.call(cbind, mat))    
      lst[[i]] = mat
      cat("\r", i*batchSize, "    ")
    }
    i = i + 1
  }
  close(myCon)

  counts = do.call(rbind, lst)
  rm(lst)
  colnames(counts) = colNames
  counts
}
