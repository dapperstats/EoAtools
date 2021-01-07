


#' @title Measure overlap between pairs of CDFs
#'
#' @description Measures overlap between pairs of CDFs. First uses
#'  \code{numdiff} (\code{\link[pracma]{numderiv}}) to differentiate the CDFs,
#'   then calculates the classic Overlap Coefficient. 
#'
#' @param Fs \code{list} of CDFs 
#'
#' @param vals \code{numeric} vector for where \code{Fs} should be 
#'  differentiated
#'
#' @param digits rounding digits
#'
#' @param ... additional arguments to pass to \code{numdiff} 
#'  (\code{\link[pracma]{numderiv}})
#' 
#' @return \code{matrix} of pairwise overlap values
#'
#' @export
#'
overlap <- function(Fs, vals, digits = 10, ...){

  nFs <- length(Fs)
  dFs <- list()
  for(i in 1:nFs){
    dFs[[i]] <- numdiff(Fs[[i]], vals, ...)
  }

  O <- matrix(NA, nFs, nFs)

  for(i in 1:nFs){
    for(j in 1:nFs){
      O[i,j] <- sum(pmin(standardize(dFs[[i]]), standardize(dFs[[j]])))
    }
  }
  if(!is.null(names(Fs))){
    colnames(O) <- rownames(O) <- names(Fs)
  }
  round(O, digits)
}




standardize <- function(x){
  x/sum(x)
}

