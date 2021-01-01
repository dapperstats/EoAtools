    #
    # Function used for the processing
    #  
    #  posterior = which of the three posteriors 
    #    (FSFY = full site full year, FSMP = full site monitored period, 
    #     SAMP = searched area monitored period)
    #  PCIp = posterior credibility interval of interest
    #  digits = decimal places to leave 
    #  (1: in raw posterior, 2: in summary stats)
    #
#' @title Process EoA Lambda Posteriors
#' @description Collect and process Lambda posteriors from EoA.
#'   \cr
#'   Run this function after conducting a single year analysis.
#' @return \code{list} of posteriors of lambda.
#'
#' @export
#'
EoApostL <- function(){

  # grab the results list

    element_name <- "syresult"

    results <- getData(element_name)
    prob_obs <- results$prob_obs
    X <- results$X

    # Full Site Full Year
    # Full Site Monitored Period
    # Searched Area Monitored Period

      BabAnn <- results$BabAnn
      Bab <- results$Bab
      BabRaw <- results$BabRaw
    
  # return the list of functions created

  list(postL_Ann = posteriorL.ab(X, BabAnn[[1]], BabAnn[[2]]),
       postLpdf_Ann = posteriorLpdf.ab(X, BabAnn[[1]], BabAnn[[2]]),
       postL_Raw = posteriorL.ab(X, BabRaw[[1]], BabRaw[[2]]),
       postLpdf_Raw = posteriorLpdf.ab(X, BabRaw[[1]], BabRaw[[2]]),
       postL = posteriorL.ab(X, Bab[[1]], Bab[[2]]),
       postLpdf = posteriorLpdf.ab(X, Bab[[1]], Bab[[2]]))


}