#' @title Posterior of M
#' @description Collect and process the posterior M from EoA.
#'
#' @param posterior Which of the three posteriors
#' @param PCIp posterior credibility interval of interest
#' @param digits decimal places to leave
#' @param ... ...
#' @return \code{list} of EmpiricalPosterior and SummaryStats
#'
#' @export
#'
EoApostM <- function(posterior = "FSFY", PCIp = c(0.025, 0.975), 
                     digits = c(5, 2), ...){

  # extract the three posterior distributions:
  # Full Site Full Year
    
    pm_FSFY <- getData("pMgX.ann")

  # Full Site Monitored Period

    pm_FSMP <- getData("pMgX")

  # Searched Area Monitored Period

    pm_SAMP <- getData("pMgX.raw")

  # expand them so that they all cover the same mortality values  
  # the maximum number of mortality values that have a non-0 probability
  # across the three posteriors

    mlv <- max(c(length(pm_FSFY), length(pm_FSMP), length(pm_SAMP)))

  # expand the posteriors according to mlv

    pm_FSFY <- c(pm_FSFY, rep(0, mlv - length(pm_FSFY)))
    pm_FSMP <- c(pm_FSMP, rep(0, mlv - length(pm_FSMP)))
    pm_SAMP <- c(pm_SAMP, rep(0, mlv - length(pm_SAMP)))

  # combine into a matrix

    pmmat <- matrix(c(pm_FSFY, pm_FSMP, pm_SAMP), ncol = 3, byrow = F)  

  # create a vector of possible mortalities according to mlv

    nM <- 0 : (mlv - 1)

  # select out the specific posterior of interest 
  #  pm is the point probability, cpm is the cumulative probabilities

    posteriors <- c("FSFY", "FSMP", "SAMP")
    pp <- which(posteriors == posterior)
    pm <- pmmat[ , pp]
    cpm <- cumsum(pm)

    # summarize

  # expected value

    EV <- round(sum(nM * pm), digits[2])

  # highest posterior probability value(s)
  #   if there is a tie, all values are reported with comma separation

    HPPV <- nM[which(pm == max(pm))]
    HPPV[which(length(HPPV) > 1)] <- paste(HPPV, collapse = ", ")
    HPPV <- HPPV[1]

  # median

    MV <- nM[which(cpm >= 0.5)][1]
    
  # PCI

    PCI1 <- nM[which(cpm >= PCIp[1])][1]
    PCI2 <- nM[which(cpm >= PCIp[2])][1]
    PCI <- paste(PCI1, PCI2, sep = ", ")

  # combine and return

    outp <- c(EV, HPPV, MV, PCI)
    names(outp) <- c(ExpectedValue = EV, 
                   "HighestPosteriorProbabilityValue(s)" = HPPV,
                   Median = MV,
                   PosteriorCredibilityInterval = PCI)

    OUTP <- vector("list", 2)
    OUTP[[1]] <- data.frame(Mortalities = nM, 
                            Probability = round(pm, digits[1]))
    OUTP[[2]] <- outp
    names(OUTP) <- c("EmpiricalPosterior", "SummaryStats")

    OUTP
}