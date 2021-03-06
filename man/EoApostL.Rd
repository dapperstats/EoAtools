% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EoApostL.R
\name{EoApostL}
\alias{EoApostL}
\alias{EoApostL_parameters}
\title{Process EoA Lambda Posteriors}
\usage{
EoApostL()

EoApostL_parameters(X = NULL, Ba = NULL, Bb = NULL, ids)
}
\arguments{
\item{X}{X}

\item{Ba}{Ba}

\item{Bb}{Bb}

\item{ids}{ids}
}
\value{
\code{list} of posteriors of lambda.
}
\description{
Collect and process Lambda posteriors from EoA.
  \cr
  \code{EoApostL} is for running after conducting a single year analysis 
  and grabs the resulting output in the active R session.
  \cr
  \code{EoApostL_parameters} is for working with already fitted parameters
  from potentially multiple models.
}
