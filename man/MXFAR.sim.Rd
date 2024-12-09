% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Codes.R
\name{MXFAR.sim}
\alias{MXFAR.sim}
\title{Multivariate MXFAR simulation}
\usage{
MXFAR.sim(SeriesNum, Tlength, d, Y_d, Fmat, reff.Sigma, eps.Sigma, nburn = 500)
}
\arguments{
\item{SeriesNum}{A numeric scalar stating the number of individual series to be simulated}

\item{Tlength}{A numeric scalar stating the number of time points per series}

\item{d}{The lag of the reference signal}

\item{Y_d}{The index of the time series to be used as endogenous reference signal. Put '0' for an exogenous standard normal reference signal}

\item{Fmat}{A list containing all functional coefficient matrices. The length of 'Fmat' dictates the order of the FAR model and the dimension of each matrix should match the dimension of the multivariate time series.}

\item{reff.Sigma}{A numeric vector of variances for the random effects in 'Fmat'.}

\item{eps.Sigma}{A numeric vector of random error variances for the individual series.}

\item{nburn}{The number of time points that will be discarded after simulation to ensure convergence of the Markov Chain.}
}
\value{
A list including the simulated time series and reference signal series.
}
\description{
This function simulates multiple multivariate time series from a multivariate MXFAR model with a chosen reference signal.
}
