% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Codes.R
\name{APE}
\alias{APE}
\title{Multi-fold cross validation for the MXFAR model}
\usage{
APE(SeriesNum, Tlength, y, u, p, d, r, Q, bwp = 0.1, numpoints = 50)
}
\arguments{
\item{SeriesNum}{A numeric vector stating the number of series in a 'group' data}

\item{Tlength}{A numeric scalar stating the number of time points per series}

\item{y}{A matrix containing the 'stacked' multivariate time series data from multiple series}

\item{u}{A vector containing the 'stacked' reference signal series}

\item{p}{The order of MXFAR(p) model}

\item{d}{The lag of the reference signal}

\item{r}{The number of time points for predictions}

\item{Q}{The number of segments considered for the time series}

\item{bwp}{The bandwidth proportion from 0 to 1 (with respect to the range of observations)}

\item{numpoints}{The number of discretized points in the range of the reference signal 'u' that will be used to evaluate the local linear approximation}
}
\value{
The calculated APE metric.
}
\description{
This function calculates the accumulated prediction error (APE) from a fitted MXFAR model
}
