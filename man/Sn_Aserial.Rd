% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Sn_Aserial.R
\name{Sn_Aserial}
\alias{Sn_Aserial}
\title{Computes the Moebius Cramer-von Mises statistics for the test of randomness}
\usage{
Sn_Aserial(x, p, trunc.level)
}
\arguments{
\item{x}{Time series.}

\item{p}{Number of consecutive observations for the test.}

\item{trunc.level}{Only subsets of cardinality <= trunc.level (default=2) are considered for the Moebius statistics.}
}
\value{
\item{stats}{Cramer-von Mises Moebius statistics}

\item{cardA}{Cardinality of subsets}

\item{M}{Matrix for multitpliers bootstrap for stats}

\item{Asets}{Vector of (0,1) for Moebius subsets}

\item{Sn}{Cramer-von Mises Sn statistic}

\item{J}{Matrix for multipliers bootstrap for Sn}
}
\description{
This function computes Moebius Cramer-von Mises statistics for a tests of randomness for observations X(1), ..., X(p).
}
\references{
B.R Nasri (2022). Tests of serial dependence for arbitrary distributions

#'@examples
X <- SimAR1Poisson(c(5,0.2),100)
out <- Sn_Aserial(X,5,3)
}
\keyword{internal}
