% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Sn_A.R
\name{Sn_A}
\alias{Sn_A}
\title{Computes the Moebius Cramer-von Mises statistics for the test of independence between random variables}
\usage{
Sn_A(x, trunc.level)
}
\arguments{
\item{trunc.level}{Only subsets of cardinality <= trunc.level (default=2) are considered for the Moebius statistics.}

\item{data}{Matrix of observations.}
}
\value{
\item{stats}{Cramer-von Mises Moebius statistics}

\item{cardA}{cardinality of subsets}

\item{M}{Matrix for multitpliers bootstrap for stats}

\item{Asets}{vector of (0,1) for Moebius subsets}

\item{Sn}{Cramer-von Mises Sn statistic}

\item{J}{Matrix for multipliers bootstrap for Sn}
}
\description{
This function he Moebius Cramer-von Mises statistics for a tests of randomnes for observations X(1), ...l, X(p).
}
\references{
Genest, Neslehova, Remillard & Murphy (2019). Testing for independence in arbitrary distributions

#'@examples
X <- rnorm(100,5)
out <- Sn_A(X,3)
}
\keyword{internal}