% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{calcSWFCCentroids}
\alias{calcSWFCCentroids}
\title{Calculate the centroids of SFCM}
\usage{
calcSWFCCentroids(data, wdata, belongmatrix, m, alpha)
}
\arguments{
\item{data}{A matrix representing the observed data with n rows and p columns}

\item{wdata}{A matrix representing the lagged observed data with nrows and p columns}

\item{belongmatrix}{A n X k matrix giving for each observation n, its
probability to belong to the cluster k}

\item{m}{An integer representing the fuzziness degree}

\item{alpha}{A float representing the weight of the space in the analysis (0
is a typical fuzzy-c-mean algorithm, 1 is balanced between the two
dimensions, 2 is twice the weight for space)}
}
\value{
A n X k matrix representing the belonging probabilities of each
  observation to each cluster
}
\description{
Calculate the new centroids of the clusters based on the membership matrix for SFCM
}
\keyword{internal}
