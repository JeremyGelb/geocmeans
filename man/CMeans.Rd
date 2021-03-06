% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main_functions.R
\name{CMeans}
\alias{CMeans}
\title{C-means}
\usage{
CMeans(
  data,
  k,
  m,
  maxiter = 500,
  tol = 0.01,
  standardize = TRUE,
  verbose = TRUE,
  seed = NULL
)
}
\arguments{
\item{data}{A dataframe with only numerical variable}

\item{k}{An integer describing the number of cluster to find}

\item{m}{A float for the fuzziness degree}

\item{maxiter}{A float for the maximum number of iteration}

\item{tol}{The tolerance criterion used in the evaluateMatrices function for
convergence assessment}

\item{standardize}{A boolean to specify if the variables must be centered and
reduce (default = True)}

\item{verbose}{A boolean to specify if the messages should be displayed}

\item{seed}{An integer used for random number generation. It ensures that the
start centers will be the same if the same integer is selected.}
}
\value{
A named list with :
 \itemize{
        \item Centers: a dataframe describing the final centers of the groups
        \item Belongings: the final belonging matrix
        \item Groups: a vector with the names of the most likely group for each observation
        \item Data: the dataset used to perform the clustering (might be standardized)
}
}
\description{
The clasical c-mean algorithm
}
\examples{
data(LyonIris)
AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
"TxChom1564","Pct_brevet","NivVieMed")
dataset <- LyonIris@data[AnalysisFields]
result <- CMeans(dataset,k = 5, m = 1.5, standardize = TRUE)
}
