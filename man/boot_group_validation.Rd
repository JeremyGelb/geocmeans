% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boostrap_clust_validation.R
\name{boot_group_validation}
\alias{boot_group_validation}
\title{Bootstrap check the robustness of a classification}
\usage{
boot_group_validation(
  object,
  nsim = 1000,
  maxiter = 1000,
  tol = 0.01,
  init = "random",
  verbose = TRUE,
  seed = NULL
)
}
\arguments{
\item{object}{A FCMres object, typically obtained from functions CMeans,
GCMeans, SFCMeans, SGFCMeans}

\item{nsim}{The number of replications to do for the bootstrap evaluation}

\item{maxiter}{An integer for the maximum number of iteration}

\item{tol}{The tolerance criterion used in the evaluateMatrices function for
convergence assessment}

\item{init}{A string indicating how the initial centers must be selected.
"random" indicates that random observations are used as centers. "kpp" use
a distance based method resulting in more dispersed centers at the
beginning. Both of them are heuristic.}

\item{verbose}{A boolean to specify if the progress bar should be displayed.}

\item{seed}{An integer to control randomness, default is NULL}
}
\value{
A list of two values: group_consistency, a dataframe indicating for
  each cluster its consistency across simulations. group_centers a list with
  a dataframe for each cluster. The values in the dataframes are the centers
  of the clusters at each simulation.
}
\description{
Check that the obtained groups are stable by bootstrap
}
\details{
Considering that the classification produced by a FCM like algorithm
  depends on its initial state, it is important to check if the groups
  obtained are stable. This function use a bootstrap method to do so. During
  a selected number of iterations (at least 1000), a sample of size n (with
  replacement) is drawn from the original dataset. For each sample, the same
  classification algorithm is applied and the results are compared with the
  reference results. For each original group, the mots similar group is
  identified by calculating the Jaccard similarity index between the columns
  of the two membership matrices. This index is comprised between 0 (exact
  difference) and 1 (perfect similarity) and a value is calculated for each
  group at each iteration. One can investigate the values obtained to
  determine if the groups are stable. Values under 0.5 are worrisome and
  indicate that the group is dissolving. Values between 0.6 and 0.75 indicate
  a pattern in the data, but an important uncertainty. Values above 0.8
  indicate strong groups. The values of the centers obtained at each
  iteration are also returned, it is important to check the they follow
  approximately a normal distribution and are at least unimodal.
}
\examples{
\dontrun{
data(LyonIris)

#selecting the columns for the analysis
AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14",
                   "Pct_65","Pct_Img","TxChom1564","Pct_brevet","NivVieMed")

#rescaling the columns
Data <- LyonIris@data[AnalysisFields]
for (Col in names(Data)){
  Data[[Col]] <- as.numeric(scale(Data[[Col]]))
}

Cmean <- CMeans(Data,4,1.5,500,standardize = FALSE, seed = 456,
    tol = 0.00001, verbose = FALSE)

validation <- boot_group_validation(Cmean, nsim = 1000, maxiter = 1000,
    tol = 0.01, init = "random")
}
}