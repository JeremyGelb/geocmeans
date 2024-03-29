% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/side_functions.R
\name{adjustSpatialWeights}
\alias{adjustSpatialWeights}
\title{Semantic adjusted spatial weights}
\usage{
adjustSpatialWeights(data, listw, style, mindist = 1e-11)
}
\arguments{
\item{data}{A dataframe with numeric columns}

\item{listw}{A nb object from spdep}

\item{style}{A letter indicating the weighting scheme (see spdep doc)}

\item{mindist}{A minimum value for distance between two observations. If two
neighbours have exactly the same values, then the euclidean distance between
them is 0, leading to an infinite spatial weight. In that case, the minimum
distance is used instead of 0.}
}
\value{
A listw object (spdep like)
}
\description{
Function to adjust the spatial weights so that they represent semantic
distances between neighbours
}
\examples{
data(LyonIris)
AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
"TxChom1564","Pct_brevet","NivVieMed")
dataset <- sf::st_drop_geometry(LyonIris[AnalysisFields])
queen <- spdep::poly2nb(LyonIris,queen=TRUE)
Wqueen <- spdep::nb2listw(queen,style="W")
Wqueen2 <- adjustSpatialWeights(dataset,queen,style="C")
}
