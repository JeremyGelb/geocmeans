% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FCMclass.R
\name{is.FCMres}
\alias{is.FCMres}
\title{is method for FCMres}
\usage{
\method{is}{FCMres}(object, class2 = "FCMres")
}
\arguments{
\item{object}{A FCMres object, typically obtained from functions CMeans, GCMeans, SFCMeans, SGFCMeans}

\item{class2}{Character string giving the names of the classe to test (usually "FCMres")}
}
\value{
A boolean, TRUE if x can be considered as a FCMres object, FALSE otherwise
  group
}
\description{
Check if an object can be considered as a FCMres object
}
\examples{
data(LyonIris)
AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
"TxChom1564","Pct_brevet","NivVieMed")
dataset <- sf::st_drop_geometry(LyonIris[AnalysisFields])
queen <- spdep::poly2nb(LyonIris,queen=TRUE)
Wqueen <- spdep::nb2listw(queen,style="W")
result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
is(result, "FCMres")
}
