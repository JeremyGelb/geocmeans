#' @title Undecided observations
#'
#' @description Identify the observation for which the classification is uncertain
#'
#' @param belongmatrix The membership matrix obtained at the end of the algorithm
#' @param tol A float indicating the minimum required level of membership to be
#'   not considered as undecided
#' @param out The format of the output vector. Default is "character". If
#' "numeric", then the undecided units are set to -1.
#' @return A vector indicating the most likely group for each observation or
#'   "Undecided" if the maximum probability for the observation does not reach
#'   the value of the tol parameter
#' @export
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' undecidedUnits(result$Belongings, tol = 0.45)
undecidedUnits <- function(belongmatrix, tol = 0.1, out = "character") {
  belongmatrix <- as.data.frame(belongmatrix)
  if(out == "character"){
    choose_from <- colnames(belongmatrix)
    alt <- "Undecided"
  }else if(out == "numeric"){
    choose_from <- 1:ncol(belongmatrix)
    alt <- -1
  }else{
    stop("The argument out must be one of c('numeric', 'character')")
  }
  groups <- choose_from[max.col(belongmatrix, ties.method = "first")]
  rowMax <- do.call(pmax, belongmatrix)
  DF <- data.frame(groups = groups, maxprob = rowMax)
  return(ifelse(DF$maxprob < tol, alt , DF$groups))
}


#' @title Diversity index
#'
#' @description Calculate the diversity (or entropy) index.
#'
#' @details
#' The diversity (or entropy) index \insertCite{theil1972statistical}{geocmeans} is calculated for each observation an vary between 0 and 1. When the value is close to 0, the
#' observation belong to only one cluster (as in hard clustering). When the value is
#' close to 1, the observation is undecided and tends to belong to each cluster. Values
#' above 0.9 should be investigated. The formula is:
#'
#' \deqn{H2_{i} = \frac{-\sum[u_{ij}\ln(u_{ij})]}{\ln(k)}.}
#'
#' with \emph{i} and observation, \emph{j} a cluster, \emph{k} the number of clusters and
#' \emph{u} the membership matrix.
#'
#' It is a simplified formula because the sum of each row of a membership matrix
#' is 1.
#'
#' @references
#' \insertAllCited{}
#'
#' @param belongmatrix The membership matrix obtained at the end of the algorithm
#' @return A vector with the values of the diversity (entropy) index
#' @importFrom Rdpack reprompt
#' @export
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' calcUncertaintyIndex(result$Belongings)
calcUncertaintyIndex <- function(belongmatrix){
  k <- ncol(belongmatrix)
  return(-1 * rowSums(belongmatrix * log(belongmatrix)) / log(k))
}

