# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Indices of clustering quality #####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Fuzzy Silhouette index
#' @description Calculate the Silhouette index of clustering quality.
#' @details
#' The index is calculated with the function SIL.F from the package fclust.
#' When the dataset is too large, an approach by subsampling is used to avoid
#' crash.
#' @param data The original dataframe used for the clustering (n*p)
#' @param belongings A membership matrix (n*k)
#' @return A float, the fuzzy Silhouette index
calcSilhouetteIdx <- function(data, belongings){ #nocov start
  FS <- tryCatch(
    fclust::SIL.F(data, belongings, alpha = 1),
    error = function(e){
        warning("impossible to calculate Silhouette index with fclust::SIL.F, This is
 most likely due to a large dataset. We use here an approximation by subsampling...")
        FS <- stats::median(sapply(1:10, function(i){
          ids <- sample(1:nrow(data), size = 1000)
          fclust::SIL.F(data[ids,], belongings[ids,], alpha = 1)
        }))
        return(FS)
      }
  )
  return(FS)
} #nocov end
# not tested, come simply from fclust


#' @title Negentropy Increment index
#'
#' @description Calculate the Negentropy Increment index of clustering quality.
#'
#' @details
#' The Negentropy Increment index \insertCite{da2020incremental}{geocmeans} is based on the assumption that a normally shaped cluster is more
#' desirable. It uses the difference between the average negentropy
#' of all the clusters in the partition, and that of the  whole partition.
#' A smaller value indicates a better partition. The formula is:
#'
#' \deqn{NI=\frac{1}{2} \sum_{j=1}^{k} p_{i} \ln \left|{\boldsymbol{\Sigma}}_{j}\right|-\frac{1}{2} \ln \left|\boldsymbol{\Sigma}_{d a t a}\right|-\sum_{j=1}^{k} p_{j} \ln p_{j}}
#'
#' with  a cluster, \emph{|.|} the determinant of a matrix,
#' \itemize{
#'  \item \emph{j} a cluster
#'  \item \emph{|.|} the determinant of a matrix
#'  \item \eqn{\left|{\boldsymbol{\Sigma}}_{j}\right|} the covariance matrix of the dataset weighted by the membership values to cluster \emph{j}
#'  \item \eqn{\left|\boldsymbol{\Sigma}_{d a t a}\right|} the covariance matrix of the dataset
#'  \item \eqn{p_{j}} the sum of the membership values to cluster \emph{j} divided by the number of observations.
#' }
#' @references
#' \insertAllCited{}
#'
#' @param data The original dataframe used for the clustering (n*p)
#' @param belongmatrix A membership matrix (n*k)
#' @param centers The centres of the clusters
#' @return A float: the Negentropy Increment index
#' @export
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' calcNegentropyI(result$Data, result$Belongings, result$Centers)
calcNegentropyI <- function(data, belongmatrix, centers){
  # https://ieeexplore.ieee.org/document/8970493
  data <- as.matrix(data)
  n <- nrow(belongmatrix)
  ni <- colSums(belongmatrix)
  pi <- ni / n

  # calculer les matrices de covariances
  Si <- sapply(1:ncol(belongmatrix), function(i){
    det(stats::cov.wt(data, wt = belongmatrix[,i])$cov)
  })
  Sdata <- stats::cov(data)

  p1 <- (sum(pi * log(Si))) / 2
  p2 <- log(det(Sdata)) / 2
  p3 <- sum(pi * log(pi))
  return(p1 - p2 - p3)

}


#' @title Generalized Dunn’s index (53)
#'
#' @description Calculate the Generalized Dunn’s index (v53) of clustering quality.
#'
#' @details
#' The Generalized Dunn’s index  \insertCite{da2020incremental}{geocmeans} is a
#' ratio of the worst pair-wise separation of clusters and the worst compactness
#' of clusters. A higher value indicates a better clustering. The formula
#' is:
#'
#' \deqn{GD_{r s}=\frac{\min_{i \neq j}\left[\delta_{r}\left(\omega_{i}, \omega_{j}\right)\right]}{\max_{k}\left[\Delta_{s}\left(\omega_{k}\right)\right]}}
#'
#' The numerator is a measure of the minimal separation between all the clusters
#' \emph{i} and \emph{j} given by the formula:
#'
#' \deqn{\delta_{r}\left(\omega_{i}, \omega_{j}\right)=\frac{\sum_{l=1}^{n}\left\|\boldsymbol{x_{l}}-\boldsymbol{c_{i}}\right\|^{\frac{1}{2}} . u_{il}+\sum_{l=1}^{n}\left\|\boldsymbol{x_{l}}-\boldsymbol{c_{j}}\right\|^{\frac{1}{2}} . u_{jl}}{\sum{u_{i}} + \sum{u_{j}}}}
#'
#' where \emph{u} is the membership matrix and \eqn{u_{i}} is the column of
#' \emph{u} describing the membership of the \emph{n} observations to cluster
#' \emph{i}. \eqn{c_{i}} is the center of the cluster \emph{i}.
#'
#' The denominator is a measure of the maximal dispersion of all clusters, given
#' by the formula:
#'
#' \deqn{\frac{2*\sum_{l=1}^{n}\left\|\boldsymbol{x}_{l}-\boldsymbol{c_{i}}\right\|^{\frac{1}{2}}}{\sum{u_{i}}}}
#'
#'@references
#' \insertAllCited{}
#'
#' @param data The original dataframe used for the clustering (n*p)
#' @param belongmatrix A membership matrix (n*k)
#' @param centers The centres of the clusters
#' @return A float: the  Generalized Dunn’s index (53)
#' @export
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' calcGD53(result$Data, result$Belongings, result$Centers)
calcGD53 <- function(data, belongmatrix, centers){
  # https://ieeexplore.ieee.org/document/8970493
  data <- as.matrix(data)

  # calculate the separation based on separated compactness
  cols <- 1:ncol(belongmatrix)
  cps <- sapply(cols, function(i){
    ct <- centers[i,]
    sum(sweep(data, 2, ct, `-`)**2 * belongmatrix[,i]) ** (1/2)
  })
  ns <- colSums(belongmatrix)
  seps <- do.call(rbind,lapply(cols, function(i){
    others <- cols[cols!=i]
    values <- sapply(others, function(j){
      Mij <- (cps[[i]] + cps[[j]]) / (ns[[i]] + ns[[j]])
      return(Mij)
    })
  }))

  num <- min(seps)

  # calculate the compactness
  compacts <- sapply(cols, function(i){
    ct <- centers[i,]
    cp <- sum(sweep(data, 2, ct, `-`)**2 * belongmatrix[,i]) ** (1/2)
    2*cp / sum(belongmatrix[,i])
  })
  denom <- max(compacts)

  return(num/denom)
}


#' @title Generalized Dunn’s index (43)
#'
#' @description Calculate the Generalized Dunn’s index (v43) of clustering quality.
#'
#' @details
#' The Generalized Dunn’s index  \insertCite{da2020incremental}{geocmeans} is a
#' ratio of the worst pair-wise separation of clusters and the worst compactness
#' of clusters. A higher value indicates a better clustering. The formula
#' is:
#'
#' \deqn{GD_{r s}=\frac{\min_{i \neq j}\left[\delta_{r}\left(\omega_{i}, \omega_{j}\right)\right]}{\max_{k}\left[\Delta_{s}\left(\omega_{k}\right)\right]}}
#'
#' The numerator is a measure of the minimal separation between all the clusters
#' \emph{i} and \emph{j} given by the formula:
#'
#' \deqn{\delta_{r}\left(\omega_{i}, \omega_{j}\right)=\left\|\boldsymbol{c}_{i}-\boldsymbol{c}_{j}\right\|}
#'
#' which is basically the Euclidean distance between the centres of clusters \eqn{c_{i}} and \eqn{c_{j}}
#'
#' The denominator is a measure of the maximal dispersion of all clusters, given
#' by the formula:
#'
#' \deqn{\frac{2*\sum_{l=1}^{n}\left\|\boldsymbol{x}_{l}-\boldsymbol{c_{i}}\right\|^{\frac{1}{2}}}{\sum{u_{i}}}}
#'
#'@references
#' \insertAllCited{}
#'
#' @param data The original dataframe used for the clustering (n*p)
#' @param belongmatrix A membership matrix (n*k)
#' @param centers The centres of the clusters
#' @return A float: the  Generalized Dunn’s index (43)
#' @export
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' calcGD43(result$Data, result$Belongings, result$Centers)
calcGD43 <- function(data, belongmatrix, centers){
  # https://ieeexplore.ieee.org/document/8970493
  data <- as.matrix(data)

  # calculate the separation based on euclidean distance
  cols <- 1:ncol(belongmatrix)
  seps <- do.call(rbind,lapply(cols, function(i){
    others <- cols[cols!=i]
    values <- sapply(others, function(j){
      Mij <- sum((centers[i,] - centers[j,])**2) ** (1/2)
      return(Mij)
    })
  }))
  num <- min(seps)

  # calculate the compactness
  compacts <- sapply(cols, function(i){
    ct <- centers[i,]
    cp <- sum(sweep(data, 2, ct, `-`)**2 * belongmatrix[,i]) ** (1/2)
    2*cp / sum(belongmatrix[,i])
  })
  denom <- max(compacts)

  return(num/denom)
}



#' @title Davies-Bouldin index
#'
#' @description Calculate the Davies-Bouldin index of clustering quality.
#'
#' @details
#' The Davies-Bouldin index \insertCite{da2020incremental}{geocmeans} can be seen as the ratio of the within cluster dispersion and the
#' between cluster separation. A lower value indicates a higher cluster compacity
#' or a higher cluster separation. The formula is:
#'
#' \deqn{DB = \frac{1}{k}\sum_{i=1}^k{R_{i}}}
#'
#' with:
#'
#' \deqn{R_{i} =\max_{i \neq j}\left(\frac{S_{i}+S_{j}}{M_{i, j}}\right)}
#' \deqn{S_{l} =\left[\frac{1}{n_{l}} \sum_{l=1}^{n}\left\|\boldsymbol{x_{l}}-\boldsymbol{c_{i}}\right\|*u_{i}\right]^{\frac{1}{2}}}
#' \deqn{M_{i, j} =\sum\left\|\boldsymbol{c}_{i}-\boldsymbol{c}_{j}\right\|}
#'
#' So, the value of the index is an average of \eqn{R_{i}} values. For each cluster, they represent
#' its worst comparison with all the other clusters, calculated
#' as the ratio between the compactness of the two clusters and the separation
#' of the two clusters.
#'
#'@references
#' \insertAllCited{}
#'
#' @param data The original dataframe used for the clustering (n*p)
#' @param belongmatrix A membership matrix (n*k)
#' @param centers The centres of the clusters
#' @return A float: the Davies-Bouldin index
#' @export
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' calcDaviesBouldin(result$Data, result$Belongings, result$Centers)
calcDaviesBouldin <- function(data, belongmatrix, centers){
  # using this formula : https://en.wikipedia.org/wiki/Davies%E2%80%93Bouldin_index
  c <- apply(data,2,mean)
  data <- as.matrix(data)

  Si <- sapply(1:ncol(belongmatrix), function(i){
    ct <- centers[i,]
    (sum((sweep(data, 2, ct, `-`))**2) / sum(belongmatrix[,i]))**(1/2)
  })

  cols <- 1:ncol(belongmatrix)
  Rij <- do.call(rbind,lapply(cols, function(i){
    others <- cols[cols!=i]
    S <- Si[[i]]
    values <- sapply(others, function(j){
      Sj <- Si[[i]]
      Mij <- sum((centers[i,] - centers[j,])**2) ** (1/2)
      return( (S + Sj)/Mij )
    })
  }))

  Di <- apply(Rij, 1, max)
  DB <- mean(Di)

  return(DB)
}


#' @title Calinski-Harabasz index
#'
#' @description Calculate the Calinski-Harabasz index of clustering quality.
#'
#' @details
#' The Calinski-Harabasz index \insertCite{da2020incremental}{geocmeans} is the ratio between the clusters separation (between groups sum of squares) and the clusters cohesion (within groups sum of squares). A greater
#' value indicates either more separated clusters or more cohesive clusters.
#'
#'@references
#' \insertAllCited{}
#'
#' @param data The original dataframe used for the clustering (n*p)
#' @param belongmatrix A membership matrix (n*k)
#' @param centers The centres of the clusters
#' @return A float: the Calinski-Harabasz index
#' @export
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' calcCalinskiHarabasz(result$Data, result$Belongings, result$Centers)
calcCalinskiHarabasz <- function(data, belongmatrix, centers){
  # using this formula : https://www.geeksforgeeks.org/calinski-harabasz-index-cluster-validity-indices-set-3/
  data <- as.matrix(data)
  c <- apply(data,2,mean)
  k <- ncol(belongmatrix)
  nk <- colSums(belongmatrix)
  p1 <- (sum(rowSums((centers - c) **2) * nk)) / (k-1)
  p2 <- sum(apply(centers, 1, function(ck){
    sum((sweep(data, 2, ck, `-`))**2)
  })) / (nrow(data)-k)
  return(p1/p2)
}


#' @title Fukuyama and Sugeno index
#'
#' @description Calculate Fukuyama and Sugeno index of clustering quality
#'
#' @details
#' The Fukuyama and Sugeno index \insertCite{fukuyama1989new}{geocmeans} is the difference between the compacity of clusters and the separation of clusters. A smaller value indicates a better clustering.
#' The formula is:
#'
#' \deqn{S(c)=\sum_{k=1}^{n} \sum_{i=1}^{c}\left(U_{i k}\right)^{m}\left(\left\|x_{k}-v_{i}\right\|^{2}-\left\|v_{i}-\bar{x}\right\|^{2}\right) 2}
#'
#' with \emph{n} the number of observations, \emph{k} the number of clusters and \eqn{\bar{x}} the mean of the dataset.
#'
#' @references
#' \insertAllCited{}
#'
#' @param data The original dataframe used for the clustering (n*p)
#' @param belongmatrix A membership matrix (n*k)
#' @param centers The centres of the clusters
#' @param m The fuzziness parameter
#' @return A float: the Fukuyama and Sugeno index
#' @export
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' calcFukuyamaSugeno(result$Data,result$Belongings, result$Centers, 1.5)
calcFukuyamaSugeno <- function(data,belongmatrix,centers,m){
  v_hat <- colMeans(data)
  um <- (belongmatrix)**m
  values <- sapply(1:ncol(belongmatrix),function(i){
    w <- um[,i]
    v <- centers[i,]
    t1 <- calcEuclideanDistance2(data,v)
    t2 <- sum((v - v_hat)**2)
    dists <- (t1 - t2) * w * 2
    return(sum(dists))
  })
  return(sum(values))
}



#' @title Explained inertia index
#'
#' @description Calculate the explained inertia by a classification
#'
#' @param data The original dataframe used for the classification (n*p)
#' @param belongmatrix A membership matrix (n*k)
#' @return A float: the percentage of the total inertia explained
#' @export
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' calcexplainedInertia(result$Data,result$Belongings)
calcexplainedInertia <- function(data,belongmatrix){
  #step1 : calculating the centers
  centers <- t(apply(belongmatrix, MARGIN=2,function(column){
    values <- apply(data,MARGIN=2,function(var){
      return(weighted.mean(var,column))
    })
    return(values)
  }))
  means <- apply(data, 2, mean)
  baseinertia <- sum(calcEuclideanDistance2(data, means))
  restinertia <- sapply(1:nrow(centers), function(i) {
    point <- centers[i, ]
    dists <- calcEuclideanDistance2(data, point) * belongmatrix[, i]
    return(sum(dists))
  })
  explainedinertia <- 1 - (sum(restinertia) / baseinertia)
  return(explainedinertia)
}


#' @title calculate the quality index required
#'
#' @description A selector function to get the right quality index
#'
#' @param name The name of the index to calculate
#' @param ... The parameters needed to calculate the index
#' @return A float: the value of the index
#' @keywords internal
#' @examples
#' # this is an internal function, no example provided
calcQualIdx <- function(name, ...){
  dots <- list(...)
  val <- NULL
  if(name == "Silhouette.index"){
    #val <- fclust::SIL.F(dots$data, dots$belongmatrix, alpha=1)
    val <- calcSilhouetteIdx(dots$data, dots$belongmatrix)
  }

  if(name == "Partition.entropy"){
    val <- fclust::PE(dots$belongmatrix)
  }

  if(name == "Partition.coeff"){
    val <- fclust::PC(dots$belongmatrix)
  }

  if(name == "Modified.partition.coeff"){
    val <- fclust::MPC(dots$belongmatrix)
  }

  if(name == "XieBeni.index"){
    if(dots$m > 1){
      val <- fclust::XB(dots$data, dots$belongmatrix, dots$centers, dots$m)
    }else{
      val <- NA
    }
  }
  if(name == "FukuyamaSugeno.index"){
    val <- calcFukuyamaSugeno(dots$data, dots$belongmatrix, dots$centers, dots$m)
  }

  if(name == "CalinskiHarabasz.index"){
    val <- calcCalinskiHarabasz(dots$data, dots$belongmatrix, dots$centers)
  }

  if(name == "DaviesBoulin.index"){
    val <- calcDaviesBouldin(dots$data, dots$belongmatrix, dots$centers)
  }

  if(name == "Explained.inertia"){
    val <- calcexplainedInertia(dots$data, dots$belongmatrix)
  }

  if(name == "GD43.index"){
    val <- calcGD43(dots$data, dots$belongmatrix, dots$centers)
  }

  if(name == "GD53.index"){
    val <- calcGD53(dots$data, dots$belongmatrix, dots$centers)
  }

  if(name == "Negentropy.index"){
    val <- calcNegentropyI(dots$data, dots$belongmatrix, dots$centers)
  }

  if(is.null(val)){
    stop('The index name must be one of "Silhouette.index", "Partition.entropy",
    "Partition.coeff", "XieBeni.index", "FukuyamaSugeno.index","Explained.inertia",
    "DaviesBoulin.index", "CalinskiHarabasz.index", "GD43.index", "GD53.index",
    "Negentropy.index"')
  }
  return(val)
}



#' @title Quality indexes
#'
#' @description calculate several clustering quality indexes (some of them come from fclust
#' package)
#'
#' @param data The original dataframe used for the classification (n*p)
#' @param belongmatrix A membership matrix (n*k)
#' @param m The fuzziness parameter used for the classification
#' @param indices A character vector with the names of the indices to calculate, default is :
#' c("Silhouette.index", "Partition.entropy", "Partition.coeff", "XieBeni.index", "FukuyamaSugeno.index",
#' "Explained.inertia"). Other available indices are : "DaviesBoulin.index", "CalinskiHarabasz.index",
#' "GD43.index", "GD53.index" and "Negentropy.index"
#' @return A named list with with the values of the required indices
#' @export
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' calcqualityIndexes(result$Data,result$Belongings, m=1.5)
calcqualityIndexes <- function(data, belongmatrix, m, indices = c("Silhouette.index", "Partition.entropy", "Partition.coeff", "XieBeni.index", "FukuyamaSugeno.index", "Explained.inertia")) {

  centers <- t(apply(belongmatrix, MARGIN=2,function(column){
    values <- apply(data,MARGIN=2,function(var){
      return(weighted.mean(var,column))
    })
    return(values)
  }))

  data <- as.matrix(data)

  idxs <- lapply(indices, function(name){
    val <- calcQualIdx(name, data = data, belongmatrix = belongmatrix,
                       centers = centers, m = m)
  })
  names(idxs) <- indices
  return(idxs)
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Indices of spatial clustering quality #####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Spatial diagnostic
#'
#' @description Utility function to facilitate the spatial diagnostic of a classification
#'
#' Calculate the following indicators: Moran I index (spdep::moranI) for each
#' column of the membership matrix, Join count test (spdep::joincount.multi) for
#' the most likely groups of each datapoint, Spatial consistency index (see
#' function spConsistency) and the Elsa statistic (see function calcElsa). Note
#' that if the FCMres object given was constructed with rasters, the joincount
#' statistic is not calculated and no p-values are provided for the Moran I
#' indices.
#'
#' @template FCMresobj-arg
#' @template nblistw2-arg
#' @param window If rasters were used for the classification, the window must be
#' specified instead of a list.w object. Can also be NULL if object is a FCMres object.
#' @param undecided A float giving the threslhod to detect undecided observations. An
#' observation is undecided if its maximum membership value is bellow this float. If
#' null, no observations are undecided.
#' @template matdist-arg
#' @param nrep An integer indicating the number of permutation to do to simulate
#'   the random distribution of the spatial inconsistency
#' @return A named list with :
#' \itemize{
#'         \item MoranValues : the moran I values for each column of the membership
#'          matrix (spdep::MoranI)
#'         \item JoinCounts : the result of the join count test calculated with
#'          the most likely group for each datapoint (spdep::joincount.multi)
#'         \item SpConsist : the mean value of the spatial consistency index
#'         (the lower, the better, see ?spConsistency for details)
#'}
#' @export
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' spatialDiag(result, undecided=0.45, nrep=30)
spatialDiag <- function(object, nblistw = NULL, window = NULL, undecided = NULL, matdist = NULL, nrep = 50) { #nocov start

  cls <- class(object)[[1]]
  ## prior check of parameters
  if(cls != "FCMres") {
    if(is.null(nblistw)){
      stop("if object is not a FCMres object, nblistw must be provided")
    }
    if(is.null(matdist)){
      warning("if object is not a FCMres object, matdist must be provided, otherwise ELSA will not be calculated")
    }else{
      check_matdist(matdist)
    }
    if(is.null(window) == FALSE){
      warning("When object is not a FCMres object, it is assumed that the data used are not rasters, window is not used in that context")
    }
  }else{
    if(object$isRaster){
      if(is.null(object$window) & is.null(window)){
        stop("The FCMres object was created with rasters, but does not contain a window slot and the window parameter is NULL. Please specify a window")
      }
    }else{
      if(is.null(object$nblistw) & is.null(nblistw)){
        stop("The FCMres object was created with a SpatialDataFrame or a dataframe, but does not contain a nblistw slot and the nblistw parameter is NULL. Please specify a nblistw")
      }
    }
    if(is.null(matdist)){
      matdist <- as.matrix(stats::dist(object$Centers))
    }

  }

  ## check if we are in rasterMode here is a mistake !
  rasterMode <- FALSE
  if(cls == "FCMres"){
    if(object$isRaster){
      rasterMode <- TRUE
    }
  }

  ## WE ARE NOT IN RASTERMODE ##
  if(rasterMode == FALSE){
    if(cls != "FCMres"){
      membership <- object
      Groups <- (1:ncol(membership))[max.col(membership, ties.method = "first")]
    }else{
      membership <- object$Belongings
      Groups <- as.numeric(gsub("V","",object$Groups, fixed = TRUE))
      if(is.null(nblistw)){
        nblistw <- object$nblistw
      }
    }

    # calculating spconsistency
    Consist <- spConsistency(membership, nblistw = nblistw, nrep = nrep)
    # calculating ELSA
    if(is.null(matdist) == FALSE){
      Elsa <- calcELSA(Groups, nblistw = nblistw, matdist = matdist)
    }else{
      Elsa <- NULL
    }
    # calculating MORAN I
    Values <- sapply(1:ncol(membership), function(i) {
      x <- membership[, i]
      morantest <- spdep::moran.mc(x, nblistw, nsim = 999)
      return(list(MoranI = morantest$statistic, pvalue = morantest$p.value,
                  Cluster = paste("Cluster_", i, sep = "")))
    })
    morandf <- as.data.frame(t(Values))
    for(col in colnames(morandf)){
      morandf[[col]] <- unlist(morandf[[col]])
    }
    # calculating join count test
    Groups <- as.character(Groups)
    if (is.null(undecided) == FALSE) {
      Maximums <- do.call(pmax, as.data.frame(membership))
      Groups[Maximums < undecided] <- "undecided"
    }
    Groups <- as.factor(Groups)
    # calcul des join count test
    spjctetst <- spdep::joincount.multi(Groups, nblistw, zero.policy = TRUE)

    #returning the results
    return(list(MoranValues = morandf, JoinCounts = spjctetst,
                SpConsist = Consist$Mean, SpConsistSamples = Consist$samples,
                Elsa = Elsa))

  }else{
    ## WE ARE IN RASTERMODE ##
    if(is.null(window)){
      window <- object$window
    }
    ## calculating spconsistency
    Consist <- spConsistency(object, window = window, nrep = nrep)
    ## calculating ELSA
    if(is.null(matdist) == FALSE){
      Elsa <- calcELSA(object, window = window, matdist = matdist)
    }else{
      Elsa <- NULL
    }
    ## calculating Moran I
    MoranI_vals <- sapply(object$rasters[1:object$k], function(rast){
      calc_moran_raster(rast,window)
    })
    morandf <- data.frame(
      Cluster = paste("Cluster_", 1:length(MoranI_vals), sep = ""),
      MoranI = MoranI_vals
    )
    ## returning the results
    return(list(MoranValues = morandf, SpConsist = Consist$Mean,
                SpConsistSamples = Consist$samples, Elsa = Elsa))
  }
} #nocov end
