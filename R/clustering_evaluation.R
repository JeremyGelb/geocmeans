# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Indices of clustering quality #####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# a version of the silhouette index when the dataset is really too big.
calcSilhouetteIdx <- function(data, belongings){
  FS <- tryCatch(
    fclust::SIL.F(data, belongings, alpha = 1),
    error = function(e){
        warning("impossible to calculate Silhouette index with fclust::SIL.F, This is
 most likely due to a large dataset. We use here an approximation by subsampling...")
        FS <- median(sapply(1:10, function(i){
          ids <- sample(1:nrow(data), size = 1000)
          fclust::SIL.F(data[ids,], belongings[ids,], alpha = 1)
        }))
        return(FS)
      }
  )
  return(FS)
}


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
#' @param centers The centers of the clusters
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
#' @param centers The centers of the clusters
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
#' which is basically the Euclidean distance between the centers of clusters \eqn{c_{i}} and \eqn{c_{j}}
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
#' @param centers The centers of the clusters
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
#' The Davies-Bouldin index \insertCite{da2020incremental}{geocmeans} can be seen as a ratio between the within cluster dispersion and the
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
#' So, the value of the index is an average of \eqn{R_{i}} values. They represent
#' for each cluster its worst comparison with all the other clusters, calculated
#' as the ratio between the compactness of the two clusters and the separation
#' of the two clusters.
#'
#'@references
#' \insertAllCited{}
#'
#' @param data The original dataframe used for the clustering (n*p)
#' @param belongmatrix A membership matrix (n*k)
#' @param centers The centers of the clusters
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
#' @param centers The centers of the clusters
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
#' \deqn{J_{m}=\sum_{j=1}^{n} \sum_{i=1}^{k}\left(\mu_{i j}\right)^{m}\left\|x_{j}-c_{i}\right\|^{2}-\sum_{i=1}^{k}\left\|c_{i}-\bar{c}\right\|^{2}}
#'
#' with \emph{n} the number of observations, \emph{k} the number of clusters and \eqn{\bar{c}} the mean of the centers of the clusters.
#'
#' @references
#' \insertAllCited{}
#'
#' @param data The original dataframe used for the clustering (n*p)
#' @param belongmatrix A membership matrix (n*k)
#' @param centers The centers of the clusters
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
  v_hat <- apply(centers,2,mean)
  um <- (belongmatrix)**m
  values <- sapply(1:ncol(belongmatrix),function(i){
    w <- um[,i]
    v <- centers[i,]
    t1 <- calcEuclideanDistance2(data,v)
    t2 <- sum((v - v_hat)**2)
    dists <- (t1 - t2)*w
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
    val <- calcDaviesBouldin(dots$data, dots$belongmatrix, dots$centers)
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
#' Calculate the following indicators : Moran I index (spdep::moranI) for each
#' column of the membership matrix, Join count test (spdep::joincount.multi) for
#' the most likely groups of each datapoint, Spatial consistency index (see
#' function spConsistency) and the Elsa statistic (see function calcElsa). Note
#' that if the FCMres object given was constructed with rasters, the joincount
#' statistic is not calculated and no p-values are provided for the Moran I
#' indices.
#'
#' @template FCMresobj_arg
#' @template nblistw-arg. Can also be NULL if object is a FCMres object
#' @param window if rasters were used for the classification, the window must be
#' specified instead of a list.w object. Can also be NULL if object is a FCMres object.
#' @param centers If object is not a FCMres object, then centers of the clusters must
#' be given. A hard partition will be created by using the membership matrix given as
#' obejct and the ELSA statistic will be calculate.
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
spatialDiag <- function(object, nblistw = NULL, window = NULL, undecided = NULL, centers = NULL, nrep = 50) {

  ## prior check of parameters
  if(class(object)[[1]] != "FCMres") {
    if(is.null(nblistw)){
      stop("if object is not a FCMres object, nblistw must be provided")
    }
    if(is.null(window) == FALSE){
      stop("When object is not a FCMres object, window can not be specified by hand")
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

  }

  ## check if we are in rasterMode
  rasterMode <- FALSE
  if(class(object)[[1]] == "FCMres"){
    rasterMode <- object$isRaster
  }

  # calcul de la coherence spatiale et de ELSA
  if(is.null(nblistw == FALSE)){
    Consist <- spConsistency(object, nblistw = nblistw, nrep = nrep)
    if(class(object)[[1]] == "FCMres"){
      Elsa <- calcELSA(object, nblistw = nblistw)
    }else{
      if(is.null(centers) == FALSE){
        Groups <- (1:ncol(object))[max.col(object, ties.method = "first")]
        Elsa <- calcELSA(Groups, nblistw = nblistw, centers = centers)
      }else{
        warning("impossible to calculate ELSA if object is not from FCMres class and centers are not given.")
        Elsa <- NULL
      }
    }
  }else if (is.null(window) == FALSE){
    Consist <- spConsistency(object, window = window, nrep = nrep)
    Elsa <- calcELSA(object, w = window)
  }else{
    Consist <- spConsistency(object, nrep = nrep)
    Elsa <- calcELSA(object)
  }

  ## IF WE ARE NOT IN RASTERMODE
  if(rasterMode == FALSE){
    belongmatrix <- as.data.frame(belongmatrix)
    # calcul des I de Moran pour les appartenances
    Values <- sapply(1:ncol(belongmatrix), function(i) {
      x <- belongmatrix[, i]
      morantest <- spdep::moran.mc(x, nblistw, nsim = 999)
      return(list(MoranI = morantest$statistic, pvalue = morantest$p.value,
                  Cluster = paste("Cluster_", i, sep = "")))
    })
    morandf <- as.data.frame(t(Values))
    for(col in colnames(morandf)){
      morandf[[col]] <- unlist(morandf[[col]])
    }
    # attribution des groupes
    groups <- colnames(belongmatrix)[max.col(belongmatrix, ties.method = "first")]
    if (is.null(undecided) == FALSE) {
      Maximums <- do.call(pmax, belongmatrix)
      groups[Maximums < undecided] <- "undecided"
    }
    groups <- as.factor(groups)
    # calcul des join count test
    spjctetst <- spdep::joincount.multi(groups, nblistw, zero.policy = TRUE)
    return(list(MoranValues = morandf, JoinCounts = spjctetst,
                SpConsist = Consist$Mean, SpConsistSamples = Consist$samples,
                Elsa = Elsa))

  }else{
    # IF WE ARE IN rasterMODE
    # calculating the morans I for membership values
    name <- names(object$rasters)
    ok_names <- name[grepl("group",name, fixed = TRUE)]

    MoranI_vals <- sapply(object$rasters[ok_names], function(rast){
      calc_moran_raster(rast,window)
    })
    morandf <- data.frame(
      Cluster = paste("Cluster_", 1:length(MoranI_vals), sep = ""),
      MoranI = MoranI_vals
    )
    return(list(MoranValues = morandf, SpConsist = Consist$Mean,
                SpConsistSamples = Consist$samples, Elsa = Elsa))

  }

}


#' @title Spatial consistency index
#'
#' @description Calculate a spatial consistency index
#'
#' @details This index is experimental, it aims to measure how much a clustering solution
#' is spatially consistent. A classification is spatially inconsistent if
#' neighbouring observation do not belong to the same group. See detail for
#' a description of its calculation
#'
#' The total spatial inconsistency (*Scr*) is calculated as follow
#'
#' \deqn{isp = \sum_{i}\sum_{j}\sum_{k} (u_{ik} - u_{jk})^{2} * W_{ij}}
#'
#' With U the membership matrix, i an observation, k the neighbours of i and W
#' the spatial weight matrix This represents the total spatial inconsistency of
#' the solution (true inconsistency) We propose to compare this total with
#' simulated values obtained by permutations (simulated inconsistency). The
#' values obtained by permutation are an approximation of the spatial
#' inconsistency obtained in a random context Ratios between the true
#' inconsistency and simulated inconsistencies are calculated A value of 0
#' depict a situation where all observations are identical to their neighbours
#' A value of 1 depict a situation where all observations are as much different
#' as their neighbours that what randomness can produce A classification
#' solution able to reduce this index has a better spatial consistency
#'
#' @template FCMresobj_arg
#' @template nblistw-arg. Can also be NULL if object is a FCMres object
#' @param window if rasters were used for the classification, the window must be
#' specified instead of a list.w object. Can also be NULL if object is a FCMres object.
#' @param nrep An integer indicating the number of permutation to do to simulate
#' spatial randomness. Note that if rasters are used, each permutation can be very long.
#' @return A named list with
#'  \itemize{
#'         \item Mean : the mean of the spatial consistency index
#'         \item prt05 : the 5th percentile of the spatial consistency index
#'         \item prt95 : the 95th percentile of the spatial consistency index
#'         \item samples : all the value of the spatial consistency index
#' }
#' @export
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' spConsistency(result$Belongings, nblistw = Wqueen, nrep=50)
spConsistency <- function(object, nblistw = NULL, window = NULL, nrep = 999) {

  if(class(object)[[1]] == "FCMres"){
    belongmat <- as.matrix(object$Belongings)
    if(object$isRaster & is.null(window)){
      window <- object$window
      if(is.null(window)){
        stop("impossible to find a window in the given object, please
             specify one by hand.")
      }
    }
  }else{
    belongmat <- as.matrix(object)
  }

  # if we are not in raster mode

  if(is.null(window)){

    if(is.null(nblistw)){
      stop("The nblistw must be provided if spatial vector data is used")
    }
    weights <- nblistw$weights
    neighbours <- nblistw$neighbours
    ## calcul de l'inconsistence spatiale actuelle
    obsdev <- sapply(1:nrow(belongmat), function(i) {
      row <- belongmat[i, ]
      idneighbour <- neighbours[[i]]
      neighbour <- belongmat[idneighbour, ]
      if (length(idneighbour) == 1){
        neighbour <- t(as.matrix(neighbour))
      }
      W <- weights[[i]]
      diff <- (neighbour-row[col(neighbour)])**2
      tot <- sum(rowSums(diff) * W)
      return(tot)
    })

    totalcons <- sum(obsdev)

    ## simulation de l'inconsistance spatiale
    belongmat <- t(belongmat)
    n <- ncol(belongmat)
    simulated <- vapply(1:nrep, function(d) {
      belong2 <- belongmat[,sample(n)]
      simvalues <- vapply(1:ncol(belong2), function(i) {
        row <- belong2[,i]
        idneighbour <- neighbours[[i]]
        neighbour <- belong2[,neighbours[[i]]]
        if (length(idneighbour) == 1){
          neighbour <- t(as.matrix(neighbour))
        }
        W <- weights[[i]]
        diff <- (neighbour-row)
        tot <- sum(diff^2 * W)
        return(tot)
      }, FUN.VALUE = 1)
      return(sum(simvalues))
    },FUN.VALUE = 1)
    ratio <- totalcons / simulated

  # if we are using a raster mode.
  }else{

    # we must calculate for each pixel its distance to its neighbours
    # on the membership matrix. So we will calculate the distance for each
    # raster in object$rasters and then sum them all group
    rastnames <- names(object$rasters)
    ok_names <- rastnames[grepl("group",rastnames, fixed = TRUE)]
    rasters <- object$rasters[ok_names]
    matrices <- lapply(rasters, raster::as.matrix)
    totalcons <- calc_raster_spinconsistency(matrices,window)

    # we must now do the same thing but with resampled values
    warning("Calculating the permutation for the spatial inconsistency
            when using raster can be long, depending on the raster size.
            Note that the high number of cell in a raster reduces the need of
            a great number of replications.")
    all_ids <- 1:raster::ncell(rasters[[1]])

    datavecs <- lapply(rasters, function(rast){
      mat <- raster::as.matrix(rast)
      dim(mat) <- NULL
      return(mat)
    })
    mat_dim <- dim(raster::as.matrix(rasters[[1]]))
    simulated <- sapply(1:nrep, function(i){
      Ids <- sample(all_ids)
      # resampling the matrices
      new_matrices <- lapply(datavecs, function(vec){
        new_vec <- vec[Ids];
        dim(new_vec) <- mat_dim
        return(new_vec)
      })

      # calculating the index value
      inconsist <-  calc_raster_spinconsistency(new_matrices, window)
      return(inconsist)
    })
    ratio <- totalcons / simulated
  }

  return(list(Mean = mean(ratio), Median = quantile(ratio, probs = c(0.5)),
              prt05 = quantile(ratio, probs = c(0.05)),
              prt95 = quantile(ratio, probs = c(0.95)),
              samples = ratio))
}


#' @title calculate ELSA statistic for a hard partition
#'
#' @description Calculate ELSA statistic for a hard partition. This local indicator of
#' spatial autocorrelation can be used to determine where observations belong to different
#' clusters.
#'
#' @details The ELSA index \insertCite{naimi2019elsa}{geocmeans} can be used to measure
#' local autocorrelation for a categorical variable. It varies between 0 and 1, 0 indicating
#' a perfect positive spatial autocorrelation and 1 a perfect heterogeneity. It is based on
#' the Shanon entropy index, and uses a measure of difference between categories. Thus it
#' can reflect that proximity of two similar categories is still a form of positive
#' autocorellation. The authors suggest to calculate the mean of the index at several lag
#' distance to create an entrogram which quantifies global spatial structure and
#' represents as a variogram-like graph.
#'
#' @param object A FCMres object, typically obtained from functions CMeans,
#'   GCMeans, SFCMeans, SGFCMeans. Can also be a vector of categories. This vector must
#'   be filled with integers starting from 1. -1 can be used to indicate missing categories.
#' @template nblistw-arg
#' @param w The size of the padding to create the window if working with raster data or a matrix
#' representing the window. This matrix must be binary (0 and 1).
#' @param centers A dataframe or a matrix. Each row must represent the center of a cluster.
#' An euclidean distance matrix will be calculated on them to weight the ELSA statistic.
#' @return A float: the sum of spatial inconsistency
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' elsa_valus <- calcELSA(result)
calcELSA <- function(object, nblistw = NULL, w = NULL, centers = NULL){

  # testing if we have all the required parameters
  if(class(object)[[1]] != "FCMres"){

    if(class(object) != "numeric"){
      stop("if object is not a FCMres object, it must be a numeric vector")
    }

    vec <- object[[object != -1]]
    if(min(vec) != 1){
      stop("if object is a numeric vector, its lower value must be 1")
    }

    if(length(unique(vec)) != length(1:max(vec))){
      stop("if object is a numeric vector, its values must be like 1,2,3,4,...m, -1 can be used to indicate missing values")
    }

    if(!isSymmetric(centers)){
      stop("centers must be a symetric matrix")
    }

    if(nrow(centers) != length(unique(vec))){
      stop("the dimension of centers must equal the number of categories in object.")
    }

    if(is.null(nblistw)){
      stop("if object is not a FCMres object, nblistw must be provided.")
    }
    if(is.null(centers)){
      stop("if object is not a FCMres object, centers must be provided.")
    }
  }else{
    if(object$isRaster){
      if(is.null(object$window) & is.null(w)){
        stop("impossible to extract window from object, w must be given.")
      }
    }else{
      if(is.null(object$nblistw) & is.null(nblistw)){
        stop("impossible to extract nblistw from object, nblistw must be given.")
      }
    }
  }

  # case 1 : object is a simple vector of categories
  if(class(object)[[1]] != "FCMres"){
    matdist <- as.matrix(dist(centers))
    vals <- elsa_vector(object, nblistw, matdist)
  }else{
  # case 2 : the object is a FCMres object
    if(is.null(centers)){
      matdist <- as.matrix(dist(object$Centers))
    }else{
      matdist <- as.matrix(dist(Centers))
    }

    if(object$isRaster){
      if(is.null(w)){
        w <- object$window
      }
      vals <- elsa_raster(object$rasters$Groups, w, matdist)
    }else{
      vec <- as.numeric(gsub("V","",object$Groups,fixed = TRUE))
      if(is.null(nblistw)){
        vals <- elsa_vector(vec, object$nblistw, matdist)
      }else{
        vals <- elsa_vector(vec, nblistw, matdist)
      }
    }
  }
  return(vals)
}


#' @title calculate spatial inconsistency for raster
#'
#' @description Calculate the spatial inconsistenct sum for a set of rasters
#'
#' @param matrices A list of matrices
#' @param window The window to use to define spatial neighbouring
#' @return A float: the sum of spatial inconsistency
#' @keywords internal
#' @examples
#' # this is an internal function, no example provided
calc_raster_spinconsistency <- function(matrices, window){
  totalcons <- sum(focal_euclidean_list(matrices,window), na.rm = T)
  return(totalcons)
}

#' @title calculate ELSA spatial statistic for vector dataset
#'
#' @description calculate ELSA spatial statistic for vector dataset
#'
#' @param categories An integer vector representing the m categories (1,2,3,..., m),
#' -1 is used to indicate missing values.
#' @param nblistw A listw object from spdep representing neighbour relations
#' @param dist A numeric matrix (m*m) representing the distances between categories
#' @return A vector: the local values of ELSA
#' @keywords internal
#' @examples
#' # this is an internal function, no example provided
elsa_vector <- function(categories, nblistw, dist){
  d <- max(dist)
  m <- length(unique(categories))
  values <- sapply(1:length(categories), function(i){

    xi <- categories[[i]]
    if(xi == -1){
      return(-1)
    }else{
      neighbours <- nblistw$neighbours[[i]]
      w <- nblistw$weights[[i]]
      xjs <- categories[neighbours]
      Eai <- sum(dist[xi,xjs] * w) / (d * sum(w))

      nn <- length(w) # the number of observations in the area...
      if(nn  > m){
        mi <- m
      }else{
        mi <- nn
      }

      probs <- table(c(xjs,xi)) / (length(xjs)+1)
      Eci <- -1 * (sum(probs * log2(probs)) / log2(mi))
      return(Eai * Eci)
    }

  })
  return(values)

}


#' @title calculate ELSA spatial statistic for raster dataset
#'
#' @description calculate ELSA spatial statistic for vector dataset
#'
#' @param rast An integer raster or matrix representing the m categories (0,1,2,..., m)
#' @param w The size of the window to apply, or a binary matrix to define a window.
#' @param dist A numeric matrix (m*m) representing the distances between categories
#' @return A raster or a matrix: the local values of ELSA
#' @keywords internal
#' @examples
#' # this is an internal function, no example provided
elsa_raster <- function(rast, w, dist){

  if(class(w)[[1]] == "matrix"){
    u <- unique(w)
    if(length(union(c(0,1),u)) > 2){
      stop("The provided matrix to calculate ELSA is not a binary matrix (0,1)")
    }
    fun <- Elsa_categorical_matrix_window
  }else if (class(w)[[1]] == "numeric"){
    fun <- Elsa_categorical_matrix
  }else{
    stop("for calculating ELSA on raster, w must be either a binary matrix or an integer")
  }

  isRaster <- class(rast)[[1]] == "RasterLayer"

  if(isRaster){
    mat <- raster::as.matrix(rast)
  }else{
    mat <- rast
  }

  mat <- ifelse(is.na(mat),-1,mat)
  # let us check if the lowest value beside -1 is 0
  vec <- c(mat)
  m <- min(vec[vec!=-1])
  if(m > 0){
    mat<- ifelse(mat > 0, mat-1,mat)
  }
  cats <- unique(c(mat))
  refcats <- 0:max(cats)
  cats <- cats[cats != -1]
  cats <- cats[order(cats)]
  if(refcats != cats){
    stop(paste("the values of the raster used for ELSA calculation must be integers starting from 0 (or 1).
 There must be no jumps between categories. The categories in the actual raster are : ", paste(cats,collapse = ","),sep=""))
  }
  mat2 <- fun(mat, w, dist)
  mat2 <- ifelse(mat2 == -1, NA, mat2)
  if(isRaster){
    raster::values(rast) <- mat2
    return(rast)
  }else{
    return(mat2)
  }

}

#
# categories <- rep(1, times = 9)
# categories[[5]] <- 2
# dist <- rbind(c(0,1), c(1,0))
# neighmat <- matrix(0, ncol = 9, nrow = 9)
# p1 <- c(1,1,1,2,2,2,2,2,3,3,3,4,4,4,4,4,5,5,5,5,5,5,5,5,6,6,6,6,6,7,7,7,8,8,8,8,8,9,9,9)
# p2 <- c(2,4,5,1,4,5,6,3,2,5,6,1,2,5,7,8,1,2,3,4,6,7,8,9,2,3,4,8,9,4,5,8,4,5,6,7,9,5,6,8)
# for(i in 1:length(p1)){
#   neighmat[p1[[i]], p2[[i]]] <- 1
# }
# nb <- spdep::mat2listw(neighmat)
#
# elsa_vector(categories, nb, dist)
