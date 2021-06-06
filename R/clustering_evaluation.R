# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Indices of clustering quality #####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Negentropy Increment index
#'
#' @description Calculate the Negentropy Increment index of clustering quality.
#' This index based on the assumption that a normally shapped cluster is more
#' desirable. It uses the difference between the average negentropy
#' of all the clusters in the partition, and that of the  whole partition.
#' a smaller value indicates a better partition.
#'
#' @param data The original dataframe used for the clustering (n*p)
#' @param belongmatrix A membership matrix (n*k)
#' @param centers The centers of the clusters
#' @return A float : the David-Bouldin index
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
#' This index is a ratio of the worst pair-wise separation of clusters and the
#' worst pair-wise compactness of clusters. Thus a higher value indicates a
#' better clustering.
#'
#' @param data The original dataframe used for the clustering (n*p)
#' @param belongmatrix A membership matrix (n*k)
#' @param centers The centers of the clusters
#' @return A float : the David-Bouldin index
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
#' This index is a ratio of the worst pair-wise separation of clusters and the
#' worst pair-wise compactness of clusters. Thus a higher value indicates a
#' better clustering.
#'
#' @param data The original dataframe used for the clustering (n*p)
#' @param belongmatrix A membership matrix (n*k)
#' @param centers The centers of the clusters
#' @return A float : the David-Bouldin index
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



#' @title David-Bouldin index
#'
#' @description Calculate the David-Bouldin index of clustering quality. This index
#' can be seen as a ratio between the within cluster dispersion and the
#' between cluster separation. A lower value indicates a higher cluster compacity
#' or a higher cluster separation.
#'
#' @param data The original dataframe used for the clustering (n*p)
#' @param belongmatrix A membership matrix (n*k)
#' @param centers The centers of the clusters
#' @return A float : the David-Bouldin index
#' @export
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' calcDavidBouldin(result$Data, result$Belongings, result$Centers)
calcDavidBouldin <- function(data, belongmatrix, centers){
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
#' @description Calculate the Calinski-Harabasz index of clustering quality. This index
#' is the ratio between the clusters separation and the clusters cohesion. A greater
#' value indicates eithermore separated clusters or more cohesive clusters.
#'
#' @param data The original dataframe used for the clustering (n*p)
#' @param belongmatrix A membership matrix (n*k)
#' @param centers The centers of the clusters
#' @return A float : the Calinski-Harabasz index
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
#' @param data The original dataframe used for the clustering (n*p)
#' @param belongmatrix A membership matrix (n*k)
#' @param centers The centers of the clusters
#' @param m The fuzzyness parameter
#' @return A float : the Fukuyama and Sugeno index
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
    t1 <- calcEuclideanDistance(data,v)
    t2 <- sum((v - v_hat)**2)
    #t3 <- calcEuclideanDistance(wdata,v) ajouter l'espace dans l'indice ?
    #dists <- ((t1+alpha*t3)/(1+alpha) - t2)*w
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
#' @return A float : the percentage of the total inertia explained
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
  baseinertia <- sum(calcEuclideanDistance(data, means))
  restinertia <- sapply(1:nrow(centers), function(i) {
    point <- centers[i, ]
    dists <- calcEuclideanDistance(data, point) * belongmatrix[, i]
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
#' @return A float : the value of the index
#' @keywords internal
#' @examples
#' # this is an internal function, no example provided
calcQualIdx <- function(name, ...){
  dots <- list(...)
  val <- NULL
  if(name == "Silhouette.index"){
    val <- fclust::SIL.F(dots$data, dots$belongmatrix, alpha=1)
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
    val <- calcDavidBouldin(dots$data, dots$belongmatrix, dots$centers)
  }

  if(name == "DavidBoulin.index"){
    val <- calcDavidBouldin(dots$data, dots$belongmatrix, dots$centers)
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
    "DavidBoulin.index", "CalinskiHarabasz.index", "GD43.index", "GD53.index",
    "Negentropy.index"')
  }
  return(val)
}



#' @title Quality indexes
#'
#' @description calculate several clustering quality indexes (most of them come from fclust
#' package)
#'
#' @param data The original dataframe used for the classification (n*p)
#' @param belongmatrix A membership matrix (n*k)
#' @param m The fuzziness parameter used for the classification
#' @param indices A character vector with the names of the indices to calculate, default is :
#' c("Silhouette.index", "Partition.entropy", "Partition.coeff", "XieBeni.index", "FukuyamaSugeno.index",
#' "Explained.inertia"). Other available indices are : "DavidBoulin.index", "CalinskiHarabasz.index",
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

  # idxsf <- fclust::SIL.F(data, belongmatrix, alpha=1)  #look for maximum
  # idxpe <- fclust::PE(belongmatrix)  #look for minimum
  # idxpc <- fclust::PC(belongmatrix)  #loook for maximum
  # idxmpc <- fclust::MPC(belongmatrix)  #look for maximum
  # idxfukusu <- calcFukuyamaSugeno(data, belongmatrix, centers, m) # look for minimum
  # idxHC <- calcCalinskiHarabasz(data, belongmatrix, centers) # look for maximum
  # idxDB <- calcDavidBouldin(data, belongmatrix, centers) # look for minimum
  # if(m>1){
  #   idxXB <- fclust::XB(data,belongmatrix,centers,m) # look for minimum
  # }else{
  #   idxXB <- NA
  # }
  #
  # # calculating the explained inertia
  # explainedinertia <- calcexplainedInertia(data,belongmatrix)
#
#   return(list(Silhouette.index = idxsf,
#               Partition.entropy = idxpe,
#               Partition.coeff = idxpc,
#               Modified.partition.coeff = idxmpc,
#               XieBeni.index = idxXB,
#               FukuyamaSugeno.index = idxfukusu,
#               CalinskiHarabasz.index = idxHC,
#               DavidBoulin.index = idxDB,
#               Explained.inertia = explainedinertia))
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Indices of spatial clustering quality #####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Spatial diagnostic
#'
#' @description Utility function to facilitate the spatial diagnostic of a classification
#'
#' Calculate the following indicators : Moran I index (spdep::moranI) for each
#' column of the belonging matrix, Join count test (spdep::joincount.multi) for
#' the most likely groups of each datapoint, Spatial consistency index (see
#' function spConsistency)
#'
#' @param belongmatrix A membership matrix
#' @param nblistw A list.w  object describing the neighbours (spdep package)
#' @param undecided A float between 0 and 1 giving the minimum value that an
#'   observation must get in the membership matrix to not be considered as
#'   uncertain (default = NULL)
#' @param nrep An integer indicating the number of permutation to do to simulate
#'   the random distribution of the spatial inconsistency
#' @return A named list with :
#' \itemize{
#'         \item MoranValues : the moran I values fo each column of the membership
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
#' spatialDiag(result$Belongings, Wqueen, undecided=0.45, nrep=30)
spatialDiag <- function(belongmatrix, nblistw, undecided = NULL, nrep = 50) {
  # calcul de la coherence spatiale
  Consist <- spConsistency(belongmatrix, nblistw = nblistw, nrep = nrep)
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
  return(list(MoranValues = morandf, JoinCounts = spjctetst, SpConsist = Consist$Mean, SpConsistSamples = Consist$samples))
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
#' @param belongmatrix A membership matrix
#' @param nblistw A list.w object describing the neighbours (spdep package)
#'   observation must get in the membership matrix to not be considered as
#'   uncertain (default = NULL)
#' @param nrep An integer indicating the number of permutation to do to simulate
#' @return A named list with
#'  \itemize{
#'         \item Mean : the mean of the spatial consistency index
#'         \item prt05 : the 5th percentile of the spatial consistency index
#'         \item prt95 : the 95th percentule of the spatial consistency index
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
#' spConsistency(result$Belongings, Wqueen, nrep=50)
spConsistency <- function(belongmatrix, nblistw, nrep = 999) {
  belongmat <- as.matrix(belongmatrix)
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
  return(list(Mean = mean(ratio), Median = quantile(ratio, probs = c(0.5)),
              prt05 = quantile(ratio, probs = c(0.05)),
              prt95 = quantile(ratio, probs = c(0.95)),
              samples = ratio))
}
