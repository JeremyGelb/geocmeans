#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### spatial consistency index ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
#' @template FCMresobj-arg
#' @template nblistw2-arg
#' @param window if rasters were used for the classification, the window must be
#' specified instead of a list.w object. Can also be NULL if object is a FCMres object.
#' @param nrep An integer indicating the number of permutation to do to simulate
#' spatial randomness. Note that if rasters are used, each permutation can be very long.
#' @param adj A boolean indicating if the adjusted version of the indicator must be
#' calculated when working with rasters. When working with vectors, see the function
#' adjustSpatialWeights to modify the list.w object.
#' @return A named list with
#'  \itemize{
#'         \item Mean : the mean of the spatial consistency index
#'         \item prt05 : the 5th percentile of the spatial consistency index
#'         \item prt95 : the 95th percentile of the spatial consistency index
#'         \item samples : all the value of the spatial consistency index
#'         \item sum_diff : the total sum of squarred difference between observations and their neighbours
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
spConsistency <- function(object, nblistw = NULL, window = NULL, nrep = 999, adj = FALSE) {

  if(inherits(object, "FCMres")){
    belongmat <- as.matrix(object$Belongings)
    if(object$isRaster & is.null(window)){
      window <- object$window
      if(is.null(window)){
        stop("impossible to find a window in the given object, please
             specify one by hand.")
      }
    }
    if(object$isRaster == FALSE & is.null(nblistw)){
      nblistw <- object$nblistw
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
      diff <- (neighbour-row[col(neighbour)])**2 * W
      tot <- sum(rowSums(diff))
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
    mat_dim <- dim(matrices[[1]])

    if(adj){
      dataset <- lapply(1:ncol(object$Data), function(ic){
        vec1 <- object$Data[,ic]
        vec2 <- rep(NA,length(object$missing))
        vec2[object$missing] <- vec1
        rast <- object$rasters[[1]]
        raster::values(rast) <- vec2
        mat <- as.matrix(rast)
        return(mat)

      })
      totalcons <- calc_raster_spinconsistency(matrices, window, adj, dataset)

    }else{
      totalcons <- calc_raster_spinconsistency(matrices,window)
    }


    # we must now do the same thing but with resampled values
    warning("Calculating the permutation for the spatial inconsistency
            when using raster can be long, depending on the raster size.
            Note that the high number of cell in a raster reduces the need of
            a great number of replications.")
    # creating a vector of ids for each cell in raster
    all_ids <- 1:raster::ncell(rasters[[1]])

    # converting the matrices (columns of membership matrix) into 1d vectors
    mem_vecs <- lapply(rasters, function(rast){
      mat <- raster::as.matrix(rast)
      dim(mat) <- NULL
      return(mat)
    })

    # if necessary, doing the same with the original data
    if(adj){
      data_vecs <- lapply(dataset, function(mat){
        vec <- mat
        dim(vec) <- NULL
        return(vec)
      })
    }
    # extracting the dimension of the raster
    #dim(raster::as.matrix(rasters[[1]]))

    # starting the simulations
    simulated <- sapply(1:nrep, function(i){

      # resampling the ids
      Ids <- sample(all_ids)

      # resampling the matrices of memberships
      new_matrices <- lapply(mem_vecs, function(vec){
        new_vec <- vec[Ids]

        # swapping the NA at their original place
        val_na <- new_vec[!object$missing]
        loc_na <- is.na(new_vec)
        new_vec[!object$missing] <- NA
        new_vec[loc_na] <- val_na

        dim(new_vec) <- mat_dim
        return(new_vec)
      })

      if(adj){
        # resampling the matrices of the data
        new_dataset <- lapply(data_vecs, function(vec){
          new_vec <- vec[Ids];

          # swapping the NA at their original place
          val_na <- new_vec[!object$missing]
          loc_na <- is.na(new_vec)
          new_vec[!object$missing] <- NA
          new_vec[loc_na] <- val_na

          dim(new_vec) <- mat_dim
          return(new_vec)
        })
        # calculating the index value
        inconsist <-  calc_raster_spinconsistency(new_matrices, window, adj, new_dataset)

      }else{
        # calculating the index value
        inconsist <-  calc_raster_spinconsistency(new_matrices, window)
      }


      return(inconsist)
    })
    ratio <- totalcons / simulated
  }

  return(list(Mean = mean(ratio), Median = quantile(ratio, probs = c(0.5)),
              prt05 = quantile(ratio, probs = c(0.05)),
              prt95 = quantile(ratio, probs = c(0.95)),
              samples = ratio,
              sum_diff = totalcons))
}



#' @title calculate spatial inconsistency for raster
#'
#' @description Calculate the spatial inconsistency sum for a set of rasters
#'
#' @param matrices A list of matrices
#' @param window The window to use to define spatial neighbouring
#' @param adj A boolean indicating if the adjusted version of the algorithm must be
#' calculated
#' @param dataset A list of matrices with the original data (if adj = TRUE)
#' @return A float: the sum of spatial inconsistency
#' @keywords internal
#' @examples
#' # this is an internal function, no example provided
calc_raster_spinconsistency <- function(matrices, window, adj = FALSE, dataset = NULL){
  if(adj & is.null(dataset)){
    stop("When the adjusted version of spinconsistency is required, dataset must be given")
  }
  arr <- array(do.call(c,matrices), c(nrow(matrices[[1]]), ncol(matrices[[1]]), length(matrices)))
  if(adj){
    arr2 <- array(do.call(c,dataset), c(nrow(dataset[[1]]), ncol(dataset[[1]]), length(dataset)))
    totalcons <- adj_spconsist_arr_window_globstd(arr2, arr, window)
  }else{
    totalcons <- sum(focal_euclidean_arr_window(arr,window), na.rm = TRUE)
  }
  return(totalcons)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### ELSA for HARD PARTITION ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Check validity of a dissimilarity matrix
#'
#' @description Check the validity of a dissimilarity matrix
#'
#' @param matdist A dissimilarity matrix
#' @keywords internal
#' @examples
#' # this is an internal function, no example provided
check_matdist <- function(matdist){
  if(!isSymmetric(matdist)){
    stop("matdist must be a symetric matrix")
  }
  if(sum(diag(matdist)) != 0){
    stop("matdist must have a diagonal filled with zeros")
  }
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
#' autocorelation. The authors suggest to calculate the mean of the index at several lag
#' distance to create an entrogram which quantifies global spatial structure and can be
#' represented as a variogram-like graph.
#'
#' @param object A FCMres object, typically obtained from functions CMeans,
#'   GCMeans, SFCMeans, SGFCMeans. Can also be a vector of categories. This vector must
#'   be filled with integers starting from 1. -1 can be used to indicate missing categories.
#' @template nblistw-arg
#' @template window2-arg
#' @template matdist-arg
#' @return A depending of the input, a vector of ELSA values or a raster with the ELSA values.
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' elsa_valus <- calcELSA(result)
calcELSA <- function(object, nblistw = NULL, window = NULL, matdist = NULL){

  # testing if we have all the required parameters
  if(inherits(object, "FCMres") == FALSE){

    if(inherits(object, "numeric") == FALSE){
      stop("if object is not a FCMres object, it must be a numeric vector")
    }

    vec <- object[object != -1]
    if(min(vec) != 1){
      stop("if object is a numeric vector, its lower value must be 1")
    }

    if(length(unique(vec)) != length(1:max(vec))){
      stop("if object is a numeric vector, its values must be like 1,2,3,4,...m, -1 can be used to indicate missing values")
    }

    if(is.null(matdist)){
      stop("if object is not a FCMres object, matdist must be provided")
    }else{
      check_matdist(matdist)
    }


    if(nrow(matdist) != length(unique(vec))){
      stop("the dimension of matdist must equal the number of categories in object.")
    }

    if(is.null(nblistw)){
      stop("if object is not a FCMres object, nblistw must be provided.")
    }
  }else{
    if(object$isRaster){
      if(is.null(object$window) & is.null(window)){
        stop("impossible to extract window from object, window must be given.")
      }
    }else{
      if(is.null(object$nblistw) & is.null(nblistw)){
        stop("impossible to extract nblistw from object, nblistw must be given.")
      }
    }
  }

  # case 1 : object is a simple vector of categories
  if(inherits(object,'FCMres') == FALSE){
    vals <- elsa_vector(object, nblistw, matdist)
  }else{
    # case 2 : the object is a FCMres object
    if(is.null(matdist)){
      matdist <- as.matrix(stats::dist(object$Centers))
    }

    if(object$isRaster){
      if(is.null(window)){
        window <- object$window
      }
      vals <- elsa_raster(object$rasters$Groups, window, matdist)
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
      probs <- probs[probs > 0]
      if(mi == 1){
        return(0)
      }else{
        Eci <- -1 * (sum(probs * log2(probs)) / log2(mi))
        return(Eai * Eci)
      }
    }

  })
  return(values)

}


#' @title calculate ELSA spatial statistic for raster dataset
#'
#' @description calculate ELSA spatial statistic for vector dataset
#'
#' @param rast An integer raster or matrix representing the m categories (0,1,2,..., m)
#' @template window2-arg
#' @template matdist-arg
#' @return A raster or a matrix: the local values of ELSA
#' @keywords internal
#' @examples
#' # this is an internal function, no example provided
elsa_raster <- function(rast, window, matdist){

  if(inherits(window, "matrix")){
    u <- unique(window)
    if(length(union(c(0,1),u)) > 2){
      stop("The provided matrix to calculate ELSA is not a binary matrix (0,1)")
    }
    fun <- Elsa_categorical_matrix_window
  }else{
    stop("for calculating ELSA on raster, window must be a binary matrix")
  }

  isRaster <- inherits(rast, "RasterLayer")

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
  if(sum(refcats - cats)!=0){
    stop(paste("the values of the raster used for ELSA calculation must be integers starting from 0 (or 1).
 There must be no jumps between categories. The categories in the actual raster are : ", paste(cats,collapse = ","),sep=""))
  }
  mat2 <- fun(mat, window, matdist)
  mat2 <- ifelse(mat2 == -1, NA, mat2)
  if(isRaster){
    raster::values(rast) <- mat2
    return(rast)
  }else{
    return(mat2)
  }

}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### ELSA for fuzzy PARTITION ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title calculate ELSA statistic for a fuzzy partition
#'
#' @description Calculate ELSA statistic for a fuzzy partition. This local indicator of
#' spatial autocorrelation can be used to identify areas where close observations tend to
#' belong to different clusters.
#'
#' @details The fuzzy ELSA index is a generalization of the ELSA index \insertCite{naimi2019elsa}{geocmeans}. It can be used to measure
#' local autocorrelation for a membership matrix. It varies between 0 and 1, 0 indicating
#' a perfect positive spatial autocorrelation and 1 a perfect heterogeneity. It is based on
#' the Shannon entropy index, and uses a measure of dissimilarity between categories.
#'
#' @param object A FCMres object, typically obtained from functions CMeans,
#'   GCMeans, SFCMeans, SGFCMeans. Can also be a membership matrix. Each row of this matrix
#'   must sum up to 1. Can also be a list of rasters, in which case each raster must represent
#'   the membership values for one cluster and the sum of all the rasters must be a raster filled
#'   with ones.
#' @template nblistw-arg
#' @template window2-arg
#' @template matdist-arg
#' @return either a vector or a raster with the ELSA values.
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' elsa_valus <- calcFuzzyELSA(result)
calcFuzzyELSA <- function(object, nblistw = NULL, window = NULL, matdist = NULL){

  # testing if we have all the required parameters
  if(inherits(object, c("FCMres","matrix","list")) == FALSE){
    stop("object must be a FCMres object, a matrix or a list of rasters")
  }

  if(inherits(object,"matrix")){
    sums <- rowSums(object) != 1
    if(any(sums)){
      stop("if object is a matrix, the sum of each row must be 1")
    }
  }

  if(inherits(object,"FCMres") == FALSE & is.null(matdist)){
    stop("if object is not a FCMre, matdist must be specified")
  }

  if(inherits(object,"matrix") & is.null(nblistw)){
    stop("if object is a matrix, nblistw must be specified")
  }

  if(inherits(object,"list") & inherits(object,"list") == FALSE & is.null(window)){
    stop("if object is a list, window must be specified")
  }

  if(is.null(matdist)==FALSE){
    if(!isSymmetric(matdist)){
      stop("matdist must be a symetric matrix")
    }
  }


  if(inherits(object,"matrix")){
    if(nrow(matdist) != ncol(object)){
      stop("the number of columns in object (matrix) must match the number of rows in matdist")
    }
    if(is.null(nblistw)){
      stop("if object is a matrix, nblistw must be provided")
    }
  }

  if(inherits(object,"list") & (inherits(object,"FCMres") == FALSE)){
    if(nrow(matdist) != length(object)){
      stop("the number of rasters in object (list) must match the number of rows in matdist")
    }
    if(is.null(window)){
      stop("if object is a list of rasters, window must be provided")
    }
  }

  if(inherits(object,"FCMres")){
    if(object$isRaster){
      if(is.null(object$window) & is.null(window)){
        stop("impossible to extract window from object, window must be provided")
      }
    }else{
      if(is.null(object$nblistw) & is.null(nblistw)){
        stop("impossible to extract nblistw from object, nblistw must be provided")
      }
    }
  }


  if(inherits(object,"FCMres")){
    if(object$isRaster==FALSE){
      #case 1 : object is a FCMres and vector type
      mat <- object$Belongings

      if(is.null(matdist)){
        matdist <- as.matrix(stats::dist(object$Centers))
      }

      if(is.null(nblistw)){
        nblistw <- object$nblistw
      }

      return(elsa_fuzzy_vector(mat, nblistw, matdist))

    }else{
      #case 2 : object is a FCMres and raster type
      if(is.null(matdist)){
        matdist <- as.matrix(stats::dist(object$Centers))
      }
      if(is.null(window)){
        window <- object$window
      }
      mats <- lapply(object$rasters[1:object$k], raster::as.matrix)
      arr <- do.call(c,mats)
      arr <- array(arr, dim = c(nrow(mats[[1]]), ncol(mats[[1]]),object$k))
      values <- Elsa_fuzzy_matrix_window(arr, window, matdist)
      rast <- object$rasters[[1]]
      raster::values(rast) <- values
      return(rast)
    }

  }else if (inherits(object,"matrix")){
    #case 3 : object is a matrix
    values <- elsa_fuzzy_vector(object, nblistw, matdist)
    return(values)
  }else if (inherits(object,"list")){
    #case 4 : object is a list of rasters
    mats <- lapply(object, raster::as.matrix)
    arr <- do.call(c,mats)
    arr <- array(arr, dim = c(nrow(mats[[1]]), ncol(mats[[1]]), length(object)))
    values <- Elsa_fuzzy_matrix_window(arr, window, matdist)
    rast <- object[[1]]
    raster::values(rast) <- values
    return(rast)
  }else{
    stop("invalid arguments passed to calcFuzzyELSA")
  }
  return(values)
}



#' @title Local Fuzzy ELSA statistic for vector
#'
#' @description Calculate the Local Fuzzy ELSA statistic using a nblistw object
#'
#' @param memberships A membership matrix
#' @param nbslitw The spatial weight matrix (nblistw object from spdep)
#' @template matdist-arg
#' @return A vector of local ELSA values
#' @keywords internal
#' @examples
#' # this is an internal function, no example provided
elsa_fuzzy_vector <- function(memberships, nblistw, matdist){
  d <- max(matdist)
  m <- ncol(memberships)
  values <- sapply(1:nrow(memberships), function(i){
    xi <- memberships[i,]
    neighbours <- nblistw$neighbours[[i]]
    xj <- memberships[neighbours,]
    diffs <- abs(xi - t(xj))
    sis <- apply(diffs, 2, function(x){
      out <- outer(x,x)
      return(sum(out * matdist))
    })
    s1 <- sum(sis)/2.0
    nn <- length(neighbours)
    eai <- s1 / (d*nn)
    pks <- (colSums2(xj) + xi) / (nn+1)
    pks <- pks[pks>0]
    if(nn > m){
      mi <- m
    }else{
      mi <- nn
    }
    eci <- -1 * (sum(pks * log2(pks)) / log2(mi))
    return(eai * eci)
  })
  return(values)
}



#' @title Local Fuzzy ELSA statistic for raster
#'
#' @description Calculate the Local Fuzzy ELSA statistic for a numeric raster
#'
#' @param rasters A List of RasterLayer or a List of matrices, or an array
#' @template window2-arg
#' @template matdist-arg
#' @return A raster or a matrix (depending on the input): the values of
#' local fuzzy ELSA statistic
#' @keywords internal
#' @examples
#' # this is an internal function, no example provided
calcFuzzyElsa_raster <- function(rasters,window,matdist){
  if(inherits(rasters[[1]], "matrix")){
    mats <- rasters
    isRaster <- FALSE
    vals <- do.call(c,mats)
    arr <- array(vals, c(nrow(mats[[1]]),ncol(mats[[1]]),length(rasters)))
  }else if (inherits(rasters[[1]], "RasterLayer")){
    mats <- lapply(rasters, raster::as.matrix)
    isRaster <- TRUE
    vals <- do.call(c,mats)
    arr <- array(vals, c(nrow(mats[[1]]),ncol(mats[[1]]),length(rasters)))
  }else if (inherits(rasters,"array")){
    arr <- rasters
    isRaster <- FALSE
  }else{
    stop("rasters must be a list of matrix or a list of RasterLayer or an array")
  }
  elsa <- Elsa_fuzzy_matrix_window(arr, window, matdist)
  if(isRaster){
    rast <- rasters[[1]]
    raster::values(rast) <- elsa
  }else{
    dims <- dim(arr)
    rast <- matrix(elsa, nrow = dims[[1]], ncol = dims[[2]])
  }
  return(rast)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Moran for rasters ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Global Moran I for raster
#'
#' @description Calculate the global Moran I for a numeric raster
#'
#' @param rast A RasterLayer or a matrix
#' @param window The window defining the neighbour weights
#' @return A float: the global Moran I
#' @examples
#' data("Arcachon")
#' rast <- Arcachon[[1]]
#' w <- matrix(1, nrow = 3, ncol = 3)
#' calc_moran_raster(rast, w)
calc_moran_raster <- function(rast, window){
  if(inherits(rast,"RasterLayer")){
    mat <- raster::as.matrix(rast)
  }else if (inherits(rast,"matrix")){
    mat <- rast
  }else{
    stop("rast parameter must be on of matrix or RasterLayer")
  }
  if(inherits(window,"matrix")){
    fun <- moranI_matrix_window
  }else if (inherits(window, "numeric")){
    #fun <- moranI_matrix
    fun <- moranI_matrix_window
    size <- n+1+n
    window <- matrix(1, nrow = size, ncol = size)
  }else{
    stop("window parameter must be an integer or a matrix")
  }
  val <- fun(mat, window)
  return(val)
}


#' @title Local Moran I for raster
#'
#' @description Calculate the Local Moran I for a numeric raster
#'
#' @param rast A RasterLayer or a matrix
#' @param window The window defining the neighbour weights
#' @return A RasterLayer or a matrix depending on the input with the local Moran I values
#' @examples
#' data("Arcachon")
#' rast <- Arcachon[[1]]
#' w <- matrix(1, nrow = 3, ncol = 3)
#' calc_local_moran_raster(rast, w)
calc_local_moran_raster <- function(rast, window){
  if(inherits(rast,"RasterLayer")){
    mat <- raster::as.matrix(rast)
  }else if (inherits(rast,"matrix")){
    mat <- rast
  }else{
    stop("rast parameter must be on of matrix or RasterLayer")
  }
  if(inherits(window,"matrix")){
    window <- window
  }else if (inherits(window,"numeric")){
    w <- 1+2*inherits
    window <- matrix(1, nrow = w, ncol = w)
    window[ceiling(w/2),ceiling(w/2)] <- 0
  }else{
    stop("window parameter must be an integer or a matrix")
  }
  vals <- local_moranI_matrix_window(mat, window)

  if(inherits(rast,"RasterLayer")){
    raster::values(rast) <- vals
    return(rast)
  }else{
    return(matrix(vals, ncol = ncol(mat), nrow = nrow(mat)))
  }

}
