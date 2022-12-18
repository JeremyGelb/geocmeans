#' @title Check the robustness of a classification by Bootstrap
#'
#' @description Check that the obtained groups are stable by bootstrap
#'
#' @details Considering that the classification produced by a FCM like algorithm
#'   depends on its initial state, it is important to check if the groups
#'   obtained are stable. This function uses a bootstrap method to do so. During
#'   a selected number of iterations (at least 1000), a sample of size n (with
#'   replacement) is drawn from the original dataset. For each sample, the same
#'   classification algorithm is applied and the results are compared with the
#'   reference results. For each original group, the most similar group is
#'   identified by calculating the Jaccard similarity index between the columns
#'   of the two membership matrices. This index is comprised between 0 (exact
#'   difference) and 1 (perfect similarity) and a value is calculated for each
#'   group at each iteration. One can investigate the values obtained to
#'   determine if the groups are stable. Values under 0.5 are a concern and
#'   indicate that the group is dissolving. Values between 0.6 and 0.75 indicate
#'   a pattern in the data, but a significant uncertainty. Values above 0.8
#'   indicate strong groups. The values of the centres obtained at each
#'   iteration are also returned, it is important to ensure that they approximately
#'   follow a normal distribution (or are at least unimodal).
#'
#' @param object A FCMres object, typically obtained from functions CMeans,
#'   GCMeans, SFCMeans, SGFCMeans
#' @param nsim The number of replications to do for the bootstrap evaluation
#' @param maxiter An integer for the maximum number of iterations
#' @param tol The tolerance criterion used in the evaluateMatrices function for
#'   convergence assessment
#' @param init A string indicating how the initial centres must be selected.
#'   "random" indicates that random observations are used as centres "kpp" use
#'   a distance-based method resulting in more dispersed centres at the
#'   beginning. Both of them are heuristic.
#' @param verbose A boolean to specify if the progress bar should be displayed.
#' @param seed An integer used for random number generation. It ensures that the
#' starting centres will be the same if the same value is selected.
#' @return A list of two values: group_consistency: a dataframe indicating
#'   the consistency across simulations each cluster ; group_centres: a list with
#'   a dataframe for each cluster. The values in the dataframes are the centres
#'   of the clusters at each simulation.
#' @export
#' @examples
#' \dontrun{
#' data(LyonIris)
#'
#' #selecting the columns for the analysis
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14",
#'                    "Pct_65","Pct_Img","TxChom1564","Pct_brevet","NivVieMed")
#'
#' #rescaling the columns
#' Data <- sf::st_drop_geometry(LyonIris[AnalysisFields])
#' for (Col in names(Data)){
#'   Data[[Col]] <- as.numeric(scale(Data[[Col]]))
#' }
#'
#' Cmean <- CMeans(Data,4,1.5,500,standardize = FALSE, seed = 456,
#'     tol = 0.00001, verbose = FALSE)
#'
#' validation <- boot_group_validation(Cmean, nsim = 1000, maxiter = 1000,
#'     tol = 0.01, init = "random")
#' }
boot_group_validation <- function(object, nsim = 1000, maxiter = 1000, tol = 0.01, init = "random", verbose = TRUE, seed = NULL){

  if(object$algo %in% c("FCM", "GFCM", "SFCM", "SGFCM") == FALSE){
    stop('bootstrap group validation can only be performed for FCMres object
         if algo is one of "FCM", "GFCM", "SFCM", "SGFCM"')
  }

  ## calulating the lagged dataset if necessary -----------------------
  if(object$algo %in% c("SGFCM", "SFCM")){
    if(object$isRaster == FALSE){
      wdata <- calcLaggedData(object$Data, object$listw, object$lag_method)
    }else{
      dataset <- lapply(1:ncol(object$Data), function(i){
        rast <- object$rasters[[1]]
        vec <- rep(NA, times = terra::ncell(rast))
        vec[object$missing] <- object$Data[,i]
        terra::values(rast, mat = FALSE) <- vec
        return(rast)
      })
      wdata <- calcWdataRaster(object$window, dataset, object$lag_method, object$missing)
    }

  }else{
    wdata <- NULL
  }

  if(is.null(seed) == FALSE){
    set.seed(seed)
  }

  k <- ncol(object$Belongings)

  ## Starting the iteration of the boostraping -----------------------
  if(verbose){
    pb <- txtProgressBar(1, nsim, style = 3)
    print("Starting the bootstrap iterations...")
  }

  all_perm <- lapply(1:nsim, function(i){
    su_res <- boot_worker(object, wdata, tol, maxiter, init)
    if (verbose){
      setTxtProgressBar(pb, i)
    }

    return(su_res)
  })

  ## calculating the consistency of clusters -----------------------

  # the format of this table is :
  # rows are the permutated clusters
  # columns are the original clusters
  if(verbose){
    print("Calculating the Jaccard values...")
  }

  cons_values <- lapply(all_perm, function(results){
    matX <- results$Belongings
    matY <- object$Belongings[results$idx,]
    qual_mat <- calc_jaccard_mat(matY,matX)
    colnames(qual_mat) <- 1:ncol(qual_mat)
    rownames(qual_mat) <- 1:nrow(qual_mat)

    # we must pair the clusters of the original partition
    # and the permutated ones, starting by the best pair
    # and then decreasing. To do so, I edit the matrix
    # and remove the rows / columns used.
    clst_consist <- rep(-1, times = k)
    matidx <- rep(-1, times = k)
    for (j in 1:(k-1)){
      d <- j-1
      best <- which(qual_mat == max(qual_mat), arr.ind = TRUE)
      if(length(best) > 2){
        best <- best[1,]
      }
      c1 <- as.numeric(colnames(qual_mat)[[best[[2]]]])
      r1 <- as.numeric(rownames(qual_mat)[[best[[1]]]])

      clst_consist[[c1]] <- qual_mat[best[[1]], best[[2]]]
      matidx[[c1]] <- r1
      rkeep <- 1:(k-j+1)
      rkeep <- rkeep != best[[1]]
      ckeep <- 1:(k-j+1)
      ckeep <- ckeep != best[[2]]
      qual_mat <- qual_mat[rkeep,ckeep]
    }

    missing <- (1:k)[(1:k) %in% matidx == F]
    clst_consist[clst_consist == -1] <- qual_mat
    matidx[matidx == -1] <- missing

    return(list(clst_consist,matidx))

  })
  mat_valid <- do.call(rbind, lapply(cons_values, function(v){v[[1]]}))
  mat_idx <- do.call(rbind, lapply(cons_values, function(v){v[[2]]}))

  ## creating a list with the centres values boostraped of clusters -----------------------
  print("Extracting the centres of the clusters...")

  clust_table <- do.call(rbind,lapply(1:nsim, function(j){
    all_perm[[j]]$Centers[mat_idx[j,],]
  }))

  clust_table <- as.data.frame(clust_table)
  colnames(clust_table) <- colnames(object$Data)

  gp <- rep(1:k, times = nsim)
  clust_centers <- split(clust_table,f = as.factor(gp))
  names(clust_centers) <- paste("group",1:k, sep = "")
  mat_valid <- data.frame(mat_valid)
  names(mat_valid) <- paste("group",1:k, sep = "")
  return(
    list("group_consistency" = mat_valid,
         "group_centers" = clust_centers
          )
  )
}



#' @title Check that the obtained groups are stable by bootstrap (multicore)
#'
#' @description Check that the obtained groups are stable by bootstrap with
#'   multicore support
#'
#' @details For more details, see the documentation of the function
#'   boot_group_validation
#'
#' @param object A FCMres object, typically obtained from functions CMeans,
#'   GCMeans, SFCMeans, SGFCMeans
#' @param nsim The number of replications to do for the bootstrap evaluation
#' @param maxiter An integer for the maximum number of iterations
#' @param tol The tolerance criterion used in the evaluateMatrices function for
#'   convergence assessment
#' @param init A string indicating how the initial centres must be selected.
#'   "random" indicates that random observations are used as centres. "kpp" use
#'   a distance based method resulting in more dispersed centres at the
#'   beginning. Both of them are heuristic.
#' @param verbose A boolean to specify if the progress bar should be displayed.
#' @param seed An integer to control randomness, default is NULL
#' @return A list of two values: group_consistency: a dataframe indicating
#'   the consistency across simulations each cluster ; group_centres: a list with
#'   a dataframe for each cluster. The values in the dataframes are the centres
#'   of the clusters at each simulation.
#' @export
#' @examples
#' \dontrun{
#' data(LyonIris)
#'
#' #selecting the columns for the analysis
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14",
#'                    "Pct_65","Pct_Img","TxChom1564","Pct_brevet","NivVieMed")
#'
#' #rescaling the columns
#' Data <- sf::st_drop_geometry(LyonIris[AnalysisFields])
#' for (Col in names(Data)){
#'   Data[[Col]] <- as.numeric(scale(Data[[Col]]))
#' }
#'
#' Cmean <- CMeans(Data,4,1.5,500,standardize = FALSE, seed = 456,
#'     tol = 0.00001, verbose = FALSE)
#'
#' future::plan(future::multisession(workers=2))
#'
#' validation <- boot_group_validation.mc(Cmean, nsim = 1000, maxiter = 1000,
#'     tol = 0.01, init = "random")
#' ## make sure any open connections are closed afterward
#' if (!inherits(future::plan(), "sequential")) future::plan(future::sequential)
#' }
boot_group_validation.mc <- function(object, nsim = 1000, maxiter = 1000, tol = 0.01, init = "random", verbose = TRUE, seed = NULL){

  if(object$algo %in% c("FCM", "GFCM", "SFCM", "SGFCM") == FALSE){
    stop('bootstrap group validation can only be performed for FCMres object
         if algo is one of "FCM", "GFCM", "SFCM", "SGFCM"')
  }

  ## calulating the lagged dataset if necessary -----------------------
  if(object$algo %in% c("SGFCM", "SFCM")){
    if(object$isRaster == FALSE){
      wdata <- calcLaggedData(object$Data, object$listw, object$lag_method)
    }else{
      dataset <- lapply(1:ncol(object$Data), function(i){
        rast <- object$rasters[[1]]
        vec <- rep(NA, times = terra::ncell(rast))
        vec[object$missing] <- object$Data[,i]
        terra::values(rast, mat = FALSE) <- vec
        return(rast)
      })
      wdata <- calcWdataRaster(object$window, dataset, object$lag_method, object$missing)
    }

  }else{
    wdata <- NULL
  }

  k <- ncol(object$Belongings)

  ## Starting the iteration of the boostraping -----------------------

  iseq <- 1:nsim
  if(verbose){
    progressr::with_progress({
      p <- progressr::progressor(along = iseq)
      all_perm <- future.apply::future_lapply(iseq, function(i) {

        su_res <- boot_worker(object, wdata, tol, maxiter, init)

        p(sprintf("i=%g", i))

        return(su_res)

      }, future.seed = seed)
    })
  }else{
    all_perm <- future.apply::future_lapply(iseq, function(i) {
      su_res <- boot_worker(object, wdata, tol, maxiter, init)
      return(su_res)
    },future.seed = seed)
  }

  ## calculating the consistency of clusters -----------------------

  # the format of this table is :
  # rows are the permutated clusters
  # columns are the original clusters
  if(verbose){
    print("Calculating the Jaccard values...")
  }

  cons_values <- future.apply::future_lapply(all_perm, function(results){
    matX <- results$Belongings
    matY <- object$Belongings[results$idx,]
    qual_mat <- calc_jaccard_mat(matY,matX)
    colnames(qual_mat) <- 1:ncol(qual_mat)
    rownames(qual_mat) <- 1:nrow(qual_mat)

    # we must pair the clusters of the original partition
    # and the permutated ones, starting by the best pair
    # and then decreasing. To do so, I edit the matrix
    # and remove the rows / columns used.
    clst_consist <- rep(-1, times = k)
    matidx <- rep(-1, times = k)
    for (j in 1:(k-1)){
      d <- j-1
      best <- which(qual_mat == max(qual_mat), arr.ind = TRUE)
      if(length(best) > 2){
        best <- best[1,]
      }
      c1 <- as.numeric(colnames(qual_mat)[[best[[2]]]])
      r1 <- as.numeric(rownames(qual_mat)[[best[[1]]]])

      clst_consist[[c1]] <- qual_mat[best[[1]], best[[2]]]
      matidx[[c1]] <- r1
      rkeep <- 1:(k-j+1)
      rkeep <- rkeep != best[[1]]
      ckeep <- 1:(k-j+1)
      ckeep <- ckeep != best[[2]]
      qual_mat <- qual_mat[rkeep,ckeep]
    }

    missing <- (1:k)[(1:k) %in% matidx == F]
    clst_consist[clst_consist == -1] <- qual_mat
    matidx[matidx == -1] <- missing

    return(list(clst_consist,matidx))

  })

  mat_valid <- do.call(rbind, lapply(cons_values, function(v){v[[1]]}))
  mat_idx <- do.call(rbind, lapply(cons_values, function(v){v[[2]]}))

  ## creating a list with the centers values boostraped of clusters -----------------------
  print("Extracting the centres of the clusters...")

  clust_table <- do.call(rbind,lapply(1:nsim, function(j){
    all_perm[[j]]$Centers[mat_idx[j,],]
  }))

  clust_table <- as.data.frame(clust_table)
  colnames(clust_table) <- colnames(object$Data)
  gp <- rep(1:k, times = nsim)
  clust_centers <- split(clust_table,f = as.factor(gp))
  names(clust_centers) <- paste("group",1:k, sep = "")
  mat_valid <- data.frame(mat_valid)
  names(mat_valid) <- paste("group",1:k, sep = "")

  return(
    list("group_consistency" = mat_valid,
         "group_centers" = clust_centers
    )
  )
}


#' @title Worker function for cluster bootstrapping
#'
#' @description Worker function for cluster bootstrapping
#'
#' @details The worker function for the functions boot_group_validation and boot_group_validation.mc
#'
#' @param object A FCMres object, typically obtained from functions CMeans, GCMeans, SFCMeans, SGFCMeans
#' @param wdata The lagged dataset if necessary, can be NULL if not required
#' @param tol The tolerance criterion used in the evaluateMatrices function for
#'   convergence assessment
#' @param maxiter An integer for the maximum number of iteration
#' @param init A string indicating how the initial centres must be selected. "random"
#' indicates that random observations are used as centres. "kpp" use a distance based method
#' resulting in more dispersed centres at the beginning. Both of them are heuristic.
#' @return A list, similar to a FCMres object, but with only necessary slots for cluster bootstraping.
#' @keywords internal
#' @examples
#' # this is an internal function, no example provided
boot_worker <- function(object, wdata, tol, maxiter, init){

  n <- nrow(object$Data)
  idx <- sample(1:n, size = n, replace = TRUE)

  # step2 : running the classification algorithm
  params <- list(tol = tol,
                 maxiter = maxiter,
                 verbose = FALSE,
                 init = init,
                 wdata = wdata)
  params <- c(object,params)

  params$data <- params$Data[idx,]
  params$wdata <- params$wdata[idx,]

  results <- do.call(main_worker, params)
  su_res <- list(
    Centers = results$Centers,
    Belongings = results$Belongings,
    idx = idx
  )
  return(su_res)

}


#' @title Match the groups obtained from two classifications
#'
#' @description Match the groups obtained from two classifications based on
#' the Jaccard index calculated on the membership matrices.
#'
#' @details We can not expect to obtain the groups in the same order in each run
#' of a classification algorithm. This function can be used match the clusters of a first
#' classification with the most similar clusters in a second classification. Thus
#' it might be easier to compare the results of two algorithms or two runs of the
#' same algorithm.
#'
#' @param object.x A FCMres object, or a simple membership matrix. It is used as the reference
#' for the ordering of the groups
#' @param object.y A FCMres object, or a simple membership matrix. The order of its groups will
#' be updated to match with the groups of object.x
#' @return The FCMres object or the membership matrix provided for the parameter object.y with
#' the order of the groups updated.
#' @export
#' @examples
#' data(LyonIris)
#'
#' #selecting the columns for the analysis
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14",
#'                    "Pct_65","Pct_Img","TxChom1564","Pct_brevet","NivVieMed")
#'
#' #rescaling the columns
#' Data <- sf::st_drop_geometry(LyonIris[AnalysisFields])
#' for (Col in names(Data)){
#'   Data[[Col]] <- as.numeric(scale(Data[[Col]]))
#' }
#'
#' Cmean <- CMeans(Data,4,1.5,500,standardize = FALSE, seed = 456, tol = 0.00001, verbose = FALSE)
#' Cmean2 <- CMeans(Data,4,1.5,500,standardize = FALSE, seed = 789, tol = 0.00001, verbose = FALSE)
#' ordered_Cmean2 <- groups_matching(Cmean,Cmean2)
groups_matching <- function(object.x,object.y){

  if(is.matrix(object.x)){
    matX <- object.x
    matY <- object.y
  }else{
    matX <- object.x$Belongings
    matY <- object.y$Belongings
  }

  k <- ncol(matX)
  qual_mat <- calc_jaccard_mat(matX,matY)
  colnames(qual_mat) <- 1:ncol(qual_mat)
  rownames(qual_mat) <- 1:nrow(qual_mat)

  clst_consist <- rep(-1, times = k)
  matidx <- rep(-1, times = k)

  for (j in 1:(k-1)){
    d <- j-1
    best <- which(qual_mat == max(qual_mat), arr.ind = TRUE)
    if(length(best) > 2){
      best <- best[1,]
    }
    c1 <- as.numeric(colnames(qual_mat)[[best[[2]]]])
    r1 <- as.numeric(rownames(qual_mat)[[best[[1]]]])

    clst_consist[[c1]] <- qual_mat[best[[1]], best[[2]]]
    matidx[[c1]] <- r1
    rkeep <- 1:(k-j+1)
    rkeep <- rkeep != best[[1]]
    ckeep <- 1:(k-j+1)
    ckeep <- ckeep != best[[2]]
    qual_mat <- qual_mat[rkeep,ckeep]
  }

  missing <- (1:k)[(1:k) %in% matidx == F]
  clst_consist[clst_consist == -1] <- qual_mat
  matidx[matidx == -1] <- missing

  if(is.matrix(object.x)){
    return(matY[,matidx])
  }else{
    object.y$Centers <- object.y$Centers[matidx,]
    object.y$Belongings <- object.y$Belongings[,matidx]
    DF <- as.data.frame(object.y$Belongings)
    names(DF) <- paste("group",1:k,sep = "")
    object.y$Groups <- colnames(DF)[max.col(DF, ties.method = "first")]
    if(object.y$isRaster){
      oldnames <- names(object.y$rasters)
      #if we have rasters, we must update their order too
      membership_rasts <- object.y$rasters[1:object.y$k]
      membership_rasts <- membership_rasts[matidx]
      i <- length(object.y$rasters)
      membership_rasts[[i]] <- object.y$rasters[[i]]
      object.y$rasters <- membership_rasts
      names(object.y$rasters) <- oldnames
    }
    return(object.y)
  }
}
