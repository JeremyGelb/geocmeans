#' @title Boostrap check the robustness of a classification
#'
#' @description Check that the obtained groups are stable by boostrap
#'
#' @details Considering that the classification produced by a FCM like algorithm depends on
#' its initial state, it is important to check if the groups obtained are stable. This function
#' use a bootstrap method to do so. During a selected number of iterations (at least 1000), a sample
#' of size n (with replacement) is drawn from the original dataset. For each sample, the same classification
#' algorithm is applied and the results are compared with the reference results. For each original group, the
#' mots similar group is identified by calculating the Jaccard similarity index between the columns
#' of the two membership matrices. This index is comprised between 0 (exact difference) and 1 (perfect
#' similarity) and a value is calculated for each group at each iteration. One can investigate the values
#' obtained to determine if the groups are stable. Values under 0.5 are worrisome and indicate that the
#' group is disolving. Values between 0.6 and 0.75 indicate a pattern in the data, but an important
#' uncertainty. Values above 0.8 indicate strong groups. The values of the centers obtained at
#' each iteration are also returned, it is important to check the they follow approximately a normal
#' distribution and are at least unimodal.
#'
#' @param object A FCMres object, typically obtained from functions CMeans, GCMeans, SFCMeans, SGFCMeans
#' @param nsim The number of replications to do for the boostrap evaluation
#' @param maxiter An integer for the maximum number of iteration
#' @param tol The tolerance criterion used in the evaluateMatrices function for
#'   convergence assessment
#' @param init A string indicating how the initial centers must be selected. "random"
#' indicates that random observations are used as centers. "kpp" use a distance based method
#' resulting in more dispersed centers at the beginning. Both of them are heuristic.
#' @return  A list of two values: group_consistency, a dataframe indicating for each cluster its
#' consistency accross simulations. group_centers a list with a dataframe for each cluster. The values
#' in the dataframes are the centers of the clusters at each simulation.
#' @examples
#' \dontrun{
#' data(LyonIris)
#'
#' #selecting the columns for the analysis
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14",
#'                    "Pct_65","Pct_Img","TxChom1564","Pct_brevet","NivVieMed")
#'
#' #rescaling the columns
#' Data <- LyonIris@data[AnalysisFields]
#' for (Col in names(Data)){
#'   Data[[Col]] <- as.numeric(scale(Data[[Col]]))
#' }
#'
#' Cmean <- CMeans(Data,4,1.5,500,standardize = FALSE, seed = 456, tol = 0.00001, verbose = FALSE)
#'
#' validation <- bstp_group_validation(Cmean, nsim = 1000, maxiter = 1000, tol = 0.01, init = "random")
#' }
bstp_group_validation <- function(object, nsim = 1000, maxiter = 1000, tol = 0.01, init = "random"){

  ## calulating the lagged dataset if necessary -----------------------
  if(object$algo %in% c("SGFCM", "SFCM")){
    wdata <- calcLaggedData(object$Data, object$listw, object$lag_method)
  }else{
    wdata <- NULL
  }

  k <- ncol(object$Belongings)

  ## Starting the iteration of the boostraping -----------------------

  all_perm <- lapply(1:nsim, function(i){

    # step1 : selecting the ids of the samples to keep in this iteration
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
  })

  ## calculating the consistency of clusters -----------------------

  # the format of this table is :
  # rows are the permutated clusters
  # columns are the original clusters
  cons_values <- lapply(all_perm, function(results){
    qual_mat <- matrix(-1, nrow = k, ncol = k)
    for(ii in 1:k){
      for(jj in 1:k){
        if(qual_mat[ii,jj] == -1){
          x <- results$Belongings[,ii]
          y <- object$Belongings[,jj][results$idx]
          val <- calc_jaccard_idx(x,y)
          qual_mat[ii,jj] <- val
        }
      }
    }
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

  ## creatin a list with the centers values boostraped of clusters -----------------------
  clust_centers <- lapply(1:k, function(i){
    idx <- mat_idx[,i]
    clust_table <- do.call(rbind,lapply(1:nsim, function(j){
      all_perm[[j]]$Centers[idx[[j]],]
    }))
    clust_table <- data.frame(clust_table)
    names(clust_table) <- colnames(object$Data)
  })
  names(clust_centers) <- paste("group ",1:k, sep = "")
  mat_valid <- data.frame(mat_valid)
  names(mat_valid) <- paste("group ",1:k, sep = "")


  return(
    list("group_consistency" = mat_valid,
         "group_centers" = clust_centers
          )
  )
}





#' @title Jaccard similarity coefficient
#'
#' @description Calculate the Jaccard similarity coefficient
#'
#' @details Considering that the classification produced by
#'
#' @param x A vector of positive reals
#' @param y A vector of positive reals
#' @return A float: the Jaccard similarity coefficient
#' @keywords internal
#' @examples
#' # this is an internal function, no example provided
calc_jaccard_idx <- function(x,y){
  mat <- cbind(x,y)
  num <- sum(apply(mat,1,min))
  denom <- sum(apply(mat,1,max))
  return(num/denom)
}
