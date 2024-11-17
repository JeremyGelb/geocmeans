# # data(LyonIris)
# # AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
# #                    "TxChom1564","Pct_brevet","NivVieMed")
# # dataset <- sf::st_drop_geometry(LyonIris[AnalysisFields])
# # queen <- spdep::poly2nb(LyonIris,queen=TRUE)
# # Wqueen <- spdep::nb2listw(queen,style="W")
# # result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
# #
# # library(compositions)
# #
# # x <- result$Belongings
# # listw <- Wqueen
# #
# #
# #
# # cmp_geary <- function(x, listw){
# #
# #   cmp <- acomp(x)
# #   x_clr <- as.matrix(clr(cmp))
# #
# #   sums_ij <- sapply(1:nrow(x_clr), function(i){
# #     ri <- x_clr[i,]
# #     diffs <- sapply(listw$neighbours[[i]], function(j){
# #       rj <- x_clr[j,]
# #       sum((ri - rj) ** 2)
# #     })
# #     sub_sum <- sum(listw$weights[[i]] * diffs)
# #     return(sub_sum)
# #   })
# #
# #   numerator <- (nrow(x_clr) - 1) * sum(sums_ij)
# #   S0 <- sum(sapply(listw$weights, sum))
# #   x_bar <- colMeans(x_clr)
# #   diff_mean <-  sum(c(sweep(x_clr, 2, x_bar, "-")**2))
# #
# #   C <- numerator / (2*S0 * diff_mean)
# #
# # }
#
#
# spConsistency <- function(object, nblistw = NULL, window = NULL, nrep = 999, adj = FALSE, mindist = 1e-11, use_clr = FALSE) {
#
#   if(inherits(object, "FCMres")){
#     belongmat <- as.matrix(object$Belongings)
#     if(object$isRaster & is.null(window)){
#       window <- object$window
#       if(is.null(window)){
#         stop("impossible to find a window in the given object, please
#              specify one by hand.")
#       }
#     }
#     if(object$isRaster == FALSE & is.null(nblistw)){
#       nblistw <- object$nblistw
#     }
#   }else{
#     belongmat <- as.matrix(object)
#   }
#
#   # if we are not in raster mode
#
#   if(is.null(window)){
#
#     if(is.null(nblistw)){
#       stop("The nblistw must be provided if spatial vector data is used")
#     }
#     weights <- nblistw$weights
#     neighbours <- nblistw$neighbours
#     ## calcul de l'inconsistence spatiale actuelle
#
#     if(use_clr){
#       belong_mat <- clr(acomp(belongmat))
#     }
#     # we could aslo use the Aitchison distance (https://ima.udg.edu/~barcelo/index_archivos/Measures_of_difference__Clustering.pdf)
#     # this is simply done by calculating the clr transformation on the original data
#     # belongmat <- log(belongmat)
#     # belongmat <- sweep(belongmat, 1, rowMeans(belongmat), "-")
#     # belongmat <- as.matrix(belongmat)
#
#     obsdev <- sapply(1:nrow(belongmat), function(i) {
#       row <- belongmat[i, ]
#       idneighbour <- neighbours[[i]]
#       neighbour <- belongmat[idneighbour, ]
#       if (length(idneighbour) == 1){
#         neighbour <- t(as.matrix(neighbour))
#       }
#       W <- weights[[i]]
#       # we are using here the euclidean distance
#       diff <- (neighbour-row[col(neighbour)])**2 * W
#       tot <- sum(rowSums(diff))
#       return(tot)
#     })
#
#     totalcons <- sum(obsdev)
#
#     ## simulation de l'inconsistance spatiale
#     belongmat <- t(belongmat)
#     n <- ncol(belongmat)
#     simulated <- vapply(1:nrep, function(d) {
#       belong2 <- belongmat[,sample(n)]
#       simvalues <- vapply(1:ncol(belong2), function(i) {
#         row <- belong2[,i]
#         idneighbour <- neighbours[[i]]
#         neighbour <- belong2[,neighbours[[i]]]
#         if (length(idneighbour) == 1){
#           neighbour <- t(as.matrix(neighbour))
#         }
#         W <- weights[[i]]
#         diff <- (neighbour-row)
#         tot <- sum(diff^2 * W)
#         return(tot)
#       }, FUN.VALUE = 1)
#       return(sum(simvalues))
#     },FUN.VALUE = 1)
#     ratio <- totalcons / simulated
#
#     # if we are using a raster mode.
#   }else{
#     # we must calculate for each pixel its distance to its neighbours
#     # on the membership matrix. So we will calculate the distance for each
#     # raster in object$rasters and then sum them all group
#     rastnames <- names(object$rasters)
#     ok_names <- rastnames[grepl("group",rastnames, fixed = TRUE)]
#     rasters <- object$rasters[ok_names]
#     matrices <- lapply(rasters, terra::as.matrix, wide = TRUE)
#     mat_dim <- dim(matrices[[1]])
#
#     if(use_clr){
#       big_mat <- do.call(cbind,lapply(matrices, c))
#       big_mat <- clr(acomp(big_mat))
#       matrices <- lapply(1:ncol(big_mat), function(i){
#         mat <- big_mat[,i]
#         dim(mat) <- mat_dim
#         return(mat)
#       })
#     }
#
#     # applying the
#
#     if(adj){
#       dataset <- lapply(1:ncol(object$Data), function(ic){
#         vec1 <- object$Data[,ic]
#         vec2 <- rep(NA,length(object$missing))
#         vec2[object$missing] <- vec1
#         rast <- object$rasters[[1]]
#         terra::values(rast) <- vec2
#         mat <- terra::as.matrix(rast, wide = TRUE)
#         return(mat)
#
#       })
#       totalcons <- calc_raster_spinconsistency(matrices, window, adj, dataset, mindist = mindist)
#
#     }else{
#       totalcons <- calc_raster_spinconsistency(matrices,window)
#     }
#
#
#     # we must now do the same thing but with resampled values
#     warning("Calculating the permutation for the spatial inconsistency
#             when using raster can be long, depending on the raster size.
#             Note that the high number of cell in a raster reduces the need of
#             a great number of replications.")
#     # creating a vector of ids for each cell in raster
#     all_ids <- 1:terra::ncell(rasters[[1]])
#
#     # converting the matrices (columns of membership matrix) into 1d vectors
#     mem_vecs <- lapply(rasters, function(rast){
#       mat <- terra::as.matrix(rast, wide = TRUE)
#       dim(mat) <- NULL
#       return(mat)
#     })
#
#     # if necessary, doing the same with the original data
#     if(adj){
#       data_vecs <- lapply(dataset, function(mat){
#         vec <- mat
#         dim(vec) <- NULL
#         return(vec)
#       })
#     }
#     # extracting the dimension of the raster
#     #dim(terra::as.matrix(rasters[[1]]))
#
#     # starting the simulations
#     simulated <- sapply(1:nrep, function(i){
#
#       # resampling the ids
#       Ids <- sample(all_ids)
#
#       # resampling the matrices of memberships
#       new_matrices <- lapply(mem_vecs, function(vec){
#         new_vec <- vec[Ids]
#
#         # swapping the NA at their original place
#         val_na <- new_vec[!object$missing]
#         loc_na <- is.na(new_vec)
#         new_vec[!object$missing] <- NA
#         new_vec[loc_na] <- val_na
#
#         dim(new_vec) <- mat_dim
#         return(new_vec)
#       })
#
#       if(adj){
#         # resampling the matrices of the data
#         new_dataset <- lapply(data_vecs, function(vec){
#           new_vec <- vec[Ids];
#
#           # swapping the NA at their original place
#           val_na <- new_vec[!object$missing]
#           loc_na <- is.na(new_vec)
#           new_vec[!object$missing] <- NA
#           new_vec[loc_na] <- val_na
#
#           dim(new_vec) <- mat_dim
#           return(new_vec)
#         })
#         # calculating the index value
#         inconsist <-  calc_raster_spinconsistency(new_matrices, window, adj, new_dataset, mindist = mindist)
#
#       }else{
#         # calculating the index value
#         inconsist <-  calc_raster_spinconsistency(new_matrices, window)
#       }
#
#
#       return(inconsist)
#     })
#     ratio <- totalcons / simulated
#   }
#
#   return(list(Mean = mean(ratio), Median = quantile(ratio, probs = c(0.5)),
#               prt05 = quantile(ratio, probs = c(0.05)),
#               prt95 = quantile(ratio, probs = c(0.95)),
#               samples = ratio,
#               sum_diff = totalcons))
# }
