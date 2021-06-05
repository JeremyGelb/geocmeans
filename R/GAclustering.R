# # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# #### THIS SECTION IS UNDER DEVELOPMENT ####
# # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
# Gaclustering <- function(dataset, k, listw, m, inter.loss = 0.1, max.alpha = 4, nepochs = 100, nspawns = 100,
#                          tol = 0.001, maxiter = 300, standardize = TRUE){
#
#   data(LyonIris)
#   AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#   "TxChom1564","Pct_brevet","NivVieMed")
#   dataset <- LyonIris@data[AnalysisFields]
#   queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#   Wqueen <- spdep::nb2listw(queen,style="W")
#
#   listw <- Wqueen
#   nspawns <- 100
#   max.alpha <- 4
#   tol = 0.001
#   maxiter = 300
#   standardize = TRUE
#   k = 4
#   m = 1.5
#   ## step1 : calculating the spatial moran eigen values
#   spVectors <- adespatial::mem(listw)
#
#   ## step 2 : define some prior parameters
#   n <- ncol(spMat)
#   nl <- 1:n
#   nlprobs <- (1/log(nl+1)) / sum(1/log(nl+1))
#
#   l <- 1:n
#   lprobs <- rep(1/n, times = n)
#
#   ## defining the original situation
#   #select for each spawn the number of vectors to use
#   sel_nl <- sample(nl, size = nspawns, prob = nlprobs, replace = TRUE)
#   total_vecs <- sum(sel_nl)
#
#   #now select the ids of the vectors to use
#   sel_vecs <- sample(l, size = total_vecs, prob = lprobs, replace = TRUE)
#   sel_vecs <- split(sel_vecs,
#                     f = rep(1:nspawns,times=sel_nl))
#
#   #and the values of alpha
#   alphas <- seq(0.01,max.alpha, 0.01)
#   alphas_probs <- rep(1/length(alphas), times = length(alphas))
#   alpha_sel <- sample(alphas, size = nspawns, prob = alphas_probs, replace = TRUE)
#
#   situation <- list(
#     nb_factors = sel_nl,
#     factors = sel_vecs,
#     alphas = alpha_sel
#   )
#
#   ## now launching the GA
#   for(i in 1:nepochs){
#     print(paste("epoch number ", i, sep = ""))
#
#     ## calculating the results from the situation
#     results <- lapply(1:length(situation$nb_factors), function(j){
#
#       spMat <- spVectors[,situation$factors[[j]]]
#       values <- main_worker("SFCM", data = dataset, wdata = spMat,
#                             nblistw = listw, k = k, m = m, alpha = situation$alpha[[j]],
#                             lag_method="mean", maxiter = maxiter, tol = tol,
#                             standardize = standardize, verbose = FALSE,
#                             seed = NULL, init = "random")
#
#     })
#
#   }
# }
