# EN DEVELOPPEMENT
#
# calc_muik_and_uik <- function(mu, listw, centers, data, m){
#   # step1 : calculating fik
#   # should be a matrix with nrow = number of observation and ncol = number of groups
#   # it is basically a lagged version of the mu (belonging) matrix
#   fik <- apply(mu,2,function(mu_i){
#     spdep::lag.listw(listw, mu_i,)
#   })
#
#   # calculating the euclidean distance between each observation and each cluster
#   obs_dists <- do.call(cbind,lapply(1:ncol(mu), function(i){
#     ct <- centers[i,]
#     s1 <- rowSums(sweep(data, 2, ct, `-`)**2)**(2/(m-1))
#     return(s1)
#   }))
#
#   s2 <- rowSums(obs_dists)
#
#   uik <- do.call(cbind,lapply(1:ncol(mu), function(i){
#     s1 <- obs_dists[,i]
#     num <- fik[,i]
#     denom <- s1 / s2
#     return(num/denom)
#   }))
#
#   muik <- do.call(cbind,lapply(1:ncol(mu), function(i){
#     s1 <- obs_dists[,i]
#     num <- 1
#     denom <- s1 / s2
#     return(num/denom)
#   }))
#
#   return(list("uik" = uik,
#               "muik" = muik))
# }
#
#
# ## TODO :
# # defining a function to calculate belongings
# # defining a function to recalculate centers
#
#
# calcCSFCMbelongings <- function(belongings, centers, data, p, q, m, listw){
#   els <- calc_muik(belongings, listw, centers, data, m)
#
#
#   num <- (els$muik**p) * (els$uik**q)
#
# }
