# library(fclust)
# library(reldist)
# library(SDMTools)
# library(fmsb)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Fonctions de diagnostic #####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' calculate several clustering quality indexes (most of them come from fclust
#' package)
#'
#' @param data the original dataframe used fot the classification (n*p)
#' @param belongmatrix A belonging matrix (n*k)
#' @param centers A dataframe representing the centers of the clusters (k*p)
#' @return a named list with
#' \itemize{
#'         \item Silhouette.index: the silhouette index (fclust::SIL.F)
#'         \item Partition.entropy: the partition entropy index (fclust::PE)
#'         \item Partition.coeff: the partition entropy coefficient (fclust::PC)
#'         \item Modified.partition.coeff: the modified partition entropy coefficient (fclust::MPC)
#'         \item Explained.inertia: the percentage of total inertia explained by the solution
#'}
#' @export
#' @examples
#' data(LyonIris)
#' library(spdep)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' calcqualityIndexes(result$Data,result$Belongings, result$Centers)
calcqualityIndexes <- function(data, belongmatrix, centers) {
    idxsf <- fclust::SIL.F(data, belongmatrix)  #look for maximum
    idxpe <- fclust::PE(belongmatrix)  #look for minimum
    idxpc <- fclust::PC(belongmatrix)  #loook for maximum
    idxmpc <- fclust::MPC(belongmatrix)  #look for maximum
    # calculating the explained inertia
    means <- apply(data, 2, mean)
    baseinertia <- sum(calcEuclideanDistance(data, means))
    restinertia <- sapply(1:nrow(centers), function(i) {
        point <- centers[i, ]
        dists <- calcEuclideanDistance(data, point) * belongmatrix[, i]
        return(sum(dists))
    })
    explainedinertia <- 1 - (sum(restinertia) / baseinertia)

    return(list(Silhouette.index = idxsf, Partition.entropy = idxpe,
                Partition.coeff = idxpc, Modified.partition.coeff = idxmpc,
                Explained.inertia = explainedinertia))
}


#' Utility function to facilitate the spatial diagnostic of a classification
#'
#' Calculate the followinf indicators : Moran I index (spdep::moranI) for each
#' column of the belonging matrix Join count test (spdep::joincount.multi) for
#' the most likely groups of each datapoint Spatial consistency index (see
#' function spConsistency)
#'
#' @param belongmatrix A belonging matrix
#' @param nblistw A list.w  object describing the neighbours (spdep package)
#' @param undecided A float between 0 and 1 giving the minimum value that an
#'   observation must get in the belonging matrix to not be considered as
#'   uncertain (default = NULL)
#' @param nrep a integer indicating the number of permutation to do to simulate
#'   the random distribution of the spatial inconsistency
#' @return a named list with :
#' \itemize{
#'         \item MoranValues : the moran I values fo each column of the belonging
#'          matrix (spdep::MoranI)
#'         \item JoinCounts : the result of the join count test calculated with
#'          the most likely group for each datapoint (spdep::joincount.multi)
#'         \item SpConsist : the mean value of the spatial consistency index
#'         (the lower, the better, see ?spConsistency for details)
#'}
#' @export
#' @examples
#' data(LyonIris)
#' library(spdep)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- nb2listw(queen,style="W")
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
    # attribution des groupes
    groups <- colnames(belongmatrix)[max.col(belongmatrix, ties.method = "first")]
    if (is.null(undecided) == FALSE) {
        Maximums <- do.call(pmax, belongmatrix)
        groups[Maximums < undecided] <- "undecided"
    }
    groups <- as.factor(groups)
    # calcul des join count test
    spjctetst <- spdep::joincount.multi(groups, nblistw, zero.policy = TRUE)
    return(list(MoranValues = morandf, JoinCounts = spjctetst, SpConsist = Consist$Mean))
}


#' Calculate a spatial consistency index
#'
#' This index is experimental, it aims to measure how much a clustering solution
#' is spatially consistent. A classification is spatially inconsistent if
#' neighbouring observation do not belong to the same group. See detail for
#' a description of its calculation
#'
#'
#' The total spatial inconsistency (*isp*) is calculated as follow
#'
#' \deqn{isp = \sum_{i}\sum_{j}\sum_{k} (A_{ik} - A_{jk})^{2} * W_{ij}}
#'
#' With A the belonging matrix, i an observation, k the neighbours of i and W
#' the spatial weight matrix This represents the total spatial inconsistency of
#' the solution (true inconsistency) We propose to compare this total with
#' simulated values obtained by permutations (simulated inconsistency). The
#' values obtained by permutation are an approximation of the spatial
#' inconsistency obtained in a random context Ratios between the true
#' inconsistency and simulated inconsistencies are calculated A value of 0
#' depict a situation where all observations are identitical to their neighbours
#' A value of 1 depict a situation where all observations are as much different
#' as their neighbours that what randomness can produce A classification
#' solution able to reduce this index has a better spatial consistency
#'
#' @param belongmatrix A belonging matrix
#' @param nblistw A list.w object describing the neighbours (spdep package)
#'   observation must get in the belonging matrix to not be considered as
#'   uncertain (default = NULL)
#' @param nrep a integer indicating the number of permutation to do to simulate
#' @return a named list with
#'  \itemize{
#'         \item Mean : the mean of the spatial consistency index
#'         \item prt05 : the 5th percentile of the spatial consistency index
#'         \item prt95 : the 95th percentule of the spatial consistency index
#'         \item samples : all the value of the spatial consistency index
#' }
#' @export
#' @examples
#' data(LyonIris)
#' library(spdep)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' spConsistency(result$Belongings, Wqueen, nrep=50)
spConsistency <- function(belongmatrix, nblistw, nrep = 999) {
    ## calcul de l'inconsistence spatiale actuelle
    obsdev <- sapply(1:nrow(belongmatrix), function(i) {
        row <- belongmatrix[i, ]
        neighbours <- belongmatrix[nblistw$neighbours[[i]], ]
        W <- nblistw$weights[[i]]
        diff <- t((t(neighbours) - row)^2)
        tot <- sum(rowSums(diff) * W)
        return(tot)
    })
    totalcons <- sum(obsdev)
    ## simulation de l'inconsistance spatiale
    simulated <- sapply(1:nrep, function(d) {
        belong2 <- belongmatrix[sample(nrow(belongmatrix)), ]
        simvalues <- sapply(1:nrow(belong2), function(i) {
            row <- belong2[i, ]
            neighbours <- belong2[nblistw$neighbours[[i]], ]
            W <- nblistw$weights[[i]]
            diff <- t((t(neighbours) - row)^2)
            tot <- sum(rowSums(diff) * W)
            return(tot)
        })
        return(sum(simvalues))
    })
    ratio <- totalcons / simulated
    return(list(Mean = mean(ratio), Median = quantile(ratio, probs = c(0.5)),
                prt05 = quantile(ratio, probs = c(0.05)),
                prt95 = quantile(ratio, probs = c(0.95)),
                samples = ratio))
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Fonctions de visualisation #####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' build some maps to visualize the results of the clustering
#'
#' @param geodata A object of class spatialpolygonesdataframe /
#' spatiallinesdataframe or spatialpointsdataframe ordered
#' like the original data used for the clustering
#' @param belongmatrix The belonging matrix obtained at the end of the algorithm
#' @param undecided A float between 0 and 1 giving the minimum value that an
#'   observation must get in the belonging matrix to not be considered as
#'   uncertain (default = NULL)
#' @return a named list with :
#' \itemize{
#'         \item ProbaMaps : a list of ggplot maps showing for each group the
#'         probability of the observations to belong to that group
#'         \item ClusterMap : a ggplot map showing the most likely group for
#'          observation
#' }
#' @export
#' @examples
#' data(LyonIris)
#' library(spdep)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' MyMaps <- mapClusters(LyonIris, result$Belongings)
mapClusters <- function(geodata, belongmatrix, undecided = NULL) {
    if(class(geodata)[[1]]=="SpatialPolygonsDataFrame"){
        return(mapPolygons(geodata, belongmatrix, undecided))
    }else if(class(geodata)[[1]]=="SpatialPointsDataFrame"){
        return(mapPoints(geodata, belongmatrix, undecided))
    }else if(class(geodata)[[1]]=="SpatialLinesDataFrame"){
        return(mapLines(geodata, belongmatrix, undecided))
    }else {
        stop("The object passed in geodata argument is not supported.
              Supported classes are : SpatialPolygonsDataFrame,
              SpatialPointsDataFrame and SpatialLinesDataFrame")
    }
}

#' Internal function to realize maps based on spatialpolygondataframe
#'
#' @param geodata A spatialpolygonsdataframe ordered like the original data used
#'   for the clustering
#' @param belongmatrix The belonging matrix obtained at the end of the algorithm
#' @param undecided A float between 0 and 1 giving the minimum value that an
#'   observation must get in the belonging matrix to not be considered as
#'   uncertain (default = NULL)
#' @return a named list with :
#' \itemize{
#'         \item ProbaMaps : a list of ggplot maps showing for each group the
#'         probability of the observations to belong to that group
#'         \item ClusterMap : a ggplot map showing the most likely group for each observation
#' }
#' @examples
#' #No example provided, this is an internal function, use the general wrapper function mapClusters
mapPolygons <- function(geodata, belongmatrix, undecided = NULL){
    belongmatrix <- as.data.frame(belongmatrix)
    names(belongmatrix) <- gsub(" ", "", names(belongmatrix), fixed = T)
    geodata@data <- cbind(geodata@data, belongmatrix)
    geodata$OID <- 1:nrow(geodata@data)

    # attribution des groupes
    Groups <- names(belongmatrix)[max.col(belongmatrix, ties.method = "first")]
    if (is.null(undecided) == FALSE) {
        Maximums <- do.call(pmax, belongmatrix)
        Undecided <- ifelse(Maximums < undecided, "Undecided", "Ok")
    } else {
        Undecided <- rep("Ok", length(Groups))
    }
    geodata$Cluster <- Groups
    geodata$Undecided <- Undecided
    FortiData <- broom::tidy(geodata, region = "OID")
    FortiData <- merge(FortiData, geodata@data, by.x = "id", by.y = "OID")
    # realisation des cartes de probabilites
    ProbaPlots <- lapply(names(belongmatrix), function(Name) {
        Plot <- ggplot2::ggplot(FortiData) +
            ggplot2::geom_polygon(ggplot2::aes_string(x = "long", y = "lat", group = "group", fill = Name),
                                  colour = "black", size = 0.1) +
            ggplot2::scale_fill_gradient(low = "white", high = "blue") +
            ggplot2::coord_fixed(ratio = 1)+
            ggplot2::theme(
                axis.title = ggplot2::element_blank(),
                axis.text = ggplot2::element_blank(),
                axis.ticks = ggplot2::element_blank()
            )
        return(Plot)
    })
    # realisation de la carte des groupes
    ClusterMap <- ggplot2::ggplot(FortiData) +
        ggplot2::geom_polygon(ggplot2::aes(x = long, y = lat, group = group, fill = Cluster), size = 0.1) +
        ggplot2::scale_fill_brewer(palette = "Set1") +
        ggplot2::geom_polygon(ggplot2::aes(x = long, y = lat, group = group),
                              fill = rgb(1, 1, 1, 0.7),
                              size = 0.1,
                              data = subset(FortiData,FortiData$Undecided == "Undecided")) +
        ggplot2::coord_fixed(ratio = 1)+
        ggplot2::theme(
            axis.title = ggplot2::element_blank(),
            axis.text = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank()
        )
    return(list(ProbaMaps = ProbaPlots, ClusterPlot = ClusterMap))
}


#' Internal function to realize maps based on spatiallinesdataframe
#'
#' @param geodata A spatiallinesdataframe ordered like the original data used
#'   for the clustering
#' @param belongmatrix The belonging matrix obtained at the end of the algorithm
#' @param undecided A float between 0 and 1 giving the minimum value that an
#'   observation must get in the belonging matrix to not be considered as
#'   uncertain (default = NULL)
#' @return a named list with :
#' \itemize{
#'         \item ProbaMaps : a list of ggplot maps showing for each group the
#'         probability of the observations to belong to that group
#'         \item ClusterMap : a ggplot map showing the most likely group for each observation
#' }
#' @examples
#' #No example provided, this is an internal function, use the general wrapper function mapClusters
mapLines <- function(geodata, belongmatrix, undecided = NULL){
    belongmatrix <- as.data.frame(belongmatrix)
    names(belongmatrix) <- gsub(" ", "", names(belongmatrix), fixed = T)
    geodata@data <- cbind(geodata@data, belongmatrix)
    geodata$OID <- 1:nrow(geodata@data)

    # attribution des groupes
    Groups <- names(belongmatrix)[max.col(belongmatrix, ties.method = "first")]
    if (is.null(undecided) == FALSE) {
        Maximums <- do.call(pmax, belongmatrix)
        Undecided <- ifelse(Maximums < undecided, "Undecided", "Ok")
    } else {
        Undecided <- rep("Ok", length(Groups))
    }
    geodata$Cluster <- Groups
    geodata$Undecided <- Undecided
    FortiData <- broom::tidy(geodata, region = "OID")
    FortiData <- merge(FortiData, geodata@data, by.x = "id", by.y = "OID")
    # realisation des cartes de probabilites
    ProbaPlots <- lapply(names(belongmatrix), function(Name) {
        Plot <- ggplot2::ggplot(FortiData) +
            ggplot2::geom_path(ggplot2::aes_string(x = "long", y = "lat", group = "group", color = Name),
                                  size = 0.1) +
            ggplot2::scale_fill_gradient(low = "white", high = "blue") +
            ggplot2::coord_fixed(ratio = 1)+
            ggplot2::theme(
                axis.title = ggplot2::element_blank(),
                axis.text = ggplot2::element_blank(),
                axis.ticks = ggplot2::element_blank()
            )
        return(Plot)
    })
    # realisation de la carte des groupes
    ClusterMap <- ggplot2::ggplot(FortiData) +
        ggplot2::geom_path(ggplot2::aes(x = long, y = lat, group = group, color = Cluster), size = 0.1) +
        ggplot2::scale_fill_brewer(palette = "Set1") +
        ggplot2::geom_path(ggplot2::aes(x = long, y = lat, group = group),
                              color = rgb(1, 1, 1, 0.7),
                              size = 0.1,
                              data = subset(FortiData,FortiData$Undecided == "Undecided")) +
        ggplot2::coord_fixed(ratio = 1)+
        ggplot2::theme(
            axis.title = ggplot2::element_blank(),
            axis.text = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank()
        )
    return(list(ProbaMaps = ProbaPlots, ClusterPlot = ClusterMap))
}


#' Internal function to realize maps based on spatialpolygondataframe
#'
#' @param geodata A spatialpointsdataframe ordered like the original data used
#'   for the clustering
#' @param belongmatrix The belonging matrix obtained at the end of the algorithm
#' @param undecided A float between 0 and 1 giving the minimum value that an
#'   observation must get in the belonging matrix to not be considered as
#'   uncertain (default = NULL)
#' @return a named list with :
#' \itemize{
#'         \item ProbaMaps : a list of ggplot maps showing for each group the probability
#'          of the observations to belong to that group
#'         \item ClusterMap : a ggplot map showing the most likely group for each observation
#' }
#' @examples
#' #No example provided, this is an internal function, use the general wrapper function mapClusters
mapPoints <- function(geodata, belongmatrix, undecided = NULL){
    belongmatrix <- as.data.frame(belongmatrix)
    names(belongmatrix) <- gsub(" ", "", names(belongmatrix), fixed = T)
    geodata@data <- cbind(geodata@data, belongmatrix)
    Coords <- sp::coordinates(geodata)
    geodata$X__ <- Coords[,1]
    geodata$Y__ <- Coords[,2]

    # attribution des groupes
    Groups <- names(belongmatrix)[max.col(belongmatrix, ties.method = "first")]
    if (is.null(undecided) == FALSE) {
        Maximums <- do.call(pmax, belongmatrix)
        Undecided <- ifelse(Maximums < undecided, "Undecided", "Ok")
    } else {
        Undecided <- rep("Ok", length(Groups))
    }
    geodata$Cluster <- Groups
    geodata$Undecided <- Undecided
    # realisation des cartes de probabilites
    ProbaPlots <- lapply(names(belongmatrix), function(Name) {
        Plot <- ggplot2::ggplot(geodata@data) +
            ggplot2::geom_point(ggplot2::aes_string(x = "X__", y = "Y__", color = Name),
                                  size = 1) +
            ggplot2::scale_fill_gradient(low = "white", high = "blue") +
            ggplot2::coord_fixed(ratio = 1)+
            ggplot2::theme(
                axis.title = ggplot2::element_blank(),
                axis.text = ggplot2::element_blank(),
                axis.ticks = ggplot2::element_blank()
            )
        return(Plot)
    })
    # realisation de la carte des groupes
    ClusterMap <- ggplot2::ggplot(geodata@data) +
        ggplot2::geom_point(ggplot2::aes(x = X__, y = Y__, color = Cluster), size = 1) +
        ggplot2::scale_fill_brewer(palette = "Set1") +
        ggplot2::geom_point(ggplot2::aes(x = X__, y = Y__),
                              fill = rgb(1, 1, 1, 0.7),
                              size = 1,
                              data = subset(geodata@data,geodata@data$Undecided == "Undecided")) +
        ggplot2::coord_fixed(ratio = 1)+
        ggplot2::theme(
            axis.title = ggplot2::element_blank(),
            axis.text = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank()
            )
    return(list(ProbaMaps = ProbaPlots, ClusterPlot = ClusterMap))
}


#' @importFrom dplyr %>%
NULL
#' Calculate some descriptive statistics of each group
#'
#' @param data the original dataframe used fot the classification
#' @param belongmatrix A belonging matrix
#' @param weighted A boolean indicating if the summary statistics must use the
#'   belonging matrix columns as weights (TRUE) or simply assign each
#'   observation to its most likely cluster and compute the statistics on each
#'   subset (default=True)
#' @param dec a integer indicating the number of digits to keep when rouding (default is 3)
#' @return a list of length k (the number of group). Each element of the list is
#'   a dataframe with summary statistics for the variables of data for each
#'   group
#' @export
#' @examples
#' data(LyonIris)
#' library(spdep)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' summarizeClusters(dataset, result$Belongings)
summarizeClusters <- function(data, belongmatrix, weighted = TRUE, dec = 3) {
    belongmatrix <- as.data.frame(belongmatrix)
    if (weighted) {
        Summaries <- lapply(1:ncol(belongmatrix), function(c) {
            W <- belongmatrix[, c]
            Values <- apply(data, 2, function(x) {
                Q5 <- round(reldist::wtd.quantile(x, q = 0.05, na.rm = TRUE, weight = W),dec)
                Q10 <- round(reldist::wtd.quantile(x, q = 0.1, na.rm = TRUE, weight = W),dec)
                Q25 <- round(reldist::wtd.quantile(x, q = 0.25, na.rm = TRUE, weight = W),dec)
                Q50 <- round(reldist::wtd.quantile(x, q = 0.5, na.rm = TRUE, weight = W),dec)
                Q75 <- round(reldist::wtd.quantile(x, q = 0.75, na.rm = TRUE, weight = W),dec)
                Q90 <- round(reldist::wtd.quantile(x, q = 0.9, na.rm = TRUE, weight = W),dec)
                Q95 <- round(reldist::wtd.quantile(x, q = 0.95, na.rm = TRUE, weight = W),dec)
                Mean <- round(weighted.mean(x, W),dec)
                Std <- round(sqrt(reldist::wtd.var(x, weight=W)),dec)
                return(list(Q5 = Q5, Q10 = Q10, Q25 = Q25, Q50 = Q50,
                            Q75 = Q75, Q90 = Q90, Q95 = Q95,
                            Mean = Mean, Std = Std))
            })
            DF <- do.call(cbind, Values)
            print(paste("Statistic summary for cluster ", c, sep = ""))
            print(DF)
            return(DF)
        })
        names(Summaries) <- paste("Cluster_", c(1:ncol(belongmatrix)), sep = "")
        return(Summaries)

    } else {
        Groups <- colnames(belongmatrix)[max.col(belongmatrix, ties.method = "first")]
        data$Groups <- Groups
        Summaries <- lapply(unique(data$Groups), function(c) {
            DF <- subset(data, data$Groups == c)
            DF$Groups <- NULL
            Values <- apply(DF, 2, function(x) {
                Q5 <- round(quantile(x, probs = 0.05, na.rm = TRUE),dec)
                Q10 <- round(quantile(x, probs = 0.1, na.rm = TRUE),dec)
                Q25 <- round(quantile(x, probs = 0.25, na.rm = TRUE),dec)
                Q50 <- round(quantile(x, probs = 0.5, na.rm = TRUE),dec)
                Q75 <- round(quantile(x, probs = 0.75, na.rm = TRUE),dec)
                Q90 <- round(quantile(x, probs = 0.9, na.rm = TRUE),dec)
                Q95 <- round(quantile(x, probs = 0.95, na.rm = TRUE),dec)
                Mean <- round(mean(x),dec)
                Std <- round(sd(x),dec)
                return(list(Q5 = Q5, Q10 = Q10, Q25 = Q25, Q50 = Q50,
                            Q75 = Q75, Q90 = Q90, Q95 = Q95,
                            Mean = Mean, Std = Std))
            })
            DF <- do.call(cbind, Values)
            print(paste("Statistic summary for cluster ", c, sep = ""))
            print(DF)
            return(DF)
        })
        names(Summaries) <- paste("Cluster_", c(1:ncol(belongmatrix)), sep = "")
        return(Summaries)
    }
}


#' Identify the observation for with the classification is uncertain
#'
#' @param belongmatrix The belonging matrix obtained at the end of the algorithm
#' @param tol A float indicating the minimum required level of belonging to be
#'   not considered as uncertaint
#' @return a vector indicating the most liekly group for each observation or
#'   "Undecided" if the maximum probability for the observation does not reach
#'   the value of the tol parameter
#' @export
#' @examples
#' data(LyonIris)
#' library(spdep)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' undecidedUnits(result$Belongings, tol = 0.45)
undecidedUnits <- function(belongmatrix, tol = 0.1) {
    belongmatrix <- as.data.frame(belongmatrix)
    groups <- colnames(belongmatrix)[max.col(belongmatrix, ties.method = "first")]
    rowMax <- do.call(pmax, belongmatrix)
    DF <- data.frame(groups = groups, maxprob = rowMax)
    return(ifelse(DF$maxprob < tol, "Undecided", DF$groups))
}


#' @importFrom dplyr %>%
NULL
#' display some descriptive informations about fixed groups
#'
#' @param data a dataframe with numeric columns
#' @param groupvar the name of the categorical variable used to split the
#'   dataset
#' @param vars the names of the numerical variables
#' @param dec the number of digits to keep
#' @return a named list with :
#' \itemize{
#'         \item Violin : a list with a length equal to the number  of columns in data.
#'         Each element of the list in a violin plot usefull to compare the distributions
#'         of each group on that variable
#'         \item Tableau : a list with a length equal to the number of group (k). Each
#'         element in a dataframe with summary statistics describing the group for the
#'         columns of data.
#' }
#' @export
#' @examples
#' data(LyonIris)
#' library(spdep)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' dataset$group <- result$Groups
#' describGroups(dataset, "group", AnalysisFields)
describGroups <- function(data, groupvar, vars, dec=5) {
    #%%%%%%%%%%% sortir des violin plots %%%%%%%%%%%
    ViolinPlots <- list()
    for (Var in vars) {
        Plot <- ggplot2::ggplot(data, ggplot2::aes_string(x = groupvar, y = Var, fill = groupvar)) +
            ggplot2::geom_violin(show.legend = FALSE) +
            ggplot2::geom_boxplot(width = 0.1, fill = "white", show.legend = FALSE) +
            ggplot2::scale_color_brewer(palette = "Dark2")
        ViolinPlots[[length(ViolinPlots) + 1]] <- Plot
    }

    # %%%%%%%%%%% sortir spyderplots %%%%%%%%%%%
    Groups <- unique(data[[groupvar]])

    Values <- do.call(rbind, lapply(Groups, function(x) {
        Sub <- subset(data, data[[groupvar]] == x)
        return(data.frame(matrix(sapply(Sub[vars], mean), ncol = length(vars))))
    }))
    Mins <- apply(Values, 2, min)
    Maxs <- apply(Values, 2, max)


    for (Gp in Groups) {
        Sub <- subset(data, data[[groupvar]] == Gp)
        Scores <- data.frame(matrix(sapply(Sub[vars], mean), ncol = length(vars)))
        names(Scores) <- vars
        datam <- rbind(Maxs, Mins, Scores)
        Chart <- fmsb::radarchart(datam, axistype = 1, pcol = rgb(0.2, 0.5, 0.5, 0.9),
                            pfcol = rgb(0.2, 0.5, 0.5, 0.5),
                            plwd = 4, cglcol = "grey", cglty = 1,
                            axislabcol = "grey", cglwd = 0.8, vlcex = 0.8,
                            title = paste("Group number : ",Gp))
    }

    # %%%%%%%%%%% sortir un tableau des stats descriptives %%%%%%%%%%%

    Q1 <- function(x) {
        return(quantile(x, probs = 0.25))
    }
    Q2 <- function(x) {
        return(quantile(x, probs = 0.5))
    }
    Q3 <- function(x) {
        return(quantile(x, probs = 0.75))
    }

    # calculating summary statistics
    TableauDes <- data %>% dplyr::group_by_at(.vars = groupvar) %>%
        dplyr::summarize_at(.vars = vars, .funs = c(min, max, mean, Q1, Q2, Q3, sd))

    # renaming columns
    NewNames <- c()
    for (fun in c("min", "max", "mean", "Q1", "Q2", "Q3", "sd")) {
        for (Var in vars) {
            NewNames[[length(NewNames) + 1]] <- paste(Var, fun, sep = " : ")
        }
    }

    names(TableauDes) <- c("Groupe", NewNames)
    for (N in NewNames) {
        TableauDes[[N]] <- round(TableauDes[[N]], dec)
    }

    # formatting the final table
    TableauDes <- as.data.frame(t(TableauDes))
    names(TableauDes) <- paste("Groupe_", 1:ncol(TableauDes), sep = "")
    TableauDes <- TableauDes[order(row.names(TableauDes)), ]
    TableauDes <- TableauDes[row.names(TableauDes) != "Groupe", ]

    return(list(Violin = ViolinPlots, Tableau = TableauDes))
}
