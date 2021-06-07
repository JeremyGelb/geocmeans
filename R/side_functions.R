# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Fonctions de visualisation #####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Mapping the clusters
#'
#' @description Build some maps to visualize the results of the clustering
#'
#' @param geodata A object of class spatialpolygonesdataframe /
#' spatiallinesdataframe or spatialpointsdataframe ordered
#' like the original data used for the clustering
#' @param belongmatrix The membership matrix obtained at the end of the algorithm
#' @param undecided A float between 0 and 1 giving the minimum value that an
#'   observation must get in the membership matrix to not be considered as
#'   uncertain (default = NULL)
#' @return A named list with :
#' \itemize{
#'         \item ProbaMaps : a list of ggplot maps showing for each group the
#'         probability of the observations to belong to that group
#'         \item ClusterMap : a ggplot map showing the most likely group for
#'          observation
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
#' MyMaps <- mapClusters(LyonIris, result$Belongings)
mapClusters <- function(geodata, belongmatrix, undecided = NULL) {
    geodata$OID <- as.character(1:nrow(geodata@data))
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


#' @title Mapping the clusters (polygons)
#'
#' @description Internal function to realize maps based on spatialpolygondataframe
#'
#' @param geodata A spatialpolygonsdataframe ordered like the original data used
#'   for the clustering
#' @param belongmatrix The membership matrix obtained at the end of the algorithm
#' @param undecided A float between 0 and 1 giving the minimum value that an
#'   observation must get in the membership matrix to not be considered as
#'   uncertain (default = NULL)
#' @return A named list with :
#' \itemize{
#'         \item ProbaMaps : a list of ggplot maps showing for each group the
#'         probability of the observations to belong to that group
#'         \item ClusterMap : a ggplot map showing the most likely group for each observation
#' }
#' @keywords internal
#' @examples
#' #No example provided, this is an internal function, use the general wrapper function mapClusters
mapPolygons <- function(geodata, belongmatrix, undecided = NULL){
    belongmatrix <- as.data.frame(belongmatrix)
    names(belongmatrix) <- gsub(" ", "", names(belongmatrix), fixed = TRUE)
    geodata@data <- cbind(geodata@data,belongmatrix)

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
    FortiData <- ggplot2::fortify(geodata, region = "OID")
    FortiData$OID <- FortiData$id
    geodata$OID <- rownames(geodata@data)
    FortiData <- dplyr::left_join(FortiData, geodata@data, by = "OID")
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
        ggplot2::geom_polygon(ggplot2::aes_string(x = "long", y = "lat", group = "group", fill = "Cluster"), size = 0.1) +
        ggplot2::scale_fill_brewer(palette = "Set1") +
        ggplot2::geom_polygon(ggplot2::aes_string(x = "long", y = "lat", group = "group"),
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


#' @title Mapping the clusters (lines)
#'
#' @description Internal function to realize maps based on spatiallinesdataframe
#'
#' @param geodata A spatiallinesdataframe ordered like the original data used
#'   for the clustering
#' @param belongmatrix The membership matrix obtained at the end of the algorithm
#' @param undecided A float between 0 and 1 giving the minimum value that an
#'   observation must get in the membership matrix to not be considered as
#'   uncertain (default = NULL)
#' @return A named list with :
#' \itemize{
#'         \item ProbaMaps : a list of ggplot maps showing for each group the
#'         probability of the observations to belong to that group
#'         \item ClusterMap : a ggplot map showing the most likely group for each observation
#' }
#' @keywords internal
#' @examples
#' #No example provided, this is an internal function, use the general wrapper function mapClusters
mapLines <- function(geodata, belongmatrix, undecided = NULL){
    belongmatrix <- as.data.frame(belongmatrix)
    names(belongmatrix) <- gsub(" ", "", names(belongmatrix), fixed = TRUE)
    geodata@data <- cbind(geodata@data, belongmatrix)

    # attribution des groupes
    Groups <- names(belongmatrix)[max.col(belongmatrix, ties.method = "first")]
    if (is.null(undecided) == FALSE) {
        Maximums <- do.call(pmax, belongmatrix)
        Undecided <- ifelse(Maximums < undecided, "Undecided", "Ok")
    } else {
        Undecided <- rep("Ok", length(Groups))
    }
    geodata$Cluster <- Groups
    FortiData <- ggplot2::fortify(geodata, region = "OID")
    FortiData$OID <- FortiData$id
    geodata$OID <- rownames(geodata@data)
    FortiData <- dplyr::left_join(FortiData, geodata@data, by = "OID")
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
        ggplot2::geom_path(ggplot2::aes_string(x = "long", y = "lat", group = "group", color = "Cluster"), size = 0.1) +
        ggplot2::scale_fill_brewer(palette = "Set1") +
        ggplot2::geom_path(ggplot2::aes_string(x = "long", y = "lat", group = "group"),
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


#' @title Mapping the clusters (points)
#'
#' @description Internal function to realize maps based on spatialpointsdataframe
#'
#' @param geodata A spatialpointsdataframe ordered like the original data used
#'   for the clustering
#' @param belongmatrix The membership matrix obtained at the end of the algorithm
#' @param undecided A float between 0 and 1 giving the minimum value that an
#'   observation must get in the membership matrix to not be considered as
#'   uncertain (default = NULL)
#' @return A named list with :
#' \itemize{
#'         \item ProbaMaps : a list of ggplot maps showing for each group the probability
#'          of the observations to belong to that group
#'         \item ClusterMap : a ggplot map showing the most likely group for each observation
#' }
#' @keywords internal
#' @examples
#' #No example provided, this is an internal function, use the general wrapper function mapClusters
mapPoints <- function(geodata, belongmatrix, undecided = NULL){
    belongmatrix <- as.data.frame(belongmatrix)
    names(belongmatrix) <- gsub(" ", "", names(belongmatrix), fixed = TRUE)
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
        ggplot2::geom_point(ggplot2::aes_string(x = "X__", y = "Y__", color = "Cluster"), size = 1) +
        ggplot2::scale_fill_brewer(palette = "Set1") +
        ggplot2::geom_point(ggplot2::aes_string(x = "X__", y = "Y__"),
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


#' @title Descriptive statistics by group
#'
#' @description Calculate some descriptive statistics of each group
#'
#' @param data The original dataframe used for the classification
#' @param belongmatrix A membership matrix
#' @param weighted A boolean indicating if the summary statistics must use the
#'   membership matrix columns as weights (TRUE) or simply assign each
#'   observation to its most likely cluster and compute the statistics on each
#'   subset (FALSE)
#' @param dec An integer indicating the number of digits to keep when rounding
#' (default is 3)
#' @param silent A boolean indicating if the results must be printed or silently returned
#' @return A list of length k (the number of group). Each element of the list is
#'   a dataframe with summary statistics for the variables of data for each
#'   group
#' @export
#' @importFrom dplyr %>%
#' @importFrom grDevices rgb
#' @importFrom stats quantile sd weighted.mean
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' summarizeClusters(dataset, result$Belongings)
summarizeClusters <- function(data, belongmatrix, weighted = TRUE, dec = 3, silent=TRUE) {
    belongmatrix <- as.data.frame(belongmatrix)
    if (weighted) {
        Summaries <- lapply(1:ncol(belongmatrix), function(c) {
            W <- belongmatrix[, c]
            Values <- apply(data, 2, function(x) {
                Q5 <- as.numeric(round(reldist::wtd.quantile(x, q = 0.05, na.rm = TRUE, weight = W),dec))
                Q10 <- as.numeric(round(reldist::wtd.quantile(x, q = 0.1, na.rm = TRUE, weight = W),dec))
                Q25 <- as.numeric(round(reldist::wtd.quantile(x, q = 0.25, na.rm = TRUE, weight = W),dec))
                Q50 <- as.numeric(round(reldist::wtd.quantile(x, q = 0.5, na.rm = TRUE, weight = W),dec))
                Q75 <- as.numeric(round(reldist::wtd.quantile(x, q = 0.75, na.rm = TRUE, weight = W),dec))
                Q90 <- as.numeric(round(reldist::wtd.quantile(x, q = 0.9, na.rm = TRUE, weight = W),dec))
                Q95 <- as.numeric(round(reldist::wtd.quantile(x, q = 0.95, na.rm = TRUE, weight = W),dec))
                Mean <- as.numeric(round(weighted.mean(x, W),dec))
                Std <- round(sqrt(reldist::wtd.var(x, weight=W)),dec)
                N <-
                return(list(Q5 = Q5, Q10 = Q10, Q25 = Q25, Q50 = Q50,
                            Q75 = Q75, Q90 = Q90, Q95 = Q95,
                            Mean = Mean, Std = Std))
            })
            DF <- do.call(cbind, Values)
            if(silent==FALSE){
                print(paste("Statistic summary for cluster ", c, sep = ""))
                print(DF)
            }

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
                Q5 <- as.numeric(round(quantile(x, probs = 0.05, na.rm = TRUE),dec))
                Q10 <- as.numeric(round(quantile(x, probs = 0.1, na.rm = TRUE),dec))
                Q25 <- as.numeric(round(quantile(x, probs = 0.25, na.rm = TRUE),dec))
                Q50 <- as.numeric(round(quantile(x, probs = 0.5, na.rm = TRUE),dec))
                Q75 <- as.numeric(round(quantile(x, probs = 0.75, na.rm = TRUE),dec))
                Q90 <- as.numeric(round(quantile(x, probs = 0.9, na.rm = TRUE),dec))
                Q95 <- as.numeric(round(quantile(x, probs = 0.95, na.rm = TRUE),dec))
                Mean <- round(mean(x),dec)
                Std <- round(sd(x),dec)
                return(list(Q5 = Q5, Q10 = Q10, Q25 = Q25, Q50 = Q50,
                            Q75 = Q75, Q90 = Q90, Q95 = Q95,
                            Mean = Mean, Std = Std))
            })
            DF <- do.call(cbind, Values)
            if(silent==FALSE){
                print(paste("Statistic summary for cluster ", c, sep = ""))
                print(DF)
            }

            return(DF)
        })
        names(Summaries) <- paste("Cluster_", c(1:ncol(belongmatrix)), sep = "")
        return(Summaries)
    }
}


#' @title Summary method for FCMres
#'
#' @description Calculate some descriptive statistics of each group of a
#' FCMres object
#'
#' @param object A FCMres object, typically obtained from functions CMeans, GCMeans, SFCMeans, SGFCMeans
#' @param data A dataframe to use for the summary statistics instead of obj$data
#' @param weighted A boolean indicating if the summary statistics must use the
#'   membership matrix columns as weights (TRUE) or simply assign each
#'   observation to its most likely cluster and compute the statistics on each
#'   subset (FALSE)
#' @param dec An integer indicating the number of digits to keep when rounding
#' (default is 3)
#' @param silent A boolean indicating if the results must be printed or silently returned
#' @param ... Not used
#' @return A list of length k (the number of group). Each element of the list is
#'   a dataframe with summary statistics for the variables of data for each
#'   group
#' @export
#' @importFrom dplyr %>%
#' @importFrom grDevices rgb
#' @importFrom stats quantile sd weighted.mean
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' summary(result)
summary.FCMres <- function(object, data = NULL, weighted = TRUE, dec = 3, silent=TRUE, ...) {
    if(is.null(data)){
        df <- object$Data
    }else{
        df <- as.matrix(data)
    }
    return(summarizeClusters(df, object$Belongings, weighted, dec, silent))
}


#' @title Spider chart
#'
#' @description Display spider charts to quickly compare values between groups
#'
#' @details For each group, the weighted mean of each variable in data is calculated
#' based on the probability of belonging to this group of each observation.
#' On the chart the exterior ring represents the maximum value obtained for
#' all the groups and the interior ring the minimum. The groups are located
#' between these two limits in a linear way.
#'
#' @param data A dataframe with numeric columns
#' @param belongmatrix A membership matrix
#' @param chartcolors A vector of color names used for the spider plot
#' @return NULL, the plots are displayed directly by the function (see fmsb::radarchart)
#' @importFrom grDevices rgb colors col2rgb
#' @importFrom graphics legend
#' @importFrom stats weighted.mean
#' @export
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' spiderPlots(dataset,result$Belongings)
spiderPlots<- function(data, belongmatrix, chartcolors=NULL){
    Groups <- ncol(belongmatrix)

    Values <- do.call(rbind, lapply(1:Groups, function(i) {
        W <- belongmatrix[,i]
        return(apply(data,2, function(row){return(weighted.mean(row,W))}))
    }))
    Mins <- apply(Values, 2, min)
    Maxs <- apply(Values, 2, max)


    for (Gp in 1:Groups) {
        Scores <- Values[Gp,]
        names(Scores) <- names(data)
        datam <- data.frame(rbind(Maxs, Mins, Scores))
        fmsb::radarchart(datam, axistype = 1, pcol = rgb(0.2, 0.5, 0.5, 0.9),
                                  pfcol = rgb(0.2, 0.5, 0.5, 0.5),
                                  plwd = 4, cglcol = "grey", cglty = 1,
                                  axislabcol = "grey", cglwd = 0.8, vlcex = 0.8,
                                  title = paste("Group number : ",Gp))
    }
    if (is.null(chartcolors)){
        selcolors <- sample(colors(),size = ncol(belongmatrix))
    }else{
        selcolors <- chartcolors
    }

    tocolors <- sapply(selcolors,function(i){
        v1 <- as.list(col2rgb(i,alpha=TRUE))
        v1[[length(v1)]] <- 100
        v1$maxColorValue <- 255
        new_color <- do.call(rgb,v1)
        return(new_color)
        })
    alldata <- data.frame(rbind(Maxs, Mins,Values))
    fmsb::radarchart(alldata,pcol=selcolors,pfcol = tocolors,
                     axislabcol = "grey",plwd = 2, cglcol = "grey",
                     cglty = 1,vlcex = 0.8,
                     plty = rep(1,ncol(belongmatrix)))
    legend("right",legend = levels(as.factor(paste("Group ",1:ncol(belongmatrix),sep=""))),
           fill=tocolors)

}



#' @title Violin plots
#'
#' @description Return violin plots to compare the distribution of each variable for each
#' group.
#'
#' @param data A dataframe with numeric columns
#' @param groups A vector indicating the group of each observation
#' @return A list of plots created with ggplot2
#' @importFrom dplyr %>%
#' @importFrom grDevices rgb
#' @export
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' violinPlots(dataset, result$Groups)
violinPlots <- function(data,groups){
    data$groups <- groups
    groupvar <- "groups"
    Plots <- list()
    vars <- names(data)
    for (Var in vars) {
        if(Var!="groups"){
            Plot <- ggplot2::ggplot(data, ggplot2::aes_string(x = groupvar, y = Var, fill = groupvar)) +
                ggplot2::geom_violin(show.legend = FALSE) +
                ggplot2::geom_boxplot(width = 0.1, fill = "white", show.legend = FALSE) +
                ggplot2::scale_color_brewer(palette = "Dark2")
            Plots[[length(Plots) + 1]] <- Plot
        }

    }
    return(Plots)
}


#' @title Bar plots
#'
#' @description Return bar plots to compare groups
#'
#' @param data A dataframe with numeric columns
#' @param belongmatrix A membership matrix
#' @param what Can be "mean" (default) or "median"
#' @param ncol An integer indicating the number of columns for the bar plot
#' @return a barplot created with ggplot2
#' @export
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' barPlots(dataset, result$Belongings)
barPlots <- function(data,belongmatrix, ncol = 3, what = "mean"){
    datasummary <- summarizeClusters(data, belongmatrix)
    if (what == "mean"){
        values <- lapply(datasummary, function(i){as.numeric(i[8,])})
    } else if (what == "median"){
        values <- lapply(datasummary, function(i){as.numeric(i[4,])})
    }else{
        warning('The parameter what is invalid, using what = "mean"')
        values <- lapply(datasummary, function(i){as.numeric(i[8,])})
    }

    values <- data.frame(do.call(rbind, values))
    names(values) <- colnames(datasummary$Cluster_1)
    values$Cluster <- rownames(values)
    values <- reshape2::melt(values, id.vars = "Cluster")


    faceplot <- ggplot2::ggplot(values) +
        ggplot2::geom_bar(ggplot2::aes_string(x = "Cluster", weight = "value", fill = "Cluster"), width = 0.7) +
        ggplot2::theme(panel.background = ggplot2::element_blank(),
              panel.grid = ggplot2::element_blank(),
              axis.text.x = ggplot2::element_blank(),
              axis.ticks.x = ggplot2::element_blank(),
              axis.title.y = ggplot2::element_blank()
        ) +
        ggplot2::facet_wrap(ggplot2::vars(values$variable), ncol=ncol, scales="free_y")
    return(faceplot)

}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Functions to select parameters #####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Worker function
#'
#' @description Worker function for select_parameters and select_parameters.mc
#'
#' @param algo A string indicating which method to use (FCM, GFCM, SFCM, SGFCM)
#' @param parameters A dataframe of parameters with columns k,m and alpha
#' @param data A dataframe with numeric columns
#' @param nblistw A list.w object describing the neighbours typically produced
#'   by the spdep package
#' @param standardize A boolean to specify if the variable must be centered and
#'   reduce (default = True)
#' @param spconsist A boolean indicating if the spatial consistency must be
#' calculated
#' @param classidx A boolean indicating if the quality of classification
#' indices must be calculated
#' @param maxiter An integer for the maximum number of iteration
#' @param tol The tolerance criterion used in the evaluateMatrices function for
#'   convergence assessment
#' @param seed An integer used for random number generation. It ensures that the
#' start centers will be the same if the same integer is selected.
#' @param verbose A boolean indicating if a progressbar should be displayed
#' @return a DataFrame containing for each combinations of parameters several clustering
#' quality indexes.
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @keywords internal
#' @examples
#' #No example provided, this is an internal function
eval_parameters <- function(algo, parameters, data, nblistw = NULL, standardize = TRUE ,spconsist = FALSE, classidx = TRUE, tol, maxiter, seed=NULL, verbose = TRUE){
    if(algo == "FCM"){
        exefun <- function(data,x, lw){
            return(CMeans(data, x$k, x$m, maxiter = maxiter, tol = tol, standardize = standardize, verbose = FALSE, seed = seed))
        }
    }else if(algo == "GFCM"){
        exefun <- function(data,x, lw){
            return(GCMeans(data, x$k, x$m, x$beta, maxiter = maxiter, tol = tol, standardize = standardize, verbose = FALSE, seed = seed))
        }
    }else if(algo == "SFCM"){
        exefun <- function(data,x, lw){
            return(SFCMeans(data, lw, x$k, x$m, x$alpha, x$lag_method, maxiter = maxiter, tol = tol, standardize = standardize, verbose = FALSE, seed = seed))
        }
    }else if(algo == "SGFCM"){
        exefun <- function(data,x, lw){
            return(SGFCMeans(data, lw, x$k, x$m, x$alpha, x$beta, x$lag_method, maxiter = maxiter, tol = tol, standardize = standardize, verbose = FALSE, seed = seed))
        }
    }else{
        stop("The algo selected must be one in FCM, GFCM, SFCM, SGFCM")
    }

    if(verbose){
        pb <- txtProgressBar(min = 0, max = nrow(parameters), style = 3)
    }
    cnt <- 1
    allIndices <- lapply(1:nrow(parameters), function(i){
        row <- parameters[i,]
        if(verbose){
            setTxtProgressBar(pb, cnt)
        }
        cnt <<- cnt+1
        templistw <- nblistw[[row$listsw]]
        result <- exefun(data,row,templistw)
        #calculating the quality indexes
        indices <- list()
        if(classidx){
            indices <- calcqualityIndexes(result$Data,result$Belongings,as.numeric(row[[2]]))
        }
        if(spconsist){
            #calculating spatial diag
            consist <- spConsistency(result$Belongings, templistw, nrep = 30)
            indices$spConsistency <- consist$Mean
            indices$spConsistency_05 <- consist$prt05
            indices$spConsistency_95 <- consist$prt95
        }

        return(unlist(indices))
    })
    dfIndices <- data.frame(do.call(rbind,allIndices))
    dfIndices$k <- parameters$k
    dfIndices$m <- parameters$m
    if(algo %in% c("SFCM","SGFCM")){
        dfIndices$alpha <- parameters$alpha
        dfIndices$listw <- parameters$listsw
        dfIndices$lag_method <- parameters$lag_method
    }
    if(algo %in% c("GFCM","SGFCM")){
        dfIndices$beta <- parameters$beta
    }
    return(dfIndices)
}


#' @title Select parameters for a clustering algorithm
#'
#' @description Function to select the parameters for a clustering algorithm.
#'
#' @param algo A string indicating which method to use (FCM, GFCM, SFCM, SGFCM)
#' @param data A dataframe with numeric columns
#' @param k A sequence of values for k to test (>=2)
#' @param m A sequence of values for m to test
#' @param alpha A sequence of values for alpha to test (NULL if not required)
#' @param beta A sequence of values for beta to test (NULL if not required)
#' @param nblistw A list of list.w objects describing the neighbours typically
#'  produced by the spdep package (NULL if not required)
#' @param lag_method A string indicating if a classical lag must be used
#' ("mean") or if a weighted median must be used ("median"). Both can be
#' tested by specifying a vector : c("mean","median")
#' @param standardize A boolean to specify if the variable must be centered and
#'   reduce (default = True)
#' @param spconsist A boolean indicating if the spatial consistency must be
#' calculated
#' @param classidx A boolean indicating if the quality of classification
#' indices must be calculated
#' @param maxiter An integer for the maximum number of iteration
#' @param tol The tolerance criterion used in the evaluateMatrices function for
#'   convergence assessment
#' @param seed An integer used for random number generation. It ensures that the
#' start centers will be the same if the same integer is selected.
#' @param verbose A boolean indicating if a progressbar should be displayed
#' @return A dataframe with indicators assessing the quality of classifications
#' @export
#' @examples
#' \donttest{
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' #set spconsist to TRUE to calculate the spatial consistency indicator
#' #FALSE here to reduce the time during package check
#' values <- select_parameters(algo = "SFCM", dataset, k = 5, m = seq(2,3,0.1),
#'     alpha = seq(0,2,0.1), nblistw = Wqueen, spconsist=FALSE)
#' }
select_parameters <- function(algo,data,k,m,alpha = NA, beta = NA, nblistw=NULL, lag_method="mean", spconsist = TRUE, classidx = TRUE, standardize = TRUE, maxiter = 500, tol = 0.01, seed=NULL, verbose = TRUE){

    if(spconsist==FALSE & classidx==FALSE){
        stop("one of spconsist and classidx must be TRUE")
    }

    if(class(nblistw)[[1]]!="list"){
        nblistw <- list(nblistw)
    }
    allcombinaisons <- expand.grid(k=k,m=m,alpha=alpha,beta = beta,listsw=1:length(nblistw),lag_method=lag_method)

    print(paste("number of combinaisons to estimate : ",nrow(allcombinaisons)))
    dfIndices <- eval_parameters(algo,allcombinaisons, data, nblistw, standardize,
        spconsist, classidx, tol, maxiter, seed, verbose)
}


#' @rdname select_parameters
#' @examples
#' \donttest{
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' #set spconsist to TRUE to calculate the spatial consistency indicator
#' #FALSE here to reduce the time during package check
#' values <- selectParameters(algo = "SFCM", dataset, k = 5, m = seq(2,3,0.1),
#'     alpha = seq(0,2,0.1), nblistw = Wqueen, spconsist=FALSE)
#' }
#' @export
selectParameters <- select_parameters


#' @title Select parameters for clustering algorithm (multicore)
#'
#' @description Function to select the parameters for a clustering algorithm.
#' This version of the function allows to use a plan defined with the package
#' future to reduce calculation time.
#'
#' @param algo A string indicating which method to use (FCM, GFCM, SFCM, SGFCM)
#' @param data A dataframe with numeric columns
#' @param k A sequence of values for k to test (>=2)
#' @param m A sequence of values for m to test
#' @param alpha A sequence of values for alpha to test (NULL if not required)
#' @param beta A sequence of values for beta to test (NULL if not required)
#' @param nblistw A list of list.w objects describing the neighbours typically
#'  produced by the spdep package (NULL if not required)
#' @param lag_method A string indicating if a classical lag must be used
#' ("mean") or if a weighted median must be used ("median"). Both can be
#' tested by specifying a vector : c("mean","median")
#' @param spconsist A boolean indicating if the spatial consistency must be
#' calculated
#' @param classidx A boolean indicating if the quality of classification
#' indices must be calculated
#' @param standardize A boolean to specify if the variable must be centered and
#'   reduce (default = True)
#' @param maxiter An integer for the maximum number of iteration
#' @param tol The tolerance criterion used in the evaluateMatrices function for
#'   convergence assessment
#' @param seed An integer used for random number generation. It ensures that the
#' start centers will be the same if the same integer is selected.
#' @param chunk_size The size of a chunk used for multiprocessing. Default is 100.
#' @param verbose A boolean indicating if a progressbar should be displayed
#' @return A dataframe with indicators assessing the quality of classifications
#' @export
#' @importFrom utils setTxtProgressBar txtProgressBar capture.output
#' @examples
#' \donttest{
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' future::plan(future::multiprocess(workers=2))
#' #set spconsist to TRUE to calculate the spatial consistency indicator
#' #FALSE here to reduce the time during package check
#' values <- select_parameters.mc("SFCM", dataset, k = 5, m = seq(1,2.5,0.1),
#'     alpha = seq(0,2,0.1), nblistw = Wqueen, spconsist=FALSE)
#' \dontshow{
#'    ## R CMD check: make sure any open connections are closed afterward
#'    if (!inherits(future::plan(), "sequential")) future::plan(future::sequential)
#' }
#'}
select_parameters.mc <- function(algo,data,k,m, alpha = NA, beta = NA, nblistw = NULL, lag_method="mean",  spconsist = TRUE, classidx = TRUE, standardize = TRUE, maxiter = 500, tol = 0.01, seed = NULL, chunk_size=100, verbose = FALSE){

    if(spconsist==FALSE & classidx==FALSE){
        stop("one of spconsist and classidx must be TRUE")
    }

    if(class(nblistw)[[1]]!="list"){
        nblistw <- list(nblistw)
    }
    if(is.null(seed)){
        seed <- FALSE
    }
    allcombinaisons <- expand.grid(k=k,m=m,alpha=alpha, beta = beta, listsw=1:length(nblistw), lag_method=lag_method)
    if (verbose){
        print(paste("number of combinaisons to estimate : ",nrow(allcombinaisons)))
    }
    chunks <- split(1:nrow(allcombinaisons), rep(1:ceiling(nrow(allcombinaisons) / chunk_size),
                                       each = chunk_size, length.out = nrow(allcombinaisons)))
    chunks <- lapply(chunks,function(x){return(allcombinaisons[x,])})
    # step2 : starting the function
    iseq <- 1:length(chunks)
    if(verbose){
        progressr::with_progress({
            p <- progressr::progressor(along = iseq)
            values <- future.apply::future_lapply(iseq, function(i) {
                sprintf(algo)
                parameters <- chunks[[i]]
                indices <- eval_parameters(algo, parameters, data, nblistw, standardize,
                                                             spconsist, classidx, tol, maxiter)
                p(sprintf("i=%g", i))
                return(indices)
            }, future.seed = seed)
        })
    }else{
        values <- future.apply::future_lapply(iseq, function(i) {
            parameters <- chunks[[i]]
            indices <- eval_parameters(algo, parameters, data, nblistw, standardize,
                                                    spconsist, classidx, tol, maxiter,
                                       verbose = FALSE)
            return(indices)
        },future.seed = seed)
    }

    dfIndices <- do.call(rbind,values)
    return(dfIndices)
}

#' @rdname select_parameters.mc
#' @examples
#' \donttest{
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' future::plan(future::multiprocess(workers=2))
#' #set spconsist to TRUE to calculate the spatial consistency indicator
#' #FALSE here to reduce the time during package check
#' values <- select_parameters.mc("SFCM", dataset, k = 5, m = seq(1,2.5,0.1),
#'     alpha = seq(0,2,0.1), nblistw = Wqueen, spconsist=FALSE)
#' \dontshow{
#'    ## R CMD check: make sure any open connections are closed afterward
#'    if (!inherits(future::plan(), "sequential")) future::plan(future::sequential)
#' }
#'}
#' @export
selectParameters.mc <- select_parameters.mc


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Functions to recalculate spatial weights #####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Semantic adjusted spatial weights
#'
#' @description Function to adjust the spatial weights so that they represent semantic
#' distances between neighbours
#'
#' @param data A dataframe with numeric columns
#' @param listw A nb object from spdep
#' @param style A letter indicating the weighting scheme (see spdep doc)
#'
#' @return A listw object (spdep like)
#' @export
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' Wqueen2 <- adjustSpatialWeights(dataset,queen,style="C")
adjustSpatialWeights <- function(data,listw,style){
    new_weights <- lapply(1:nrow(data),function(i){
        row <- data[i,]
        neighbours <- data[listw[[i]],]
        dists <- 1/calcEuclideanDistance(neighbours,row)
        weights <- dists / sum(dists)
        return(weights)
    })
    new_listw <- spdep::nb2listw(listw,glist=new_weights,style = style)
    return(new_listw)
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### prediction functions #####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @title Predict matrix membership for new observations
#'
#' @description Function predict the membership matrix of a new set of observations
#'
#' @param object An FCMres object as produced by one of the following functions: CMeans,
#' GCMeans, SFCMeans, SGFCMeans
#' @param new_data A DataFrame with the new observations
#' @param listw A nb object from spdep (for the new data)
#' @param standardize  A boolean to specify if the variable must be centered and
#'   reduced (default = True)
#' @param ... not used
#'
#' @return A numeric matrix with the membership values for each new observation
#' @export
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#'
#' # rescaling all the variables used in the analysis
#' for (field in AnalysisFields) {
#'     LyonIris@data[[field]] <- scale(LyonIris@data[[field]])
#' }
#'
#' # doing the initial clustering
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SGFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, beta = 0.5, standardize = FALSE)
#'
#' # using a subset of the original dataframe as "new data"
#' new_data <- LyonIris[c(1, 27, 36, 44, 73),]
#' new_dataset <- new_data@data[AnalysisFields]
#' new_nb <- spdep::poly2nb(new_data,queen=TRUE)
#' new_Wqueen <- spdep::nb2listw(new_nb,style="W")
#'
#' # doing the prediction
#' predictions <- predict_membership(result, new_dataset, new_Wqueen, standardize = FALSE)
predict_membership <- function(object, new_data, listw = NULL, standardize = TRUE, ...){

    results <- object
    ## standardize the data if required
    if (standardize){
        for (i in 1:ncol(new_data)) {
            new_data[, i] <- scale(new_data[, i])
        }
    }

    ## calculating the lagged dataset if required
    if(results$algo %in% c("SFCM", "SGFCM")){
        if(is.null(listw)){
            stop("With a spatial clustering method like SFCM or SGFCM, the spatial matrix listw
                 associated with the new dataset must be provided")
        }
        wdata <- calcLaggedData(new_data, listw, results$lag_method)
    }

    ## selecting the appropriate function for prediction
    if(results$algo == "FCM"){
        pred_values <- calcBelongMatrix(results$Centers, new_data, results$m)
    }else if(results$algo == "GFCM"){
        pred_values <- calcFGCMBelongMatrix(results$Centers, new_data,
                                            results$m, results$beta)
    }else if(results$algo == "SFCM"){
        pred_values <- calcSFCMBelongMatrix(results$Centers, new_data,
                                            wdata,results$m, results$alpha)
    }else if(results$algo == "SGFCM"){
        pred_values <- calcSFGCMBelongMatrix(results$Centers, new_data,
                                            wdata,results$m, results$alpha, results$beta)
    }

    return(pred_values)
}


#' @title Predict method for FCMres object
#'
#' @description Function predict the membership matrix of a new set of observations
#'
#' @param object An FCMres object as produced by one of the following functions: CMeans,
#' GCMeans, SFCMeans, SGFCMeans
#' @param new_data A DataFrame with the new observations
#' @param listw A nb object from spdep (for the new data)
#' @param standardize  A boolean to specify if the variable must be centered and
#'   reduced (default = True)
#' @param ... not used
#'
#' @return A numeric matrix with the membership values for each new observation
#' @export
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#'
#' # rescaling all the variables used in the analysis
#' for (field in AnalysisFields) {
#'     LyonIris@data[[field]] <- scale(LyonIris@data[[field]])
#' }
#'
#' # doing the initial clustering
#' dataset <- LyonIris@data[AnalysisFields]
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SGFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, beta = 0.5, standardize = FALSE)
#'
#' # using a subset of the original dataframe as "new data"
#' new_data <- LyonIris[c(1, 27, 36, 44, 73),]
#' new_dataset <- new_data@data[AnalysisFields]
#' new_nb <- spdep::poly2nb(new_data,queen=TRUE)
#' new_Wqueen <- spdep::nb2listw(new_nb,style="W")
#'
#' # doing the prediction
#' predictions <- predict(result, new_dataset, new_Wqueen, standardize = FALSE)
predict.FCMres <- predict_membership


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Utilitary functions #####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Convert categories to membership matrix
#'
#' @description Function to Convert categories to membership matrix (binary matrix)
#'
#' @param categories A vector with the categories of each observation
#'
#' @return A binary matrix
#' @export
cat_to_belongings <- function(categories){
    cats <- unique(categories)
    cols <- lapply(cats, function(i){
        return (ifelse(categories == i, 1, 0))
    })
    mat <- do.call(cbind, cols)
    return(mat)
}

#' @rdname cat_to_belongings
#' @export
catToBelongings <- cat_to_belongings


