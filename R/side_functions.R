# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Fonctions de visualisation #####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Mapping the clusters
#'
#' @description Build some maps to visualize the results of the clustering
#'
#' @param geodata An object of class features collection from sf /
#' ordered like the original data used for the clustering. Can be Null if object is
#' a FCMres and has been created with rasters.
#' @param object  A FCMres object, typically obtained from functions CMeans,
#'   GCMeans, SFCMeans, SGFCMeans. Can also be a simple membership matrix.
#' @param undecided A float between 0 and 1 giving the minimum value that an
#'   observation must get in the membership matrix to not be considered as
#'   uncertain (default = NULL)
#' @return A named list with :
#' \itemize{
#'         \item ProbaMaps : a list of tmap maps showing for each group the
#'         probability of the observations to belong to that group
#'         \item ClusterMap : a tmap map showing the most likely group for
#'          observation
#' }
#' @export
#' @examples
#' \dontrun{
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- sf::st_drop_geometry(LyonIris[AnalysisFields])
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' MyMaps <- mapClusters(LyonIris, result$Belongings)
#' }
mapClusters <- function(geodata = NULL, object, undecided = NULL) {# nocov start

    #cls <- class(object)[[1]]
    isRaster <- FALSE
    if(inherits(object,"FCMres")){
        belongmatrix <- object$Belongings
        isRaster <- object$isRaster
    }else if(inherits(object,"matrix")){
        belongmatrix <- object
    }else{
        stop("object must be one of matrix of FCMres")
    }

    if(isRaster){
        return(mapRasters(object, undecided))
    }else{
        geodata$OID <- as.character(1:nrow(geodata))
        sf_type <- sf::st_geometry_type(geodata,by_geometry = FALSE)
        if(sf_type %in% c("POLYGON","MULTIPOLYGON")){
          geom_type <- "polygon"
        }else if(sf_type %in% c("POINT","MULTIPOINT")){
          geom_type <- "point"
        }else if(sf_type %in% c("LINESTRING","MULTILINESTRING")){
          geom_type <- "line"
        }else{
            stop("The object passed in geodata argument is not supported.
              Supported classes are feature collections from sf with geometry type in
              POLYGON, MULTIPOLYGON, POINT, MULTIPOINT, LINESTRING, MULTILINESTRING")
        }
        return(mapThis(geodata, belongmatrix, undecided, geom_type))
        # if(sf::st_geometry_type(geodata,by_geometry = FALSE) %in% c("POLYGON","MULTIPOLYGON")){
        #     return(mapPolygons(geodata, belongmatrix, undecided))
        # }else if(sf::st_geometry_type(geodata,by_geometry = FALSE) %in% c("POINT","MULTIPOINT")){
        #     return(mapPoints(geodata, belongmatrix, undecided))
        #   else if(sf::st_geometry_type(geodata,by_geometry = FALSE) %in% c("LINESTRING","MULTILINESTRING")){
        #     return(mapLines(geodata, belongmatrix, undecided))
        # }else {
        #     stop("The object passed in geodata argument is not supported.
        #       Supported classes are feature collections from sf with geometry type in
        #       POLYGON, MULTIPOLYGON, POINT, MULTIPOINT, LINESTRING, MULTILINESTRING")
        # }
    }


}

#' @title Mapping the clusters (rasters)
#'
#' @description Internal function to realize maps based on rasters
#'
#' @param object A FCMres object
#' @param undecided A float between 0 and 1 giving the minimum value that an
#'   observation must get in the membership matrix to not be considered as
#'   uncertain (default = NULL)
#' @return A named list with :
#' \itemize{
#'         \item ProbaMaps : a list of ggplot maps showing for each group the
#'         probability of the observations to belong to that group
#'         \item ClusterMap : a ggplot map showing the most likely group for each observation
#' }
#' @importFrom methods as
#' @keywords internal
#' @examples
#' #No example provided, this is an internal function, use the general wrapper function mapClusters
mapRasters  <- function(object, undecided){

    # realisation des cartes de probabilites
    ProbaPlots <- lapply(1:ncol(object$Belongings), function(i) {
        rast <- object$rasters[[i]]
        #spdf <- as(rast, "SpatialPixelsDataFrame")
        #df <- as.data.frame(spdf)
        df <- terra::as.data.frame(rast, xy = TRUE)
        colnames(df) <- c("x", "y", "value")

        Plot <- ggplot2::ggplot(df) +
            ggplot2::geom_raster(ggplot2::aes_string(x="x", y="y", fill="value"), alpha=0.8) +
            ggplot2::scale_fill_gradient(low = "white", high = "blue") +
            ggplot2::coord_fixed(ratio = 1)+
            ggplot2::theme(
                axis.title = ggplot2::element_blank(),
                axis.text = ggplot2::element_blank(),
                axis.ticks = ggplot2::element_blank()
            ) +
            ggplot2::labs(title = paste("membership values for group ",i,sep=""))
        return(Plot)
    })

    #finding the uncertain pixels
    if(is.null(undecided)){
        undecided <- 0
    }
    #finding for each pixel the max probability
    uncertain_vec <- undecidedUnits(object$Belongings, tol = undecided, out = "numeric")
    rast <- object$rasters[[1]]
    vec <- rep(NA, times = terra::ncell(rast))
    vec[object$missing] <- uncertain_vec
    terra::values(rast) <- vec

    #spdf <- as(rast, "SpatialPixelsDataFrame")
    #df <- as.data.frame(spdf)
    df <- terra::as.data.frame(rast, xy = TRUE)
    colnames(df) <- c("x", "y","value")
    df$value <- as.character(df$value)
    df$value <- ifelse(df$value == "-1", "undecided", paste("group", df$value, sep = " "))

    colors <- RColorBrewer::brewer.pal(ncol(object$Belongings),"Set3")
    if("undecided" %in% df$values){
        colors <- c(colors, "black")
    }

    ClusterMap <- ggplot2::ggplot(df) +
        ggplot2::geom_raster(ggplot2::aes_string(x = "x", y = "y", fill = "value")) +
        ggplot2::scale_fill_discrete(type = colors) +
        ggplot2::coord_fixed(ratio = 1)+
        ggplot2::theme(
            legend.title = ggplot2::element_blank(),
            axis.title = ggplot2::element_blank(),
            axis.text = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank()
        )
    return(list(ProbaMaps = ProbaPlots, ClusterPlot = ClusterMap))
}


#' @title Mapping the clusters
#'
#' @description Internal function to realize maps
#'
#' @param geodata feature collections ordered like the original data used
#'   for the clustering
#' @param belongmatrix The membership matrix obtained at the end of the algorithm
#' @param undecided A float between 0 and 1 giving the minimum value that an
#'   observation must get in the membership matrix to not be considered as
#'   uncertain (default = NULL)
#' @param geom_type A string indicating the type of geometry (polygon, string or point)
#' @return A named list with :
#' \itemize{
#'         \item ProbaMaps : a list of ggplot maps showing for each group the
#'         probability of the observations to belong to that group
#'         \item ClusterMap : a ggplot map showing the most likely group for each observation
#' }
#' @importFrom tmap tm_shape tm_fill tm_borders tm_dots tm_layout tm_lines
#' @keywords internal
#' @examples
#' #No example provided, this is an internal function, use the general wrapper function mapClusters
mapThis <- function(geodata, belongmatrix, undecided = NULL, geom_type = "polygon"){

    belongmatrix <- as.data.frame(belongmatrix)
    names(belongmatrix) <- paste0("cluster",1:ncol(belongmatrix))
    geodata <- cbind(geodata,belongmatrix)

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
      if(geom_type == "polygon"){
        this_map <- tm_shape(geodata) +
          tm_fill(col = Name, palette = "Blues", style = "cont") +
          tm_borders("black") +
          tm_layout(legend.outside = TRUE, frame = FALSE)
      }else if(geom_type == "line"){
        this_map <- tm_shape(geodata) +
          tm_lines(col = Name, palette = "Blues", style = "cont") +
          tm_layout(legend.outside = TRUE, frame = FALSE)
      }else if(geom_type == "point"){
        this_map <- tm_shape(geodata) +
          tm_dots(col = Name, palette = "Blues", style = "cont", size = 0.1) +
          tm_layout(legend.outside = TRUE, frame = FALSE)
      }

        return(this_map)
    })
    # realisation de la carte des groupes
    if(geom_type == "polygon"){
      ClusterMap <- tm_shape(geodata) +
        tm_fill("Cluster", palette = "Set1") +
        tm_borders("black")+
        tm_layout(legend.outside = TRUE, frame = FALSE)
    }else if (geom_type == "line"){
      ClusterMap <- tm_shape(geodata) +
        tm_lines("Cluster", palette = "Set1") +
        tm_layout(legend.outside = TRUE, frame = FALSE)
    }else if(geom_type == "point"){
      ClusterMap <- tm_shape(geodata) +
        tm_dots("Cluster", palette = "Set1", size = 0.1) +
        tm_layout(legend.outside = TRUE, frame = FALSE)
    }


    if(is.null(undecided) == FALSE){
      if(geom_type == "polygon"){
        ClusterMap <- ClusterMap +
          tm_shape(subset(geodata, geodata$Undecided == "Undecided")) +
          tm_fill("white", alpha = 0.7)+
          tm_borders("black") +
          tm_layout(legend.outside = TRUE, frame = FALSE)
      }else if (geom_type == "line"){
        ClusterMap <- ClusterMap +
          tm_shape(subset(geodata, geodata$Undecided == "Undecided")) +
          tm_lines("white", alpha = 0.7)+
          tm_layout(legend.outside = TRUE, frame = FALSE)
      }else if (geom_type == "point"){
        ClusterMap <- ClusterMap +
          tm_shape(subset(geodata, geodata$Undecided == "Undecided")) +
          tm_dots("white", alpha = 0.7, size = 0.1)+
          tm_layout(legend.outside = TRUE, frame = FALSE)
        }

    }

    return(list(ProbaMaps = ProbaPlots, ClusterPlot = ClusterMap))
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
#' dataset <- sf::st_drop_geometry(LyonIris[AnalysisFields])
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
#' \dontrun{
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- sf::st_drop_geometrie(LyonIris[AnalysisFields])
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' violinPlots(dataset, result$Groups)
#' }
violinPlots <- function(data,groups){
    data <- as.data.frame(data)
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
#' dataset <- sf::st_drop_geometry(LyonIris[AnalysisFields])
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



#force_sp_sample <- function(spatial, n, type){
#    result <- tryCatch(
#        {sp::spsample(spatial,n,type,iter = 10)},
#        error = function(e){
#            tryCatch(
#                {sp::spsample(spatial,n,type = "regular",iter = 10)},
#                error = function(ee){
#                    stop(sprintf('Error when drawing random points: ',ee$message))
#                }
#                )
#        }
#    )
#    return(result)
#}


#' @title Uncertainty map
#'
#' @description Return a map to visualize membership matrix
#'
#' @details This function maps the membership matrix by plotting
#' random points in polygons, along lines or around points representing the
#' original observations. Each cluster is associated with a color and each
#' random point has a probability to be of that color equal to the membership
#' value of the feature it belongs itself. Thus, it is possible to
#' visualize regions with uncertainty and to identify the strongest clusters.
#'
#' @param geodata An object of class feature collection from sf ordered
#' like the original data used for the clustering.
#' @param belongmatrix A membership matrix
#' @param njit The number of points to map on each feature.
#' @param radius When mapping points, the radius indicates how far random
#' points will be plotted around the original features.
#' @param colors A vector of colors to use for the groups.
#' @param pt_size A float giving the size of the random points on the final
#' map (default is 0.05)
#' @return a map created with tmap
#' @export
#' @examples
#' \dontrun{
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#'   "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- sf::st_drop_geometry(LyonIris[AnalysisFields])
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
#' uncertaintyMap(LyonIris, result$Belongings)
#' }
uncertaintyMap <- function(geodata, belongmatrix, njit = 150, radius = NULL, colors = NULL, pt_size = 0.05){

    geodata$tmpOID <- 1:nrow(geodata)
    groups <- 1:ncol(belongmatrix)
    #cls <- class(geodata)[[1]]

    geom_type <- sf::st_geometry_type(geodata, by_geometry = FALSE)

    if(geom_type %in% c("POLYGON","MULTIPOLYGON")){
        geodata$area <- sf::st_area(geodata)

    }else if(geom_type %in% c("LINESTRING", "MULTILINESTRING")){
        geodata$area <- sf::st_length(geodata)
    }else if(geom_type %in% c("POINT","MULTIPOINT")){
        geodata$area <- 1
        coords <- sf::st_coordinates(geodata)
        geodata$X <- coords[,1]
        geodata$Y <- coords[,2]
        if(is.null(radius)){
            stop("When mapping points, the parameter radius can not be NULL")
        }
    }else{
        stop("geodata must be a feature collection (POINT, MULTIPOINT, LINESTRING, MULTILINESTRING, POLYGON or MULTIPOLYGON)")
    }
    maxA <- max(geodata$area)
    rt <- as.numeric(njit / maxA)


    if((geom_type %in% c("POINT","MULTIPOINT"))==FALSE){
      n_obs <- as.integer(as.numeric(geodata$area) * rt)
      jitted <- sf::st_sample(geodata, n_obs)
      coords <- sf::st_coordinates(jitted)
    }else{
      jitted <- sf::st_sample(
        sf::st_buffer(geodata, dist = radius),
        njit)
      coords <- sf::st_coordinates(jitted)
      n_obs <- rep(njit, nrow(geodata))
    }

    clusts <- lapply(1:nrow(geodata), function(i){
      probs <- belongmatrix[i,]
      clusters <- sample(groups, size = n_obs[[i]], prob = probs, replace = TRUE)
      return(clusters)
    })
    clusts <- do.call(c, clusts)

    all_pts <- sf::st_sf(data.frame(
      "oid" = 1:length(jitted),
      "cluster" = as.character(clusts)
    ),jitted)

    if(is.null(colors)){
        colors <- c("#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD","#8C564B",
                    "#E377C2","#7F7F7F","#BCBD22","#17BECF","#AEC7E8","#FFBB78",
                    "#98DF8A","#FF9896","#C5B0D5","#C49C94","#F7B6D2","#C7C7C7",
                    "#DBDB8D","#9EDAE5")[1:ncol(belongmatrix)]
    }

    names(colors) <- as.character(groups)

    uncertain_map <- tm_shape(all_pts) +
      tm_dots("cluster", palette = colors, size = 0.01) +
      tm_layout(frame = FALSE, legend.outside = TRUE)

    if(geom_type %in% c("POLYGON","MULTIPOLYGON")){
      uncertain_map <- uncertain_map +
        tm_shape(geodata) +
        tm_borders("black")
    }

    return(uncertain_map)

} # nocov end


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
#' @template nblistw-arg
#' @template window-arg
#' @param standardize A boolean to specify if the variable must be centered and
#'   reduce (default = True)
#' @param spconsist A boolean indicating if the spatial consistency must be
#' calculated
#' @param classidx A boolean indicating if the quality of classification
#' indices must be calculated
#' @param nrep An integer indicating the number of permutation to do to simulate
#'   the random distribution of the spatial inconsistency. Only used if spconsist
#'   is TRUE.
#' @param indices A character vector with the names of the indices to calculate, to
#' evaluate clustering quality. default is :c("Silhouette.index", "Partition.entropy",
#' "Partition.coeff", "XieBeni.index", "FukuyamaSugeno.index", "Explained.inertia").
#' Other available indices are : "DaviesBoulin.index", "CalinskiHarabasz.index",
#' "GD43.index", "GD53.index" and "Negentropy.index".
#' @param maxiter An integer for the maximum number of iteration
#' @param tol The tolerance criterion used in the evaluateMatrices function for
#'   convergence assessment
#' @param seed An integer used for random number generation. It ensures that the
#' start centers will be the same if the same integer is selected.
#' @param init A string indicating how the initial centers must be selected. "random"
#' indicates that random observations are used as centers. "kpp" use a distance based method
#' resulting in more dispersed centers at the beginning. Both of them are heuristic.
#' @param verbose A boolean indicating if a progressbar should be displayed
#' @param wrapped A boolean indicating if the data passed is wrapped or not (see wrap function of terra)
#' @return a DataFrame containing for each combinations of parameters several clustering
#' quality indexes.
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @keywords internal
#' @examples
#' #No example provided, this is an internal function
eval_parameters <- function(algo, parameters, data, nblistw = NULL, window = NULL, standardize = TRUE,
                            robust = FALSE, noise_cluster = FALSE,
                            spconsist = FALSE, classidx = TRUE, nrep = 30, indices = NULL,
                            tol, maxiter, seed = NULL, init = "random", verbose = TRUE,
                            wrapped = FALSE){
    if(wrapped){
      data <- lapply(data, terra::unwrap)
    }

    if(algo == "FCM"){
        exefun <- function(data,x, ...){
            return(CMeans(data, x$k, x$m, maxiter = maxiter, tol = tol, standardize = standardize,
                          robust = robust, noise_cluster = noise_cluster, delta = x$delta,
                          verbose = FALSE, seed = seed, init = init))
        }
    }else if(algo == "GFCM"){
        exefun <- function(data,x,... ){
            return(GCMeans(data, x$k, x$m, x$beta, maxiter = maxiter, tol = tol, standardize = standardize,
                           robust = robust, noise_cluster = noise_cluster, delta = x$delta,
                           verbose = FALSE, seed = seed, init = init))
        }
    }else if(algo == "SFCM"){
        exefun <- function(data,x,...){
            dots <- list(...)
            return(SFCMeans(data, dots$lw, x$k, x$m, x$alpha, x$lag_method, window = dots$wd,
                            maxiter = maxiter, tol = tol, standardize = standardize,
                            robust = robust, noise_cluster = noise_cluster, delta = x$delta,
                            verbose = FALSE, seed = seed, init = init))
        }
    }else if(algo == "SGFCM"){
        exefun <- function(data,x,...){
            dots <- list(...)
            return(SGFCMeans(data, dots$lw, x$k, x$m, x$alpha, x$beta, x$lag_method, window = dots$wd,
                             maxiter = maxiter, tol = tol, standardize = standardize,
                             robust = robust, noise_cluster = noise_cluster, delta = x$delta,
                             verbose = FALSE, seed = seed, init = init))
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
        cnt <<- cnt+1
        templistw <- nblistw[[row$listsw]]
        tempwindow <- window[[row$window]]
        result <- exefun(data, row, lw = templistw, wd = tempwindow)
        #calculating the quality indexes
        indices_values <- list()
        if(classidx){
            indices_values <- calcqualityIndexes(result$Data,result$Belongings,as.numeric(row[[2]]),
                                          indices = indices)
        }
        if(spconsist){
            #calculating spatial diag
            consist <- spConsistency(result, nrep = nrep)
            indices_values$spConsistency <- consist$Mean
            indices_values$spConsistency_05 <- consist$prt05
            indices_values$spConsistency_95 <- consist$prt95
        }
        if(verbose){
            setTxtProgressBar(pb, cnt)
        }

        return(unlist(indices_values))
    })
    dfIndices <- data.frame(do.call(rbind,allIndices))
    dfIndices$k <- parameters$k
    dfIndices$m <- parameters$m
    if(algo %in% c("SFCM","SGFCM")){
        dfIndices$alpha <- parameters$alpha
        dfIndices$listw <- parameters$listsw
        dfIndices$window <- parameters$window
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
#' @param data A dataframe with numeric columns or a list of rasters.
#' @param k A sequence of values for k to test (>=2)
#' @param m A sequence of values for m to test
#' @param alpha A sequence of values for alpha to test (NULL if not required)
#' @param beta A sequence of values for beta to test (NULL if not required)
#' @param nblistw A list of list.w objects describing the neighbours typically
#'  produced by the spdep package (NULL if not required)
#' @param lag_method A string indicating if a classical lag must be used
#' ("mean") or if a weighted median must be used ("median"). Both can be
#' tested by specifying a vector : c("mean","median"). When working with rasters,
#' the string must be parsable to a function like mean, min, max, sum, etc. and will
#' be applied to all the pixels values in the window designated by the parameter window
#' and weighted according to the values of this matrix.
#' @param window A list of windows to use to calculate neighbouring values if
#' rasters are used.
#' @param robust A boolean indicating if the "robust" version of the algorithm must be used (see details)
#' @param noise_cluster A boolean indicatong if a noise cluster must be added to the solution (see details)
#' @param delta A float giving the distance of the noise cluster to each observation
#' @param standardize A boolean to specify if the variable must be centered and
#'   reduce (default = True)
#' @param spconsist A boolean indicating if the spatial consistency must be
#' calculated
#' @param classidx A boolean indicating if the quality of classification
#' indices must be calculated
#' @param nrep An integer indicating the number of permutation to do to simulate
#'   the random distribution of the spatial inconsistency. Only used if spconsist
#'   is TRUE.
#' @param indices A character vector with the names of the indices to calculate, to
#' evaluate clustering quality. default is :c("Silhouette.index", "Partition.entropy",
#' "Partition.coeff", "XieBeni.index", "FukuyamaSugeno.index", "Explained.inertia").
#' Other available indices are : "DaviesBoulin.index", "CalinskiHarabasz.index",
#' "GD43.index", "GD53.index" and "Negentropy.index".
#' @param maxiter An integer for the maximum number of iteration
#' @param tol The tolerance criterion used in the evaluateMatrices function for
#'   convergence assessment
#' @param seed An integer used for random number generation. It ensures that the
#' start centers will be the same if the same integer is selected.
#' @param init A string indicating how the initial centers must be selected. "random"
#' indicates that random observations are used as centers. "kpp" use a distance based method
#' resulting in more dispersed centers at the beginning. Both of them are heuristic.
#' @param verbose A boolean indicating if a progressbar should be displayed
#' @return A dataframe with indicators assessing the quality of classifications
#' @export
#' @examples
#' \donttest{
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- sf::st_drop_geometry(LyonIris[AnalysisFields])
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' #set spconsist to TRUE to calculate the spatial consistency indicator
#' #FALSE here to reduce the time during package check
#' values <- select_parameters(algo = "SFCM", dataset, k = 5, m = seq(2,3,0.1),
#'     alpha = seq(0,2,0.1), nblistw = Wqueen, spconsist=FALSE)
#' }
select_parameters <- function(algo,data,k,m,alpha = NA, beta = NA, nblistw=NULL, lag_method="mean", window = NULL,
                              spconsist = TRUE, classidx = TRUE, nrep = 30, indices = NULL,
                              standardize = TRUE,
                              robust = FALSE, noise_cluster = FALSE, delta = NA,
                              maxiter = 500, tol = 0.01,
                              seed=NULL, init = "random", verbose = TRUE){

    if(spconsist==FALSE & classidx==FALSE){
        stop("one of spconsist and classidx must be TRUE")
    }

    if(inherits(nblistw,"list") == FALSE){
        nblistw <- list(nblistw)
    }

    if(inherits(window,"list") == FALSE ){
        window <- list(window)
    }
    if(is.null(indices)){
        indices <- c("Silhouette.index", "Partition.entropy", "Partition.coeff", "XieBeni.index", "FukuyamaSugeno.index", "Explained.inertia")
    }

    allcombinaisons <- expand.grid(k=k,m=m,alpha=alpha,beta=beta,
                                   listsw=1:length(nblistw),lag_method=lag_method,window = 1:length(window),
                                   delta = delta
                                   )

    print(paste("number of combinaisons to estimate : ",nrow(allcombinaisons)))
    dfIndices <- eval_parameters(algo, allcombinaisons, data, nblistw, window, standardize,
                                 robust, noise_cluster,
                                 spconsist, classidx, nrep, indices,
                                 tol, maxiter, seed, init = init, verbose = verbose)
    return(dfIndices)
}


#' @rdname select_parameters
#' @examples
#' \donttest{
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- sf::st_drop_geometry(LyonIris[AnalysisFields])
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
#' tested by specifying a vector : c("mean","median"). When working with rasters,
#' the string must be parsable to a function like mean, min, max, sum, etc. and will
#' be applied to all the pixels values in the window designated by the parameter window
#' and weighted according to the values of this matrix.
#' @param window A list of windows to use to calculate neighbouring values if
#' rasters are used.
#' @param robust A boolean indicating if the "robust" version of the algorithm must be used (see details)
#' @param noise_cluster A boolean indicatong if a noise cluster must be added to the solution (see details)
#' @param delta A float giving the distance of the noise cluster to each observation
#' @param spconsist A boolean indicating if the spatial consistency must be
#' calculated
#' @param classidx A boolean indicating if the quality of classification
#' indices must be calculated
#' @param nrep An integer indicating the number of permutation to do to simulate
#'   the random distribution of the spatial inconsistency. Only used if spconsist
#'   is TRUE.
#' @param indices A character vector with the names of the indices to calculate, to
#' evaluate clustering quality. default is :c("Silhouette.index", "Partition.entropy",
#' "Partition.coeff", "XieBeni.index", "FukuyamaSugeno.index", "Explained.inertia").
#' Other available indices are : "DaviesBoulin.index", "CalinskiHarabasz.index",
#' "GD43.index", "GD53.index" and "Negentropy.index".
#' @param standardize A boolean to specify if the variable must be centered and
#'   reduce (default = True)
#' @param maxiter An integer for the maximum number of iteration
#' @param tol The tolerance criterion used in the evaluateMatrices function for
#'   convergence assessment
#' @param seed An integer used for random number generation. It ensures that the
#' start centers will be the same if the same integer is selected.
#' @param init A string indicating how the initial centers must be selected. "random"
#' indicates that random observations are used as centers. "kpp" use a distance based method
#' resulting in more dispersed centers at the beginning. Both of them are heuristic.
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
#' dataset <- sf::st_drop_geometry(LyonIris[AnalysisFields])
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' future::plan(future::multisession(workers=2))
#' #set spconsist to TRUE to calculate the spatial consistency indicator
#' #FALSE here to reduce the time during package check
#' values <- select_parameters.mc("SFCM", dataset, k = 5, m = seq(1,2.5,0.1),
#'     alpha = seq(0,2,0.1), nblistw = Wqueen, spconsist=FALSE)
#' ## make sure any open connections are closed afterward
#' if (!inherits(future::plan(), "sequential")) future::plan(future::sequential)
#'}
select_parameters.mc <- function(algo,data,k,m,alpha = NA, beta = NA, nblistw=NULL, lag_method="mean", window = NULL,
                                 spconsist = TRUE, classidx = TRUE, nrep = 30, indices = NULL,
                                 standardize = TRUE,
                                 robust = FALSE, noise_cluster = FALSE, delta = NA,
                                 maxiter = 500, tol = 0.01, chunk_size = 5,
                                 seed=NULL, init = "random", verbose = TRUE){

    if(spconsist==FALSE & classidx==FALSE){
        stop("one of spconsist and classidx must be TRUE")
    }

    if(inherits(nblistw, "list")== FALSE){
        nblistw <- list(nblistw)
    }

    if(inherits(window, "list") == FALSE ){
        window <- list(window)
    }
    if(is.null(indices)){
        indices <- c("Silhouette.index", "Partition.entropy", "Partition.coeff", "XieBeni.index", "FukuyamaSugeno.index", "Explained.inertia")
    }

    allcombinaisons <- expand.grid(k=k,m=m,alpha=alpha,beta=beta,delta = delta,
                                   listsw=1:length(nblistw),lag_method=lag_method, window = 1:length(window))

    if (verbose){
        print(paste("number of combinaisons to estimate : ",nrow(allcombinaisons)))
    }
    chunks <- split(1:nrow(allcombinaisons), rep(1:ceiling(nrow(allcombinaisons) / chunk_size),
                                       each = chunk_size, length.out = nrow(allcombinaisons)))

    chunks <- lapply(chunks,function(x){return(allcombinaisons[x,])})

    if(inherits(data, "data.frame") == FALSE){
      ## NOTE HERE : the rasters from terra must be read here befond send to clusters...
      dataw <- lapply(data, terra::wrap)
      wrapped <- TRUE
    }else{
      dataw <- data
      wrapped <- FALSE
    }


    # step2 : starting the function
    iseq <- 1:length(chunks)
    if(verbose){
        progressr::with_progress({
            p <- progressr::progressor(along = iseq)
            values <- future.apply::future_lapply(iseq, function(i) {
                sprintf(algo)
                parameters <- chunks[[i]]
                indices <- eval_parameters(algo, parameters, dataw, nblistw, window, standardize,
                                           robust = robust, noise_cluster = noise_cluster,
                                           spconsist, classidx, nrep, indices,
                                           tol, maxiter, init = init, verbose = FALSE, wrapped = wrapped)
                p(sprintf("i=%g", i))
                return(indices)
            }, future.seed = seed)
        })
    }else{
        values <- future.apply::future_lapply(iseq, function(i) {
            parameters <- chunks[[i]]
            indices <- eval_parameters(algo, parameters, dataw, nblistw, window, standardize,
                                       robust = robust, noise_cluster = noise_cluster,
                                       spconsist, classidx, nrep, indices,
                                       tol, maxiter, init = init, verbose = FALSE, wrapped = wrapped)
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
#' dataset <- sf::st_drop_geometry(LyonIris[AnalysisFields])
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' future::plan(future::multisession(workers=2))
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
#' @param mindist A minimum value for distance between two observations. If two
#' neighbours have exactly the same values, then the euclidean distance between
#' them is 0, leading to an infinite spatial weight. In that case, the minimum
#' distance is used instead of 0.
#'
#' @return A listw object (spdep like)
#' @export
#' @examples
#' data(LyonIris)
#' AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img",
#' "TxChom1564","Pct_brevet","NivVieMed")
#' dataset <- sf::st_drop_geometry(LyonIris[AnalysisFields])
#' queen <- spdep::poly2nb(LyonIris,queen=TRUE)
#' Wqueen <- spdep::nb2listw(queen,style="W")
#' Wqueen2 <- adjustSpatialWeights(dataset,queen,style="C")
adjustSpatialWeights <- function(data,listw,style, mindist = 0.00000000001){
    data <- as.matrix(data)
    new_weights <- lapply(1:nrow(data),function(i){
        row <- data[i,]
        neighbours <- data[listw[[i]],]
        dists <- calcEuclideanDistance2(neighbours,row)
        err <- dists == 0
        if(any(err)){
            dists[err] <- mindist
            warning("Some observartions have exactly the same values as one of their neighbours, leading to
                    an euclidean distance of 0 and an infinite weight. The mindist value is applied instead")
        }
        indists <- 1/dists
        weights <- dists / sum(dists)
        return(weights)
    })
    new_listw <- spdep::nb2listw(listw,glist=new_weights,style = style)
    return(new_listw)
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Utilitary functions #####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Convert categories to membership matrix
#'
#' @description Function to convert a character vector to a membership matrix (binary matrix).
#' The columns of the matrix are ordered with the order function.
#'
#' @param categories A vector with the categories of each observation
#'
#' @return A binary matrix
#' @export
cat_to_belongings <- function(categories){
    cats <- unique(categories)
    cats <- cats[order(cats)]
    cols <- lapply(cats, function(i){
        return (ifelse(categories == i, 1, 0))
    })
    mat <- do.call(cbind, cols)
    colnames(mat) <- cats
    return(mat)
}

#' @rdname cat_to_belongings
#' @export
catToBelongings <- cat_to_belongings


