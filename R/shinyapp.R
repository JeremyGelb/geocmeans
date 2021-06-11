#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### LAUNCHING FUNCTION for shiny apps ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# the variables used in the shiny environment must be declared as globalVariables
globalVariables(c("spatial4326", "mapfun", "variables", "belongings", "n", "mymap",
                  "dataset", "base_violinplots", "dark", "light", "uncertainMap",
                  "base_boxplots","radarchart"
                  ))


#' @title Classification result explorer
#'
#' @description Stat a local Shiny App to explore the results of a classification
#'
#' @param belongings A matrix or a dataframe representing the membership values
#' obtained for each observation
#' @param dataset A dataframe or matrix representing the data used for the
#' classification
#' @param port A integer of length 4 indicating the port on which starting the
#' shiny app. Default is 8100
#' @param spatial A spatial object (SpatialPointsDataFrame, SpatialPolygonsDataFrame or
#' SpatialLinesDataFrame) used to map the observations
#' @param ... Other parameters passed to the function runApp
#' @importFrom leaflet colorBin leaflet addPolygons addPolylines addCircles addLayersControl hideGroup addLegend addProviderTiles colorFactor
#' @importFrom grDevices colorRamp
#' @importFrom plotly plot_ly layout add_markers add_trace
#' @importFrom utils installed.packages
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
#' Data <- LyonIris@data[AnalysisFields]
#' for (Col in names(Data)){
#'   Data[[Col]] <- as.numeric(scale(Data[[Col]]))
#' }
#'
#' Cmean <- CMeans(Data,4,1.5,500,standardize = FALSE, seed = 456, tol = 0.00001, verbose = FALSE)
#'
#' sp_clust_explorer(Cmean$Belongings, Data, LyonIris)
#' }
sp_clust_explorer <- function(belongings, dataset, spatial, port = 8100, ...) {
  print("launching the app")
  appDir <- system.file("shiny-examples", "cluster_explorer", package = "geocmeans")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
  }

  mandatory_packages <- c("shiny", "leaflet", "plotly")
  if(requireNamespace(mandatory_packages) == FALSE){
    stop("The shiny app can be used only if the packages shiny, leaflet and plotly are installed ! \n
         We also recommand to install shinyWidgets, bslib and car for an optimal experience.
         ")
  }
  secondary_packages <- c("shinyWidgets", "bslib", "car")
  if(requireNamespace(secondary_packages) == FALSE){
    warning("We recommand to install the packages shinyWidgets, bslib and car for an optimal experience with
            this shiny app")
  }

  shiny_env <- new.env()

  if(is.matrix(dataset)){
    oldnames <- colnames(dataset)
    dataset <- as.data.frame(dataset)
    if(is.null(oldnames)){
      oldnames <- paste("var",1:ncol(dataset),sep="")
    }
    names(dataset) <- oldnames
  }

  # the colors to use for the groups
  colors <- c("#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD","#8C564B",
              "#E377C2","#7F7F7F","#BCBD22","#17BECF","#AEC7E8","#FFBB78",
              "#98DF8A","#FF9896","#C5B0D5","#C49C94","#F7B6D2","#C7C7C7",
              "#DBDB8D","#9EDAE5")[1:ncol(belongings)]

  assign('colors', colors, .GlobalEnv)

  # the available themes
  if("bslib" %in% installed.packages()){
    dark <-  bslib::bs_theme(bootswatch = "darkly", version = "3")
    light <- bslib::bs_theme(version = "3")
  }else{
    dark <- NULL
    light <- NULL
  }

  assign('dark', dark, .GlobalEnv)
  assign('light', light, .GlobalEnv)

  ## enregistrement des principaux objets dans un environnement
  ## shiny env qui pourra etre utilise dans l'UI et dans le SERVER
  assign('belongings', belongings, .GlobalEnv)
  assign('dataset', dataset, .GlobalEnv)
  assign('spatial', spatial, .GlobalEnv)

  ## determiner les groupes que dont on dispose et les donnees
  groups <- paste("group ", 1:ncol(belongings), sep = "")
  variables <- names(dataset)
  assign('groups', groups, .GlobalEnv)
  assign('variables', variables, .GlobalEnv)

  ## pour leaflet, les coordonnees doivent etre en 4326
  ref <- sp::CRS("+init=epsg:4326")
  spatial4326 <- sp::spTransform(spatial, ref)

  ## preparer la carte leaflet du premier panneau ***************************************

  mymap <- leaflet(height = "600px") %>%
    addProviderTiles(leaflet::providers$Stamen.TonerBackground, group = "Toner Lite", layerId = "back1") %>%
    addProviderTiles(leaflet::providers$OpenStreetMap, group = "Open Street Map", layerId = "back2")

  if(class(spatial)[[1]] == "SpatialPolygonsDataFrame"){
    mapfun <- function(map, data, weight, group, color, fillColor, layerId, ...){

      map %>% addPolygons(
        data = data,
        weight = weight,
        group = group,
        color = color,
        fillColor = fillColor,
        layerId = layerId,
        ...
      )
    }
  }else if (class(spatial)[[1]] == "SpatialPointsDataFrame"){
    mapfun <- addCircles
  }else if (class(spatial)[[1]] == "SpatialLinesDataFrame"){
    mapfun <- function(map, data, weight, group, color, fillColor, layerId, ...){
      if(is.null(fillColor)){
        fillColor <- "red"
      }
      map %>% addPolylines(
        data = data,
        weight = 3,
        group = group,
        color = fillColor,
        layerId = layerId,
        ...
      )
    }
  }else{
    stop("spatial must be one of SpatialPolygonsDataFrame, SpatialPointsDataFrame, SpatialLinesDataFrame")
  }

  # ajouter les layers
  for (i in 1:ncol(belongings)){

    bins <- seq(0,1,0.1)
    cols <- colorRamp(c("#FFFFFF", colors[[i]]), interpolate = "spline")
    pal <- leaflet::colorBin(cols, c(0,1), bins = bins)

    mymap <- mymap %>% mapfun(data = spatial4326,
                              weight = 1,
                              group = paste("group ",i,sep=""),
                              color = "black",
                              fillColor = ~pal(belongings[,i]),
                              fillOpacity = 0.7,
                              layerId = 1:nrow(spatial4326)) %>%
      addLegend(pal = pal, values = bins, opacity = 0.7,
                title = NULL, group=paste("group ",i,sep=""),
                position = "bottomright")
  }

  ## ajouter un layer de hard clustering
  colnames(belongings) <- paste("group",1:ncol(belongings), sep = " ")
  groups <- colnames(belongings)[max.col(belongings, ties.method = "first")]
  spatial4326$group <- as.factor(groups)

  factpal <- colorFactor(colors, spatial4326$group)

  mymap <- mymap %>% mapfun(data = spatial4326,
                            weight = 1,
                            group = "Most likely group",
                            color = "black",
                            fillColor = ~factpal(spatial4326$group),
                            fillOpacity = 0.7,
                            layerId = 1:nrow(spatial4326)) %>%
    addLegend(pal = factpal, values = spatial4326$group, opacity = 0.7,
              title = NULL, group= "Most likely group",
              position = "bottomright")


  # ajouter les enjolivages
  mymap <- mymap %>%
    addLayersControl(
      position = "bottomleft",
      baseGroups = c("Toner Lite","Open Street Map"),
      overlayGroups  = c(paste("group ", 1:ncol(belongings), sep = ""),"Most likely group"),
      options = leaflet::layersControlOptions(collapsed = FALSE))

  for(i in 2:ncol(belongings)){
    mymap <- mymap %>% hideGroup(paste("group ",i,sep=""))
  }
  mymap <- mymap %>% hideGroup("Most likely group")

  assign('mymap', mymap, .GlobalEnv)
  assign('mapfun', mapfun, .GlobalEnv)
  assign('spatial4326', spatial4326, .GlobalEnv)

  ## preparer la carte leaflet du troisieme panneau ***************************************
  uncertainMap <- leaflet(height = "600px") %>%
    addProviderTiles(leaflet::providers$Stamen.TonerBackground, group = "Toner Lite", layerId = "back1") %>%
    addProviderTiles(leaflet::providers$OpenStreetMap, group = "Open Street Map", layerId = "back2")

  bins <- c(0,0.3,0.6,0.9,1)
  cols <- colorRamp(c("#FFFFFF", "#D30000"), interpolate = "linear")
  pal <- leaflet::colorBin(cols, c(0,1), bins = bins)

  # layer of uncertaintyvalues
  uncertain_values <- calcUncertaintyIndex(belongings)

  uncertainMap <- uncertainMap %>% mapfun(data = spatial4326,
                                          weight = 1,
                                          group = "UncertaintyIdx",
                                          color = "black",
                                          fillColor = ~pal(uncertain_values),
                                          fillOpacity = 0.7,
                                          layerId = 1:nrow(spatial4326)) %>%
    addLegend(pal = pal, values = bins, opacity = 0.7,
              title = NULL, group="UncertaintyIdx",
              position = "bottomright")

  # binary layer of uncertain values
  values <- apply(belongings, 1, max) < 0.45
  spdf <- subset(spatial4326, values)
  uncertainMap <- uncertainMap %>% mapfun(data = spdf,
                                          weight = 1,
                                          group = "binaryUncertain",
                                          color = "black",
                                          fillColor = "red",
                                          fillOpacity = 0.7,
                                          layerId = 1:nrow(spdf)) %>%
    addLegend(opacity = 0.7,
              colors = c("#D30000"),
              labels = "uncertain observations",
              title = NULL, group="binaryUncertain",
              position = "bottomright") %>%
    addLayersControl(
      position = "bottomleft",
      baseGroups = c("Toner Lite", "Open Street Map"),
      overlayGroups  = c("UncertaintyIdx","binaryUncertain"),
      options = leaflet::layersControlOptions(collapsed = FALSE)) %>%
    hideGroup("UncertaintyIdx")

  assign('uncertainMap', uncertainMap, .GlobalEnv)
  ##******************************************************************


  shiny::runApp(appDir, display.mode = "normal",port = 8100)
}


