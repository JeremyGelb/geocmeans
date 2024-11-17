#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### LAUNCHING FUNCTION for shiny apps ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# the variables used in the shiny environment must be declared as globalVariables
# nocov start
globalVariables(c("spatial4326", "mapfun", "variables", "belongings", "n", "mymap",
                  "dataset", "base_violinplots", "dark", "light", "uncertainMap",
                  "base_boxplots","radarchart", 'rasterMode', "object","shiny_data"
                  ))


#' @title geocmeans general environment
#'
#' @description An environment used by geocmeans to store data, functions and values
#' @keywords internal
geocmeans_env <- new.env()

#see here to remove global environment variables
#https://community.rstudio.com/t/alternative-to-global-for-passing-variable-to-shiny-app/26476/2


#' @title Classification result explorer
#'
#' @description Start a local Shiny App to explore the results of a classification
#'
#' @param object A FCMres object, typically obtained from functions CMeans,
#'   GCMeans, SFCMeans, SGFCMeans
#' @param spatial A feature collection (sf) used to map the observations. Only needed if object was not created
#' from rasters.
#' @param membership A matrix or a dataframe representing the membership values
#' obtained for each observation. If NULL, then the matrix is extracted from
#' object.
#' @param dataset A dataframe or matrix representing the data used for the
#' classification. If NULL, then the matrix is extracted from object.
#' @param port An integer of length 4 indicating the port on which to start the
#' Shiny app. Default is 8100
#' @param ... Other parameters passed to the function runApp
#' @importFrom leaflet addRasterImage colorBin leaflet addPolygons addPolylines addCircles addLayersControl hideGroup addLegend addProviderTiles colorFactor
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
#' Data <- sf::st_drop_geometry(LyonIris[AnalysisFields])
#' for (Col in names(Data)){
#'   Data[[Col]] <- as.numeric(scale(Data[[Col]]))
#' }
#'
#' Cmean <- CMeans(Data,4,1.5,500,standardize = FALSE, seed = 456, tol = 0.00001, verbose = FALSE)
#'
#' sp_clust_explorer(Cmean, LyonIris)
#' }
sp_clust_explorer <- function(object = NULL, spatial = NULL, membership = NULL, dataset = NULL, port = 8100, ...) {

  # if(object$isRaster){
  #   stop("The shiny app can not be used currently to display results from raster data, sorry...")
  # }

  # creating a list to store all the data to pass to the shiny app
  shiny_data <- list()

  # checking if the directory of hte shiny app is here  ---------------------------------------
  appDir <- system.file("shiny-examples", "cluster_explorer", package = "geocmeans")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `geocmeans`.", call. = FALSE)
  }

  # checking if the mandatory packages are installed ----------------------------------------
  mandatory_packages <- c("shiny", "leaflet", "plotly")
  if(requireNamespace(mandatory_packages) == FALSE){
    stop("The shiny app can be used only if the packages shiny, leaflet and plotly are installed ! \n
         We also recommand to install shinyWidgets, bslib and car for an optimal experience.
         ")
  }

  # checking if the not necessary but usefull packages are installed ----------------------------------------
  secondary_packages <- c("shinyWidgets", "bslib", "car", "shinyhelper")
  if(requireNamespace(secondary_packages) == FALSE){
    warning("We recommand to install the packages shinyWidgets, bslib, car and shinyhelper for an optimal experience with
            this shiny app")
  }

  # checking if the objects given have the right informations ----------------------------------------
  if(is.null(object) == FALSE){
    if(is.FCMres(object) == FALSE){
      stop("If object is not NULL, it must be an object of class FCMres (see help(is.FCMres))")
    }
  }

  if(is.null(object) & is.null(spatial)){
    stop("if object is NULL, spatial must be specified")
  }

  ok_sp <- c("sf")
  if(is.null(object) == FALSE){
    if(object$isRaster == FALSE){
      if(inherits(spatial, ok_sp)  == FALSE){
        stop('spatial must be a feature collection (sf) because object was not created with rasters')
      }
    }
  }else{
    if(inherits(spatial, ok_sp)  == FALSE){
      stop('spatial must be one a feature collection (sf)')
    }
  }

  if(is.null(object) & (is.null(membership) | is.null(dataset))){
    stop("either object or both dataset and membership must be specified")
  }

  if(is.null(object) == FALSE){
    if(is.null(dataset)){
      dataset <- object$Data
    }
    if(is.null(membership)){
      belongings <- object$Belongings
    }else{
      belongings <- membership
    }
    inertia <- calcexplainedInertia(object$Data, object$Belongings)
  }else{
    inertia <- calcexplainedInertia(dataset, belongings)
  }

  #assign('inertia', inertia, .GlobalEnv)
  shiny_data$inertia <- inertia

  # # Adjusting the membership matrix if we have a noise cluster ----------------------------------------
  noise_mode <- FALSE
  if(is.null(object) == FALSE){
    if(is.null(object$noise_cluster) == FALSE){
      #test <- matrixStats::rowMaxs(object$Belongings) < object$noise_cluster
      #object$Belongings <- cbind(object$Belongings, object$noise_cluster)
      #object$Groups[test] <- "noise"
      noise_vec <- object$noise_cluster
      noise_mode <- TRUE
    }
  }



  # Preparing some global variables for the app ----------------------------------------
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

  #assign('colors', colors, .GlobalEnv)
  shiny_data$colors <- colors

  # the available themes
  if("bslib" %in% installed.packages()){
    dark <-  bslib::bs_theme(bootswatch = "darkly", version = "3")
    light <- bslib::bs_theme(version = "3")
  }else{
    dark <- NULL
    light <- NULL
  }

  #assign('dark', dark, .GlobalEnv)
  #assign('light', light, .GlobalEnv)
  shiny_data$dark <- dark
  shiny_data$light <- light

  # check if we have to deal with rasters
  rasterMode <- FALSE
  if(is.null(object) == FALSE){
    if(object$isRaster){
      rasterMode <- TRUE

      # let me check if the installed version of leaflet is enough
      if(utils::packageVersion("leaflet") <= "2.1.1"){
        stop("To use the shiny app with raster data, you must have installed a version of leaflet > 2.1.1 (ex 2.1.1.9000 from github)")
      }

    }
  }

  if(rasterMode){
    ## saving the main objets in the global environment for them to
    ## be used in the main functions UI and SERVER
    ## but reduce there size
    Ids <- sample(1:nrow(belongings), size = 1500, replace = FALSE)
    #assign('belongings', belongings[Ids,], .GlobalEnv)
    #assign('dataset', dataset[Ids,], .GlobalEnv)
    shiny_data$belongings <- belongings
    shiny_data$dataset <- dataset

    ## creating a referencing raster with the right projection
	  ref_raster <- terra::project(object$rasters[[1]], y = "epsg:3857", method = "near")
    #assign('ref_raster', ref_raster, .GlobalEnv)
    shiny_data$ref_raster <- ref_raster

    old_names <- names(object$rasters)
    object$rasters <- lapply(object$rasters, function(rast){
	  terra::project(rast,  y = "epsg:3857", method = "near")
    })
    names(object$rasters) <- old_names

  }else{
    ## saving the main objets in the global environment for them to
    ## be used in the main functions UI and SERVER
    #assign('belongings', belongings, .GlobalEnv)
    #assign('dataset', dataset, .GlobalEnv)
    #assign('spatial', spatial, .GlobalEnv)

    shiny_data$belongings <- belongings
    shiny_data$dataset <- dataset
    shiny_data$spatial <- spatial

    ## for leaflet, the CRS must be 4326
    # ref <- sp::CRS("+init=epsg:4326")
	  ref <- sf::st_crs(4326)
    #spatial4326 <- sp::spTransform(spatial, ref)
	  spatial4326 <- sf::st_transform(spatial, ref)
    #assign('spatial4326', spatial4326, .GlobalEnv)
    shiny_data$spatial4326 <- spatial4326
  }

  #assign("rasterMode", rasterMode, .GlobalEnv)
  shiny_data$rasterMode <- rasterMode


  groups <- paste("group ", 1:ncol(belongings), sep = "")
  variables <- names(dataset)
  #assign('groups', groups, .GlobalEnv)
  #assign('variables', variables, .GlobalEnv)
  shiny_data$groups <- groups
  shiny_data$variables <- variables

  ## prepare the leaflet maps in the first pannel ***************************************

  mymap <- leaflet(height = "600px") %>%
    addProviderTiles(leaflet::providers$Stadia.StamenTonerBackground, group = "Toner Lite", layerId = "back1") %>%
    addProviderTiles(leaflet::providers$OpenStreetMap, group = "Open Street Map", layerId = "back2")

  if(rasterMode == FALSE){
    geom_type <- sf::st_geometry_type(spatial, by_geometry = FALSE)
  }else{
    geom_type <- "raster"
  }


  if(geom_type %in% c("POLYGON", "MULTIPOLYGON")){
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
  }else if (geom_type %in% c("POINT", "MULTIPOINT")){
    mapfun <- addCircles
  }else if (geom_type %in% c("LINESTRING", "MULTILINESTRING")){
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
  }else if (object$isRaster){
    #nothing to do here if we have to plot rasters
    i <- 1
  }else{
    stop("spatial must be a feature collections (POINT, MULTIPOINT, LINESTRING, MULTILINESTRING, POLYGON or MULTIPOLYGON)")
  }

  # adding the layers if we are in vector mode
  if(rasterMode == FALSE){

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
    ## if required, a layer with the noise values
    if(noise_mode){
      bins <- seq(0,1,0.1)
      cols <- colorRamp(c("#FFFFFF", "black"), interpolate = "spline")
      pal <- leaflet::colorBin(cols, c(0,1), bins = bins)

      mymap <- mymap %>% mapfun(data = spatial4326,
                                weight = 1,
                                group = "noise",
                                color = "black",
                                fillColor = ~pal(noise_vec),
                                fillOpacity = 0.7,
                                layerId = 1:nrow(spatial4326)) %>%
        addLegend(pal = pal, values = bins, opacity = 0.7,
                  title = NULL, group="noise",
                  position = "bottomright")
    }

    ## and a layer for the hard partition
    colnames(belongings) <- paste("group",1:ncol(belongings), sep = " ")
    groups <- colnames(belongings)[max.col(belongings, ties.method = "first")]
    if(noise_mode){
      groups <- ifelse(matrixStats::rowMaxs(belongings) < noise_vec, "noise", groups)
    }
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

    #assign('mapfun', mapfun, .GlobalEnv)
    shiny_data$mapfun <- mapfun

  }else{
    # IF WE ARE IN RASTER MODE

    # adding all the groups
    name <- names(object$rasters)
    i <- 1
    ok_names <- name[grepl("group",name, fixed = TRUE)]
    for (name in ok_names){
      rast <- object$rasters[[name]]
      vals <- terra::values(rast, mat = FALSE)
      pal <- leaflet::colorNumeric(c("#FFFFFF", colors[[i]]),
                                   vals, na.color = "transparent")
      mymap <- mymap %>%
        addRasterImage(rast, colors = pal, opacity = 0.8,
                       group = paste("group ",i,sep="")) %>%
        addLegend(pal = pal, values = vals, opacity = 0.7,
                  title = NULL, group = paste("group ",i,sep=""),
                  position = "bottomright")
      i <- i + 1
    }

    if(noise_mode){
      pal <- leaflet::colorNumeric(c("#FFFFFF", "black"),
                                   vals, na.color = "transparent")
      noise_rast <- terra::project(object$noise_cluster,  y = "epsg:3857", method = "near")
      mymap <- mymap %>%
        addRasterImage(noise_rast, colors = pal, opacity = 0.8,
                       group = "noise") %>%
        addLegend(pal = pal, values = vals, opacity = 0.7,
                  title = NULL, group = "noise",
                  position = "bottomright")
    }

    # adding the last layer with the most likely groups
    rast <- object$rasters$Groups
    vals <- terra::values(rast, mat = FALSE)
    pal <- leaflet::colorNumeric(colors[1:ncol(object$Belongings)],
                                 vals, na.color = "transparent")
    mymap <- mymap %>%
      addRasterImage(rast, colors = pal, opacity = 0.8,
                     group= "Most likely group",) %>%
      addLegend(pal = pal, values = vals, opacity = 0.7,
                title = NULL, group= "Most likely group",
                position = "bottomright")

    #assign('mapfun', NULL, .GlobalEnv)
    #shiny_data$mapfun <- mapfun

  }

  # adding some tools for the map
  if(noise_mode){
    all_layer <- c(paste("group ", 1:ncol(belongings), sep = ""),"noise","Most likely group")
  }else{
    all_layer <- c(paste("group ", 1:ncol(belongings), sep = ""),"Most likely group")
  }

  mymap <- mymap %>%
    addLayersControl(
      position = "bottomleft",
      baseGroups = c("Toner Lite","Open Street Map"),
      overlayGroups  = all_layer,
      options = leaflet::layersControlOptions(collapsed = FALSE))

  for(i in 2:ncol(belongings)){
    mymap <- mymap %>% hideGroup(paste("group ",i,sep=""))
  }
  mymap <- mymap %>% hideGroup("Most likely group")
  if(noise_mode){
    mymap <- mymap %>% hideGroup("noise")
  }

  #assign('mymap', mymap, .GlobalEnv)
  shiny_data$mymap <- mymap


  ## preparing the map for the third pannel ***************************************
  uncertainMap <- leaflet(height = "600px") %>%
    addProviderTiles(leaflet::providers$Stadia.StamenTonerBackground, group = "Toner Lite", layerId = "back1") %>%
    addProviderTiles(leaflet::providers$OpenStreetMap, group = "Open Street Map", layerId = "back2")

  if(rasterMode == FALSE){
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

  }else{
    # IF WE ARE IN RASTER MODE
    all_values <- lapply(object$rasters[ok_names], function(rast){
      terra::values(rast, mat = FALSE)[object$missing]
    })
    maxs <- do.call(pmax,all_values)
    rast <- object$rasters[[1]]
    vals <- rep(0, times = terra::ncell(rast))
    vals[object$missing] <- maxs
    vals <- ifelse(vals < 0.45, 1,0)
    vals[!object$missing] <- NA
    terra::values(rast) <- vals

    pal <- leaflet::colorNumeric(c("#FFFFFF","#D30000"),
                                 vals, na.color = "transparent")
    uncertainMap <- uncertainMap %>%
      addRasterImage(rast, colors = pal, opacity = 0.8,
                     group= "binaryUncertain") %>%
      addLegend(pal = pal, values = vals, opacity = 0.7,
                title = NULL, group= "binaryUncertain",
                position = "bottomright",
                labels = "uncertain observations") %>%
      addLayersControl(
        position = "bottomleft",
        baseGroups = c("Toner Lite", "Open Street Map"),
        overlayGroups  = c("binaryUncertain"),
        options = leaflet::layersControlOptions(collapsed = FALSE))

  }


  #assign('uncertainMap', uncertainMap, .GlobalEnv)
  #assign('object', object, .GlobalEnv)
  shiny_data$uncertainMap <- uncertainMap
  shiny_data$object <- object

  assign('shiny_data', shiny_data, geocmeans_env)
  ##******************************************************************


  shiny::runApp(appDir, display.mode = "normal",port = 8100)
}
# nocov end

