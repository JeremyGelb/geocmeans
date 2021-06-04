#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### EN DEVELOPPEMENT ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### PLOTING FUNCTION ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Add an horizontal line to plotly
#'
#' @description Add an horizontal line to plotly
#' @param y The y value for the line
#' @param color The color of the line
#' @examples
#' #This is an internal function, no example provided
hline <- function(y = 0, color = "blue") {
  list(
    type = "line",
    x0 = 0,
    x1 = 1,
    xref = "paper",
    y0 = y,
    y1 = y,
    line = list(color = color)
  )
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### SERVER FUNCTION ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @title Shiny App server
#'
#' @description Set the server for the Shiny App
#' @param input The shiny input object
#' @param output The shiny output object
#' @param session The shiny session object
#'
#' @importFrom shiny reactive observeEvent
#' @importFrom leaflet renderLeaflet leafletProxy removeShape
#' @importFrom plotly renderPlotly plot_ly layout
#' @importFrom grDevices colorRampPalette
#' @examples
#' #This is an internal function, no example provided
shiny_server <- function(input, output, session) {

  ## creating the map man
  output$mymap <- renderLeaflet({
    mymap
  })

  ## and the first bivariate plot

  bivar_params <- reactive({
    list(x = input$var1_biplot,
         y = input$var2_biplot,
         color = input$group_biplot)
  })


  output$bivar_plot <- renderPlotly({
    params <- bivar_params()
    gpcol <- as.numeric(strsplit(params$color," ", fixed = TRUE)[[1]][[2]])
    idx <- order(belongings[,gpcol])
    biplot <- plot_ly(
      x = dataset[[params$x]][idx],
      y = dataset[[params$y]][idx],
      color = belongings[,gpcol][idx],
      colors = colorRampPalette(c("white", colors[[gpcol]]))(10),
      type = "scatter",
      mode = "markers",
      size = 2
    )

    # adjusting the color with the theme
    if(is.null(light) == FALSE){
      if (isTRUE(input$dark_mode)){
        biplot <- biplot %>% layout(
          font = list(color = "white"),
          paper_bgcolor = "#303030",
          plot_bgcolor = "#303030")
      }
    }
    biplot
  })


  ## selecting the right number of columns and rows for the ggarrange
  nv <- length(variables)
  nc <- floor(16 / ncol(belongings))
  rest <-  nv %% nc
  nr <- (nv-rest) / nc
  if(rest > 0){
    nr <- nr + 1
  }

  ## putting the original violinplots
  lapply(1:length(base_violinplots), function(i){
    vplot <- base_violinplots[[i]]
    output[[paste("violinplots",i,sep="")]] <- renderPlotly({vplot})
  })



  firsttime <- TRUE
  ## ----------here is an observer working when we click on the map--------------
  observeEvent(input$mymap_shape_click, {
    p <- input$mymap_shape_click

    ## Step1 : rendering the plot for the belongings
    df <- data.frame(
      values = belongings[p$id,],
      groups = paste("group ", 1:length(belongings[p$id,]), sep = "")
    )

    output$barplot1 <- renderPlotly({
      plot_ly(
        x = df$groups,
        y = df$values,
        type = "bar",
        name = "Groups membership values"
      )

    })

    # ## Step3 : updating the violin plots
    new_violins <- lapply(1:length(variables), function(i){
      violin <- base_violinplots[[i]]
      varname <- variables[[i]]
      value <- dataset[p$id,varname][[1]]
      if (isTRUE(input$dark_mode)){
        violin2 <- violin %>% layout(
          font = list(color = "white"),
          paper_bgcolor = "#303030",
          plot_bgcolor = "#303030",
          shapes = list(hline(value, color = "red")))
      }else{
        violin2 <- violin %>% layout(shapes = list(hline(value, color = "red")))
      }

      return(violin2)
    })
    lapply(1:length(new_violins), function(i){
      vplot <- new_violins[[i]]
      output[[paste("violinplots",i,sep="")]] <- renderPlotly({vplot})
    })

    ## Step4 : adding the highlight on the selected feature
    feat <- spatial4326[p$id,]
    if(firsttime){
      leafletProxy('mymap') %>%
        mapfun(data = feat, weight = 2, opacity = 1.0, fillOpacity = 0, color = "red",
                    layerId = "highlighter")
      firsttime <- FALSE
    }else{
      leafletProxy('mymap') %>%
        removeShape("highlighter") %>%
        mapfun(data = feat, weight = 2, opacity = 1.0, fillOpacity = 0, color = "red",
                    layerId = "highlighter")
    }


  })

  ## ----------here is an observer listening for the ellipsis--------------
  observeEvent(input$show_ellipsis,{

    # redraw the biplot
    output$bivar_plot <- renderPlotly({
      params <- bivar_params()
      gpcol <- as.numeric(strsplit(params$color," ", fixed = TRUE)[[1]][[2]])
      idx <- order(belongings[,gpcol])
      biplot <- plot_ly(
        x = dataset[[params$x]][idx],
        y = dataset[[params$y]][idx],
        color = belongings[,gpcol][idx],
        colors = colorRampPalette(c("white", colors[[gpcol]]))(10),
        type = "scatter",
        mode = "markers",
        size = 2
      )

      # adjusting the color with the theme
      if (is.null(light) == FALSE){
        if (isTRUE(input$dark_mode)){
          biplot <- biplot %>% layout(
            font = list(color = "white"),
            paper_bgcolor = "#303030",
            plot_bgcolor = "#303030")
        }
      }

      # adding the ellipsis if required
      if(isTRUE(input$show_ellipsis)){
        for(j in 1:ncol(belongings)){
          coords <- car::dataEllipse(dataset[[params$x]],
                                     dataset[[params$y]],
                                     weights = belongings[,j],
                                     levels = 0.75,
                                     draw = FALSE,
          )
          coords <- coords[1:(nrow(coords)-1),]
          coords <- rbind(coords,coords[1,])
          biplot <- biplot %>%
            add_paths(
              x = coords[,1],
              y = coords[,2],
              line = list(width = 2),
              color = I(colors[[j]]),

            )
        }
      }

      biplot
    })

  })

  ## ----------here is an observer if we can theme switch--------------
  if(is.null(light) == FALSE){
    observe(session$setCurrentTheme(
      {
      # adjust the map (remove a previous overlay)
      if(firsttime == FALSE){
        firsttime <- TRUE
        leafletProxy('mymap') %>%
          removeShape("highlighter")
      }

      # set the right colors for the biplot
      output$bivar_plot <- renderPlotly({
        params <- bivar_params()
        gpcol <- as.numeric(strsplit(params$color," ", fixed = TRUE)[[1]][[2]])
        idx <- order(belongings[,gpcol])
        biplot <- plot_ly(
          x = dataset[[params$x]][idx],
          y = dataset[[params$y]][idx],
          color = belongings[,gpcol][idx],
          colors = colorRampPalette(c("white", colors[[gpcol]]))(10),
          type = "scatter",
          mode = "markers",
          size = 2
        )

        # adjusting the color with the theme
        if (isTRUE(input$dark_mode)){
          biplot <- biplot %>% layout(
            font = list(color = "white"),
            paper_bgcolor = "#303030",
            plot_bgcolor = "#303030")
        }

        # adding the ellipsis if required
        if(isTRUE(input$show_ellipsis)){
          for(j in 1:ncol(belongings)){
            coords <- car::dataEllipse(dataset[[params$x]],
                                       dataset[[params$y]],
                                       weights = belongings[,j],
                                       levels = 0.75,
                                       draw = FALSE,
                                       )
            coords <- coords[1:(nrow(coords)-1),]
            coords <- rbind(coords,coords[1,])
            biplot <- biplot %>%
              add_paths(
                x = coords[,1],
                y = coords[,2],
                line = list(width = 2),
                color = I(colors[[j]])
              )
          }
        }

        biplot
      })

      # set the right colors for violin plots
      new_violins <- lapply(1:length(variables), function(i){
        violin <- base_violinplots[[i]]
        if (isTRUE(input$dark_mode)){
          violin2 <- violin %>% layout(
            font = list(color = "white"),
            paper_bgcolor = "#303030",
            plot_bgcolor = "#303030")
        }else{
          violin2 <- violin
        }
        output[[paste("violinplots",i,sep="")]] <- renderPlotly({violin2})

      })

      if (isTRUE(input$dark_mode)) dark else light
      }
    ))
  }

}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### UI FUNCTIONS ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @title Set the violin plot pannel
#'
#' @description A simple function to prepare the violin plot pannel
#' @param n The number of plots
#' @param nr The number of rows
#' @param nc The number of columns
#'
#' @examples
#' #This is an internal function, no example provided
violinplots_ui <- function(n, nr, nc){
  fulltxt <- "wellPanel("
  cnt <- 1
  colwidth <- 12 / nc
  for(i in 1:nr){
    new_row <-"fluidRow("
    for(j in 1:nc){
      new_col <- paste("column(plotlyOutput('violinplots",cnt,"'), width = ",colwidth,"),",sep="")
      new_row <- paste(new_row, new_col, sep = "")
      cnt <- cnt + 1
    }
    new_row <- paste(substr(new_row,1,nchar(new_row)-1),")",sep="")
    fulltxt <- paste(fulltxt, new_row, ",", sep = "")
  }
  fulltxt <- paste(substr(fulltxt,1,nchar(fulltxt)-1),")",sep="")
  return(eval(parse(text=fulltxt)))
}



#' @title Shiny App UI
#'
#' @description Set the UI for the Shiny App
#'
#' @importFrom shiny fluidPage tabsetPanel tabPanel titlePanel fluidRow wellPanel column selectInput
#' @importFrom leaflet leafletOutput
#' @importFrom plotly plotlyOutput
#' @examples
#' #This is an internal function, no example provided
shiny_ui <- function() {

  ## selecting the right number of columns and rows for the violinplots
  nv <- length(variables)
  nc <- floor(16 / ncol(belongings))
  rest <-  nv %% nc
  nr <- (nv-rest) / nc
  if(rest > 0){
    nr <- nr + 1
  }


  fluidPage(
    tabsetPanel(
      tabPanel("Interactive map", fluid = TRUE,
        ##------------------- PANNEL 1 : for the interactive map------------------
        titlePanel('Interactive map'),

        # Grid Layout
        {if(is.null(light)){
          fluidRow(wellPanel("Welcome in the classification explorer ! Please, click on a feature on the map to start exploring the results of your classification"))
        }else{
          fluidRow(column(width = 4, wellPanel("Welcome in the classification explorer ! Please, click on a feature on the map to start exploring the results of your classification")),
                   column(width = 2, checkboxInput("dark_mode", "Dark mode"))
                   )
        }},
        fluidRow(column(width = 8, wellPanel(leafletOutput("mymap"))),
                 column(width = 4, wellPanel(plotlyOutput("barplot1")))
                 ),
        ## for the violin plots, with have to build a more complex environment
        fluidRow(violinplots_ui(n, nr, nc), height = paste(200*nc,"px",sep="")),
      ),

      tabPanel("Bivariate plot", fluid = TRUE,
               ##------------------- PANNEL 2 : for the bivariate plots------------------
               titlePanel('Interactive bivariate plots'),

               # Grid Layout
               fluidRow(wellPanel("In this panel, you can explore bivariate relationships for the different groups obtained")),
               fluidRow(column(width = 2,
                               selectInput("var1_biplot", "X axis variable", variables, selected = variables[[1]]),
                               selectInput("var2_biplot", "Y axis variable", variables, selected = variables[[2]]),
                               selectInput("group_biplot", "group membership for color", paste("group ", 1:ncol(belongings), sep = "")),
                               {
                                 if("car" %in% installed.packages()){
                                   checkboxInput("show_ellipsis", "show ellipsis")
                                 }
                               },
                               ),
                        column(width = 10,plotlyOutput("bivar_plot", height = "800px")),
               ),
      )
    ),
    theme = light
  )
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### LAUNCHING FUNCTION ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# the variables used in the shiny environment must be declared as globalVariables
globalVariables(c("spatial4326", "mapfun", "variables", "belongings", "n", "mymap",
                  "dataset", "base_violinplots", "dark", "light"
                  ))

#' @title Classification result explorer
#'
#' @description Stat a local Shiny App to explore the results of a classification
#'
#' @param belongings A matrix or a dataframe representing the membership values
#' obtained for each observation
#' @param dataset A dataframe or matrix representing the data used for the
#' classification
#' @param ... Other parameters passed to the function runApp
#' @param spatial A spatial object (SpatialPointsDataFrame, SpatialPolygonsDataFrame or
#' SpatialLinesDataFrame) used to map the observations
#' @importFrom leaflet colorBin leaflet addPolygons addCircles addLayersControl hideGroup addLegend addProviderTiles
#' @importFrom grDevices colorRamp
#' @importFrom plotly plot_ly layout
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
sp_clust_explorer <- function(belongings, dataset, spatial, ...) {

  shiny_env <- new.env()

  # the colors to use for the groups
  colors <- c("#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD","#8C564B",
              "#E377C2","#7F7F7F","#BCBD22","#17BECF","#AEC7E8","#FFBB78",
              "#98DF8A","#FF9896","#C5B0D5","#C49C94","#F7B6D2","#C7C7C7",
              "#DBDB8D","#9EDAE5")[1:ncol(belongings)]

  # the available themes
  if("bslib" %in% installed.packages()){
    dark <-  bslib::bs_theme(bootswatch = "darkly", version = "3")
    light <- bslib::bs_theme(version = "3")
  }else{
    dark <- NULL
    light <- NULL
  }

  assign('dark', dark, shiny_env)
  assign('light', light, shiny_env)

  ## enregistrement des principaux objets dans un environnement
  ## shiny env qui pourra etre utilise dans l'UI et dans le SERVER
  assign('belongings', belongings, shiny_env)
  assign('dataset', dataset, shiny_env)
  assign('spatial', spatial, shiny_env)

  ## determiner les groupes que dont on dispose et les donnees
  groups <- paste("group ", 1:ncol(belongings), sep = "")
  variables <- names(dataset)
  assign('groups', groups, shiny_env)
  assign('variables', variables, shiny_env)

  ## pour leaflet, les coordonnees doivent etre en 4326
  ref <- sp::CRS("+init=epsg:4326")
  spatial4326 <- sp::spTransform(spatial, ref)

  ## preparer la carte leaflet ***************************************

  mymap <- leaflet(height = "600px") %>%
    addProviderTiles(leaflet::providers$Stamen.TonerBackground, group = "Toner Lite")

  if(class(spatial)[[1]] == "SpatialPolygonsDataFrame"){
    mapfun <- addPolygons
  }else if (class(spatial)[[1]] == "SpatialPointsDataFrame"){
    mapfun <- addCircles
  }else if (class(spatial)[[1]] == "SpatialLinesDataFrame"){
    mapfun <- addPolygons
    spatial <- rgeos::gBuffer(spatial, width = 0.01, byid = TRUE)
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

  # ajouter les enjolivages
  mymap <- mymap %>%
    addLayersControl(
      position = "bottomleft",
      baseGroups = c("Toner Lite"),
      overlayGroups  = paste("group ", 1:ncol(belongings), sep = ""),
      options = leaflet::layersControlOptions(collapsed = FALSE))

  for(i in 2:ncol(belongings)){
    mymap <- mymap %>% hideGroup(paste("group ",i,sep=""))
  }

  assign('mymap', mymap, shiny_env)
  assign('mapfun', mapfun, shiny_env)
  assign('spatial4326', spatial4326, shiny_env)
  ##******************************************************************


  ## preparer les violinplots de base *******************************
  group_names <- paste("group ", 1:ncol(belongings))
  best_cat <-group_names[max.col(belongings, ties.method = "first")]
  dataset$tmpgrp <- as.factor(best_cat)

  base_violinplots <- lapply(variables, function(i){

    dataset %>%
      plot_ly(
        x = ~tmpgrp,
        y = dataset[[i]],
        color = ~tmpgrp,
        legendgroup = ~tmpgrp,
        type = 'violin',
        colors = colors,
        box = list(
          visible = T
        ),
        meanline = list(
          visible = T
        )
      ) %>% layout(xaxis = list(title = i), showlegend = FALSE)

  })
  assign('base_violinplots', base_violinplots, shiny_env)
  assign('colors', colors, shiny_env)
  ##******************************************************************

  environment(shiny_ui) <- shiny_env
  environment(shiny_server) <- shiny_env
  app <- shiny::shinyApp(
    ui = shiny_ui,
    server = shiny_server
  )
  shiny::runApp(app, ...)
}


# to add some ellipsis on the plots
# https://www.rdocumentation.org/packages/car/versions/3.0-10/topics/Ellipses


# #library(geocmeans)
# # library(leaflet)
# # library(ggplot2)
# # library(plotly)
# # library(ggsci)
# #
# data(LyonIris)
#
# #selecting the columns for the analysis
# AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14",
#                    "Pct_65","Pct_Img","TxChom1564","Pct_brevet","NivVieMed")
#
# #rescaling the columns
# Data <- LyonIris@data[AnalysisFields]
# for (Col in names(Data)){
#   Data[[Col]] <- as.numeric(scale(Data[[Col]]))
# }
#
# Cmean <- CMeans(Data,4,1.5,500,standardize = FALSE, seed = 456, tol = 0.00001, verbose = FALSE)
#
# sp_clust_explorer(Cmean$Belongings, Data, LyonIris)
#
# # for a chloroplet map : https://rstudio.github.io/leaflet/choropleths.html
