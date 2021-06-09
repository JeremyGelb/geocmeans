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

#' @title Draw the scatter plot of pannel 2 and 3
#'
#' @description Draw the scatter plot of pannel 2 and 3
#' @param x The values for the x axis
#' @param y The values for the y axis
#' @param col_values The values for the color
#' @param colors A vector with the colors to use
#' @param gpvol The index of the color to use
#' @param belongings the membership matrix
#' @param input the shiny app input
#' @param qual A boolean indicating if the plot must be drawn for pannel 2
#' (FALSE) or 3 (TRUE)
#' @param qual_colors The colors to use if qual = TRUE
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
draw_byplot <- function(x, y, col_values, colors, gpcol, belongings, input, qual = FALSE, qual_colors = NULL){

  if (qual == FALSE){
    idx <- order(col_values)
    tx <-  x[idx]
    ty <-  y[idx]
    tcol <- col_values[idx]
    ramp <- colorRampPalette(c("white", colors[[gpcol]]))(10)
  }else{
    tx <- x
    ty <- y
    tcol <- col_values
    ramp <- qual_colors
  }

  idx <- order(belongings[,gpcol])
  biplot <- plot_ly(
    x = tx,
    y = ty,
    color = tcol,
    colors = ramp,
    type = "scatter",
    mode = "markers",
    size = 2
  )

  if(isTRUE(input$show_ellipsis) & qual == FALSE){
    for(j in 1:ncol(belongings)){
      coords <- car::dataEllipse(x,
                                 y,
                                 weights = belongings[,j],
                                 levels = 0.75,
                                 draw = FALSE,
      )
      coords <- coords[1:(nrow(coords)-1),]
      coords <- rbind(coords,coords[1,])
      biplot <- biplot %>%
        plotly::add_paths(
          x = coords[,1],
          y = coords[,2],
          line = list(width = 2),
          color = I(colors[[j]]),

        )
    }
  }

  return(biplot)
}

#' @title Draw the box plots of pannel 3
#'
#' @description Draw the box plots of pannel 3
#' @param dataset The dataset used
#' @param variables The variables names used for clustering
#' @param values A boolean vector indicating which observations must
#' be put in red on the chart
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
draw_boxplots <- function(dataset, variables, values){

  new_boxplots <- lapply(variables, function(i){
    dataset$myx <- 0
    df1 <- dataset[!values,]
    df2 <- dataset[values,]

    bxplot <- dataset %>%
      plot_ly(
        x = dataset$myx,
        y = dataset[[i]],
        type = 'box',
        boxpoints = FALSE
      ) %>%
      layout(xaxis = list(title = i), showlegend = FALSE)

    if (nrow(df1) > 0){
      bxplot <- bxplot %>%
        add_markers(
          x = ~jitter(df1$myx, factor = 5),
          y = df1[[i]],
          marker = list(size = 4,
                        color = "grey")
        )
    }

    if(nrow(df2) > 0){
      bxplot <- bxplot %>%
        add_markers(
          x = ~jitter(df2$myx, factor = 5),
          y = df2[[i]],
          marker = list(size = 8,
                        color = "red")
        )
    }

    return(bxplot)
  })

  return(new_boxplots)
}


#' @title Adjust the background of plot
#'
#' @description Adjust the background of plot (light or dark mode)
#' @param plot The plot (plotly) to adjust
#' @param ligght The actual theme used
#' @param input The shiny app input object
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
adj_bg_color <- function(plot, light, input){
  if(is.null(light) == FALSE){
    if (isTRUE(input$dark_mode)){
      plot <- plot %>% layout(
        font = list(color = "white"),
        paper_bgcolor = "#303030",
        plot_bgcolor = "#303030")
    }
  }
  plot
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
#' @importFrom shiny reactive observeEvent observe
#' @importFrom leaflet renderLeaflet leafletProxy removeShape clearGroup
#' @importFrom plotly renderPlotly plot_ly layout add_paths
#' @importFrom grDevices colorRampPalette
#'
#' @keywords interal
#' @examples
#' #This is an internal function, no example provided
shiny_server <- function(input, output, session) {

  ## creating the map man
  output$mymap <- renderLeaflet({
    mymap
  })

  ## and the map on the third pannel
  output$uncertainmap <- renderLeaflet({
    uncertainMap
  })

  ## and the first bivariate plot
  bivar_params <- reactive({
    list(x = input$var1_biplot,
         y = input$var2_biplot,
         color = input$group_biplot)
  })

  ## and the radarchart
  output$radarchart <- renderPlotly({
    radarchart
  })

  output$bivar_plot <- renderPlotly({
    params <- bivar_params()
    gpcol <- as.numeric(strsplit(params$color," ", fixed = TRUE)[[1]][[2]])

    biplot <- draw_byplot(
      x = dataset[[params$x]],
      y = dataset[[params$y]],
      col_values = belongings[,gpcol],
      colors = colors,
      gpcol = gpcol,
      belongings = belongings,
      input = input
    )

    # adjusting the color with the theme
    adj_bg_color(biplot, light, input)
  })

  ## and the second bivariate plot
  bivar_params2 <- reactive({
    list(x = input$var1_biplot2,
         y = input$var2_biplot2,
         proba = input$uncertain1)
  })


  output$bivar_plot2 <- renderPlotly({
    params <- bivar_params2()
    test <- apply(belongings, 1, max) < params$proba
    col_values <- factor(ifelse(test, "uncertain", "classified"))
    if(length(unique(col_values)) == 1){
      qual_colors <- "grey"
    }else{
      qual_colors <- c("grey","red")
    }

    biplot <- draw_byplot(
      x = dataset[[params$x]][order(col_values)],
      y = dataset[[params$y]][order(col_values)],
      col_values = col_values[order(col_values)],
      colors = colors,
      gpcol = 1,
      belongings = belongings,
      input = input,
      qual = TRUE,
      qual_colors = qual_colors
    )

    # adjusting the color with the theme
    adj_bg_color(biplot, light, input)
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

  ## putting the original boxplots
  lapply(1:length(base_boxplots), function(i){
    vplot <- base_boxplots[[i]]
    output[[paste("boxplots",i,sep="")]] <- renderPlotly({vplot})
  })



  firsttime <- TRUE
  ## ----------here is an observer working when we click on the map of the first pannel--------------
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

      violin2 <-
        adj_bg_color(violin, light, input) %>%
        layout(shapes = list(hline(value, color = "red")))
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
                    layerId = "highlighter", fillColor = NULL, group = "")
      firsttime <- FALSE
    }else{
      leafletProxy('mymap') %>%
        removeShape("highlighter") %>%
        mapfun(data = feat, weight = 2, opacity = 1.0, fillOpacity = 0, color = "red",
               fillColor = NULL, layerId = "highlighter", group = "")
    }


  })

  ## ----------here is an observer listening for the ellipsis--------------
  observeEvent(input$show_ellipsis,{

    # redraw the biplot
    output$bivar_plot <- renderPlotly({
      params <- bivar_params()
      gpcol <- as.numeric(strsplit(params$color," ", fixed = TRUE)[[1]][[2]])

      biplot <- draw_byplot(
        x = dataset[[params$x]],
        y = dataset[[params$y]],
        col_values = belongings[,gpcol],
        colors = colors,
        gpcol = gpcol,
        belongings = belongings,
        input = input
      )

      # adjusting the color with the theme
      biplot <- adj_bg_color(biplot, light, input)

      biplot
    })

  })

  ## ----------here is an observer working when we set the slider of the third pannel--------------
  observeEvent(input$uncertain1,{

    # we have to redraw the second layer of this map
    # step1 : selecting the appropriate new features
    tol <- input$uncertain1
    values <- apply(belongings, 1, max) < tol
    spdf <- subset(spatial4326, values)

    ## Step2 : remove the previous layer and add the new one
    leafletProxy('uncertainmap') %>%
      clearGroup(group = "binaryUncertain") %>%
      mapfun(data = spdf,
             weight = 1,
             group = "binaryUncertain",
             color = "black",
             fillColor = "red",
             fillOpacity = 0.7,
             layerId = 1:nrow(spdf))

    ## we also have to redraw the box plots
    new_boxplots <- draw_boxplots(dataset, variables, values)

    lapply(1:length(new_boxplots), function(i){
      vplot <- adj_bg_color(new_boxplots[[i]], light, input)
      output[[paste("boxplots",i,sep="")]] <- renderPlotly({vplot})
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

        biplot <- draw_byplot(
          x = dataset[[params$x]],
          y = dataset[[params$y]],
          col_values = belongings[,gpcol],
          colors = colors,
          gpcol = gpcol,
          belongings = belongings,
          input = input
        )

        # adjusting the color with the theme
        biplot <- adj_bg_color(biplot, light, input)

        biplot
      })

      # set the right colors for violin plots
      new_violins <- lapply(1:length(variables), function(i){
        violin <- base_violinplots[[i]]
        violin2 <- adj_bg_color(violin, light, input)
        output[[paste("violinplots",i,sep="")]] <- renderPlotly({violin2})

      })

      # set the right colors for boxplots plots
      tol <- input$uncertain1
      values <- apply(belongings, 1, max) < tol

      new_boxplots <- draw_boxplots(dataset, variables, values)

      lapply(1:length(new_boxplots), function(i){
        vplot <- adj_bg_color(new_boxplots[[i]], light, input)
        output[[paste("boxplots",i,sep="")]] <- renderPlotly({vplot})
      })


      # set the right colors for the radar
      output$radarchart <- renderPlotly({
        adj_bg_color(radarchart, light, input)
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
#' @keywords interal
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


#' @title Set the box plot pannel
#'
#' @description A simple function to prepare the box plot pannel in uncertainty
#' analysis
#' @param n The number of plots
#' @keywords interal
#' @examples
#' #This is an internal function, no example provided
boxplots_ui <- function(n){
  nc <- 6
  nr <- ceiling(n / nc)
  # NOTE : each plot will have a width of 2
  fulltxt <- "wellPanel("
  cnt <- 1
  for(i in 1:nr){
    new_row <-"fluidRow("
    for(j in 1:nc){
      new_col <- paste("column(plotlyOutput('boxplots",cnt,"'), width = 2),",sep="")
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
#' @importFrom shiny fluidPage tabsetPanel tabPanel titlePanel fluidRow wellPanel column selectInput checkboxInput sliderInput
#' @importFrom leaflet leafletOutput
#' @importFrom plotly plotlyOutput
#' @importFrom utils installed.packages
#' @keywords interal
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
                   {
                     if("shinyWidgets" %in% installed.packages()){
                       column(width = 2, shinyWidgets::materialSwitch(inputId = "dark_mode", label = "Dark mode", status = "primary"))
                     }else{
                       column(width = 2, checkboxInput("dark_mode", "Dark mode"))
                     }
                   }
                   )
        }},
        fluidRow(column(width = 8, wellPanel(leafletOutput("mymap"))),
                 column(width = 4,
                        tabsetPanel(
                        tabPanel("membership values", fluid = TRUE,
                                 wellPanel(plotlyOutput("barplot1", height = "360px"))
                                 ),
                        tabPanel("radar chart", fluid = TRUE,
                                 wellPanel(plotlyOutput("radarchart", height = "360px"))
                        ),
                        tabPanel("general informations", fluid = TRUE,
                                 wellPanel("some information about your clustering")
                        )
                        ))
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
                                   if("shinyWidgets" %in% installed.packages()){
                                     shinyWidgets::materialSwitch(inputId = "show_ellipsis", label = "show ellipsis", status = "primary")
                                   }else{
                                     checkboxInput("show_ellipsis", "show ellipsis")
                                   }
                                 }
                               },
                               ),
                        column(width = 10,plotlyOutput("bivar_plot", height = "800px")),
               ),
      ),

      tabPanel("Uncertain observations", fluid = TRUE,
               ##------------------- PANNEL 3 : Uncertainty analysis------------------
               titlePanel('Analysis of uncertain observations'),

               # Grid Layout
               fluidRow(wellPanel("In this panel, you can explore the observations that are not well classified")),
               fluidRow(column(width = 2,sliderInput("uncertain1", "minimum probability", min = 0, max = 1, value = 0.45))),
               fluidRow(column(width = 5, leafletOutput("uncertainmap", height = "600px")),
                        column(width = 1,
                               selectInput("var1_biplot2", "X axis variable", variables, selected = variables[[1]]),
                               selectInput("var2_biplot2", "Y axis variable", variables, selected = variables[[2]]),
                               ),
                        column(width = 6, plotlyOutput("bivar_plot2", height = "600px")),
                        ),
               fluidRow(boxplots_ui(nv)),
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
#' @param ... Other parameters passed to the function runApp
#' @param spatial A spatial object (SpatialPointsDataFrame, SpatialPolygonsDataFrame or
#' SpatialLinesDataFrame) used to map the observations
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
sp_clust_explorer <- function(belongings, dataset, spatial, ...) {

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

  ## preparer la carte leaflet du premier panneau ***************************************

  mymap <- leaflet(height = "600px") %>%
    addProviderTiles(leaflet::providers$Stamen.TonerBackground, group = "Toner Lite")

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
    #spatial <- rgeos::gBuffer(spatial, width = 0.01, byid = TRUE)
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
      baseGroups = c("Toner Lite"),
      overlayGroups  = c(paste("group ", 1:ncol(belongings), sep = ""),"Most likely group"),
      options = leaflet::layersControlOptions(collapsed = FALSE))

  for(i in 2:ncol(belongings)){
    mymap <- mymap %>% hideGroup(paste("group ",i,sep=""))
  }
  mymap <- mymap %>% hideGroup("Most likely group")

  assign('mymap', mymap, shiny_env)
  assign('mapfun', mapfun, shiny_env)
  assign('spatial4326', spatial4326, shiny_env)

  ## preparer la carte leaflet du troisieme panneau ***************************************
  uncertainMap <- leaflet(height = "600px") %>%
    addProviderTiles(leaflet::providers$Stamen.TonerBackground, group = "Toner Lite")

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
      baseGroups = c("Toner Lite"),
      overlayGroups  = c("UncertaintyIdx","binaryUncertain"),
      options = leaflet::layersControlOptions(collapsed = FALSE)) %>%
    hideGroup("UncertaintyIdx")

  assign('uncertainMap', uncertainMap, shiny_env)
  ##******************************************************************

  ## preparer le radar chart *****************************************

  # step1 : calculating the values (min and max normalisation)
  clustmeans <- apply(belongings, 2, function(w){
    means <- apply(dataset, 2, function(d){
      sum(d*w) / sum(w)
    })
  })

  sc_clustmeans <- do.call(rbind,lapply(1:nrow(clustmeans), function(i){
    x <- clustmeans[i,]
    return((x-min(x))/(max(x)-min(x)))
  }))
  rownames(sc_clustmeans) <- rownames(clustmeans)
  colnames(sc_clustmeans) <- colnames(clustmeans)

  # step2 : drawing the radarchart
  radarchart <- plot_ly(
    type = 'scatterpolar',
    fill = 'toself'
  )
  for(i in 1:ncol(sc_clustmeans)){
    radarchart <- radarchart %>%
      add_trace(
        r = sc_clustmeans[,i],
        theta = rownames(sc_clustmeans),
        name = paste('Group ', i, sep=""),
        fillcolor = colors[[i]],
        opacity = 0.4,
        marker=list(color = colors[[i]])
      )
  }

  radarchart <- radarchart %>%
    layout(
      polar = list(
        radialaxis = list(
          visible = T,
          range = c(0,1)
        )
      )
    )

  assign('radarchart', radarchart, shiny_env)


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


  ## preparer les boxplot de base *******************************
  base_boxplots <- draw_boxplots(dataset, variables, values)
  assign('base_boxplots', base_boxplots, shiny_env)

  ##******************************************************************
  environment(shiny_ui) <- shiny_env
  environment(shiny_server) <- shiny_env
  app <- shiny::shinyApp(
    ui = shiny_ui,
    server = shiny_server
  )
  shiny::runApp(app, ...)
}
