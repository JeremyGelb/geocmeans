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

library(shiny)
library(leaflet)
library(plotly)

server <- function(input, output, session) {

  values <- apply(belongings, 1, max) < 0.45

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

  assign('radarchart', radarchart, .GlobalEnv)

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
  assign('base_violinplots', base_violinplots, .GlobalEnv)
  ##******************************************************************


  ## preparer les boxplot de base *******************************
  base_boxplots <- draw_boxplots(dataset, variables, values)
  assign('base_boxplots', base_boxplots, .GlobalEnv)

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

  ## ----------here is an observer when the input of the opacity slider is changed--------------
  observeEvent(input$bg_opacity,{
    op_val <- input$bg_opacity

    # apply to the first map
    leafletProxy('mymap') %>%
      removeTiles("back1") %>%
      removeTiles("back2") %>%
      addProviderTiles(leaflet::providers$Stamen.TonerBackground, group = "Toner Lite", layerId = "back1",
                       options = list(opacity = op_val)) %>%
      addProviderTiles(leaflet::providers$OpenStreetMap, group = "Open Street Map", layerId = "back2",
                       options = list(opacity = op_val))

    leafletProxy('uncertainmap') %>%
      removeTiles("back1") %>%
      removeTiles("back2") %>%
      addProviderTiles(leaflet::providers$Stamen.TonerBackground, group = "Toner Lite", layerId = "back1",
                       options = list(opacity = op_val)) %>%
      addProviderTiles(leaflet::providers$OpenStreetMap, group = "Open Street Map", layerId = "back2",
                       options = list(opacity = op_val))

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