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

## selecting the right number of columns and rows for the violinplots
nv <- length(variables)
nc <- floor(16 / ncol(belongings))
rest <-  nv %% nc
nr <- (nv-rest) / nc
if(rest > 0){
  nr <- nr + 1
}

library(shiny)
library(leaflet)
library(plotly)

welcome_message <- "Welcome in the classification explorer! Click on a feature on the map to start exploring the results of your classification,
The application is designed to help you to investigate the meaning of the groups obtained after a spatial fuzzy classification."

ui <- fluidPage(
    tabsetPanel(
      tabPanel("Interactive map", fluid = TRUE,
               ##------------------- PANNEL 1 : for the interactive map------------------
               titlePanel('Interactive map'),
               # Grid Layout
               {if(is.null(light)){
                 fluidRow(
                   column(width = 4, wellPanel(welcome_message)),
                   column(width = 1, sliderInput("bg_opacity", label = "maps background opacity", min = 0, max = 1, value = 1)),
                   column(width = 7, img(src = "images/logo.png", height = 150, style = "float: right;")),

                 )
               }else{
                 fluidRow(column(width = 4, wellPanel(welcome_message)),
                          column(width = 1,{
                                   if("shinyWidgets" %in% installed.packages()){
                                     shinyWidgets::materialSwitch(inputId = "dark_mode", label = "Dark mode", status = "primary")
                                   }else{
                                     checkboxInput("dark_mode", "Dark mode")
                                   }}
                          ),
                          column(width = 1, sliderInput("bg_opacity", label = "maps background opacity", min = 0, max = 1, value = 1)),
                          column(width = 6, img(src = "images/logo.png", height = 150, style = "float: right;")),
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
