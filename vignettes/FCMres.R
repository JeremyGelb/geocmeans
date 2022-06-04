## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(geocmeans)
library(tmap)
library(dplyr)
library(ggplot2)
library(spdep)
library(raster)

data("LyonIris")

# selecting the columns for the analysis
AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14",
                   "Pct_65","Pct_Img","TxChom1564","Pct_brevet","NivVieMed")

# rescaling the columns
Data <- LyonIris@data[AnalysisFields]
for (Col in names(Data)){
  Data[[Col]] <- scale(Data[[Col]])
}

# applying the hclust function
clust <- hclust(dist(Data), method = "ward")

# getting the groups
LyonIris$Hclust_groups <- as.character(cutree(clust, k = 4))
Data$Hclust_groups <- as.character(cutree(clust, k = 4))

# mapping the groups
tm_shape(LyonIris) + 
  tm_polygons(col = "Hclust_groups", title = "groups")


## ----message=FALSE, warning=FALSE---------------------------------------------
centers <- Data %>% 
  group_by(Hclust_groups) %>%
  summarise_all(mean)

centers <- as.matrix(centers[2:ncol(centers)])

## ----message=FALSE, warning=FALSE---------------------------------------------
member_mat <- cat_to_belongings(Data$Hclust_groups)

## ----message=FALSE, warning=FALSE---------------------------------------------
Data$Hclust_groups <- NULL

hclustres <- FCMres(list(
  "Centers" = centers,
  "Belongings" = member_mat,
  "Data" = Data,
  "m" = 1,
  "algo" = "hclust"
))

## ----message=FALSE, warning=FALSE, eval = FALSE-------------------------------
#  # quick summaries about the groups
#  summary(hclustres)
#  violinPlots(hclustres$Data, hclustres$Groups)
#  spiderPlots(hclustres$Data, hclustres$Belongings)
#  mapClusters(LyonIris, hclustres)
#  
#  # some indices about classification quality
#  calcqualityIndexes(hclustres$Data,
#                     hclustres$Belongings,
#                     hclustres$m)
#  
#  # spatial diagnostic
#  Neighbours <- poly2nb(LyonIris,queen = TRUE)
#  WMat <- nb2listw(Neighbours,style="W",zero.policy = TRUE)
#  spatialDiag(hclustres, nblistw = WMat)
#  
#  # investigation with the shiny app
#  sp_clust_explorer(hclustres, spatial = LyonIris)

## ----message=FALSE, warning=FALSE---------------------------------------------
data("Arcachon")

# loading each raster as a column in a matrix
# and scale each column
all_data <- do.call(cbind, lapply(names(Arcachon), function(n){
  rast <- Arcachon[[n]]
  return(raster::values(raster::scale(rast)))
}))

# removing the rows with missing values
missing <- complete.cases(all_data)
all_data <- all_data[missing,]

# applying the kmeans algorithm with 7 groups
kmean7 <- kmeans(all_data, 7)


## ----message=FALSE, warning=FALSE---------------------------------------------

# creating Data (do not forget the standardization)
Data <- lapply(names(Arcachon), function(n){
  rast <- Arcachon[[n]]
  return(raster::scale(rast))
})
names(Data) <- names(Arcachon)

# creating rasters
ref_raster <- Arcachon[[1]]

rasters <- lapply(1:7, function(i){
  # creating a vector with only 0 values
  vals <- rep(0, ncell(ref_raster))
  # filling it with values when the pixels are not NA
  vals[missing] <- ifelse(kmean7$cluster == i,1,0)
  # setting the values in a rasterLayer
  rast <- ref_raster
  raster::values(rast) <- vals
  return(rast)
})

# creating centers
all_data <- as.data.frame(all_data)
names(all_data) <- names(Arcachon)
all_data$kmean_groups <- as.character(kmean7$cluster)

centers <- all_data %>% 
  group_by(kmean_groups) %>%
  summarise_all(mean)

centers <- as.matrix(centers[2:ncol(centers)])

## ----message=FALSE, warning=FALSE---------------------------------------------
myFCMres <- FCMres(list(
  "Data" = Data,
  "Centers" = centers,
  "rasters" = rasters,
  "m" = 1,
  "algo" = "kmeans"
))

## ----message=FALSE, warning=FALSE, eval = FALSE-------------------------------
#  # quick summaries about the groups
#  summary(myFCMres)
#  violinPlots(myFCMres$Data, myFCMres$Groups)
#  spiderPlots(myFCMres$Data, myFCMres$Belongings)
#  mapClusters(object = myFCMres)
#  
#  # some indices about classification quality
#  calcqualityIndexes(myFCMres$Data,
#                     myFCMres$Belongings,
#                     myFCMres$m)
#  
#  # spatial diagnostic
#  w1 <- matrix(1, nrow = 3, ncol = 3)
#  spatialDiag(myFCMres, window = w1, nrep = 5)
#  
#  # investigation with the shiny app
#  sp_clust_explorer(myFCMres)

