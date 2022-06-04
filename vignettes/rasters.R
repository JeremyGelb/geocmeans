## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----include=FALSE------------------------------------------------------------
#loading the pre-calculated results
load(system.file("extdata", "results_vignette_raster.rda",
                           package = "geocmeans", mustWork = TRUE))

## ----message=FALSE, warning=FALSE, fig.width= 4-------------------------------
library(raster)
library(geocmeans)
library(ggpubr)
library(future)
library(tmap)
library(viridis)
library(RColorBrewer)

data("Arcachon")

# show the pseudo-color image
plotRGB(Arcachon, r = 3, g = 2, b = 1, stretch = "hist")

## ----message=FALSE, warning=FALSE---------------------------------------------
# sonverting the RasterBrick to a simple list of RasterLayer
dataset <- lapply(names(Arcachon), function(n){
  aband <- Arcachon[[n]]
  return(aband)
})

# giving a name to each band
names(dataset) <- names(Arcachon)

## ----message=FALSE, warning=FALSE, eval=FALSE---------------------------------
#  # finding an appropriate k and m values (using a multicore plan)
#  future::plan(future::multisession(workers = 6))
#  FCMvalues <- select_parameters.mc(algo = "FCM", data = dataset,
#                                 k = 5:10, m = seq(1.1,2,0.1), spconsist = FALSE,
#                                 indices = c("XieBeni.index", "Explained.inertia",
#                                             "Negentropy.index", "Silhouette.index"),
#                                 verbose = TRUE)

## ----message=FALSE, warning=FALSE, fig.width = 5, fig.align='center'----------
# plotting the silhouette index values
ggplot(FCMvalues) + 
  geom_raster(aes(x = m, y = k, fill = Silhouette.index)) + 
  geom_text(aes(x = m, y = k, label = round(Silhouette.index,2)), size = 2)+
  scale_fill_viridis() +
  coord_fixed(ratio=0.125) 

# plotting the explained inertia
ggplot(FCMvalues) + 
  geom_raster(aes(x = m, y = k, fill = Explained.inertia)) + 
  geom_text(aes(x = m, y = k, label = round(Explained.inertia,2)), size = 2)+
  scale_fill_viridis() +
  coord_fixed(ratio=0.125)


