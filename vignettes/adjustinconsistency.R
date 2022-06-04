## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----message=FALSE, warning=FALSE, fig.pos="H", fig.align="center"------------
library(spdep)
library(geocmeans)
library(ggplot2)

#preparing the data
data("LyonIris")

#selecting the columns for the analysis
AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14",
                   "Pct_65","Pct_Img","TxChom1564","Pct_brevet","NivVieMed")

#rescaling the columns
Data <- LyonIris@data[AnalysisFields]
for (Col in names(Data)){
  Data[[Col]] <- scale(Data[[Col]])
}

#creating the queen spatial weight matrix
Neighbours <- poly2nb(LyonIris,queen = TRUE)
WMat <- nb2listw(Neighbours,style="W",zero.policy = TRUE)

#applying SFCM algorithm
SFCM <- SFCMeans(Data, WMat, k = 4, m = 1.5, alpha = 0.7,
                 tol = 0.0001, standardize = FALSE,
                 verbose = FALSE, seed = 456)

#calculating the spatial inconsistency index
consistIndex <- spConsistency(SFCM$Belongings, WMat, nrep = 500)

ggplot() + 
  geom_histogram(aes(x = consistIndex$samples),
                 bins = 30, fill = "white", color = "black") + 
  geom_vline(aes(xintercept = consistIndex$Mean),
             color = "red", linetype="dashed", size = 1) + 
  geom_text(aes(x = consistIndex$Mean+0.0015, y = 43,
                label = round(consistIndex$Mean,2))) + 
  labs(x = "Spatial Inconsistency Index", y = "")

## ----message=FALSE, warning=FALSE, fig.pos="H", fig.align="center"------------
WMat2 <- adjustSpatialWeights(Data, WMat$neighbours, style = "C")
consistIndex2 <- spConsistency(SFCM$Belongings, WMat2, nrep = 500)

ggplot() + 
  geom_histogram(aes(x = consistIndex2$samples),
                 bins = 30, fill = "white", color = "black") + 
  geom_vline(aes(xintercept = consistIndex2$Mean),
             color = "red", linetype="dashed", size = 1) + 
  geom_text(aes(x = consistIndex2$Mean+0.0015, y = 43,
                label = round(consistIndex2$Mean,2))) + 
  labs(x = "Adjusted Spatial Inconsistency Index", y = "")

