## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----message=FALSE, warning=FALSE---------------------------------------------
#charging packages and data
library(geocmeans)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(viridis)
library(spdep)

data(LyonIris)

# selecting the columns for the analysis
AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14",
                   "Pct_65","Pct_Img","TxChom1564","Pct_brevet","NivVieMed")

# rescaling the columns
Data <- LyonIris@data[AnalysisFields]
for (Col in names(Data)){
  Data[[Col]] <- scale(Data[[Col]])
}

# preparing some elements for further mapping
LyonIris$OID <- as.character(1:nrow(LyonIris))
FortiData <- ggplot2::fortify(LyonIris,region="OID")

## ----include=FALSE------------------------------------------------------------
# loading the pre-calculated results
load(system.file("extdata", "results_vignette_intro.rda",
                           package = "geocmeans", mustWork = TRUE))

## ----warning=FALSE, fig.cap = "Impact of the number of groups on the explained variance", out.width = "50%", fig.pos="H", fig.align="center"----
# finding the best k by using the r2 of the classification
# trying for k from 2 to 10
R2s <- sapply(2:10,function(k){
  Clust <- kmeans(Data,centers=k,iter.max = 150)
  R2 <- Clust$betweenss / Clust$totss
  return(R2)
})

Df <- data.frame(K=2:10,
                 R2 = R2s)

ggplot(Df)+
  geom_line(aes(x=K,y=R2s))+
  geom_point(aes(x=K,y=R2s),color="red")+
  xlab("Number of groups")+
  ylab("R2 of classification")

## ----warning=FALSE, fig.cap="Clusters identified by classical k-means", fig.pos="H", fig.align="center"----
KMeanClust <-  kmeans(Data,centers=4,iter.max = 150)
LyonIris$Cluster <-paste("cluster",KMeanClust$cluster,sep="_")

# mapping the groups
DFmapping <- merge(FortiData,LyonIris,by.x="id",by.y="OID")

ggplot(data=DFmapping)+
  geom_polygon(aes(x=long,y=lat,group=group,fill=Cluster),color=rgb(0,0,0,0))+
  coord_fixed(ratio = 1)+
  scale_fill_manual(name="Cluster",values = c("cluster_1"="palegreen3",
                                 "cluster_2"="firebrick",
                                 "cluster_3"="lightyellow2",
                                 "cluster_4"="steelblue"))+
  theme( axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank()
    )


## ----warning=FALSE------------------------------------------------------------
Cmean <- CMeans(Data,4,1.5,500,standardize = FALSE, seed = 456, tol = 0.00001, verbose = FALSE)

## ----message=FALSE, warning=FALSE---------------------------------------------
calcqualityIndexes(Data, Cmean$Belongings, m = 1.5)

## ----warning=FALSE, eval = FALSE----------------------------------------------
#  beta_values <- selectParameters("GFCM",data = Data, k = 4, m = 1.5,
#                                   beta = seq(0,1,0.05), spconsist = FALSE,
#                                   tol = 0.00001, seed = 456)

## ----warning=FALSE------------------------------------------------------------
knitr::kable(beta_values[c("beta","Silhouette.index","XieBeni.index","Explained.inertia")],
             col.names = c("beta", "silhouette index",
                           "Xie and Beni index", "explained inertia"),digits = 3)

## ----warning=FALSE------------------------------------------------------------
GCmean <- GCMeans(Data,k = 4,m = 1.5, beta = 0.7,500,standardize = FALSE, seed=456,
                  tol = 0.00001, verbose = FALSE)
r1 <- calcqualityIndexes(Data,GCmean$Belongings,m=1.5)
r2 <- calcqualityIndexes(Data,Cmean$Belongings,m=1.5)
df <- cbind(unlist(r1), unlist(r2))

knitr::kable(df,
             digits = 3,col.names = c("GFCM", "FCM"))

## ----warning=FALSE, fig.cap = "Probability of belonging to cluster 1", out.width="80%", fig.pos="H", fig.align="center"----
cmeansMaps<- mapClusters(LyonIris,Cmean$Belongings,undecided = 0.45)
GcmeansMaps<- mapClusters(LyonIris,GCmean$Belongings,undecided = 0.45)

ggarrange(cmeansMaps$ProbaMaps[[1]],GcmeansMaps$ProbaMaps[[1]], 
          nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom")

## ----warning=FALSE, fig.cap = "Probability of belonging to cluster 2", out.width="80%", fig.pos="H", fig.align="center"----
ggarrange(cmeansMaps$ProbaMaps[[2]],GcmeansMaps$ProbaMaps[[2]], 
          nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom")

## ----warning=FALSE, fig.cap = "Probability of belonging to cluster 3", out.width="80%", fig.pos="H", fig.align="center"----
ggarrange(cmeansMaps$ProbaMaps[[3]],GcmeansMaps$ProbaMaps[[3]], 
          nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom")

## ----warning=FALSE, fig.cap = "Probability of belonging to cluster 4", out.width="80%", fig.pos="H", fig.align="center"----
ggarrange(cmeansMaps$ProbaMaps[[4]],GcmeansMaps$ProbaMaps[[4]], 
          nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom")

## ----warning=FALSE, fig.cap = "Most likely clusters and undecided units", out.width="80%", fig.pos="H", fig.align="center"----
ggarrange(cmeansMaps$ClusterPlot,GcmeansMaps$ClusterPlot,
          nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom")

## ----message=FALSE, warning=FALSE---------------------------------------------
library(spdep)

Neighbours <- poly2nb(LyonIris,queen = TRUE)
WMat <- nb2listw(Neighbours,style="W",zero.policy = TRUE)

## ----warning=FALSE,  message=FALSE, eval = FALSE------------------------------
#  DFindices_SFCM <- selectParameters(algo = "SFCM", data = Data,
#                                 k = 4, m = 1.5, alpha = seq(0,2,0.05),
#                                 nblistw = WMat, standardize = FALSE,
#                                 tol = 0.0001, verbose = FALSE, seed = 456)
#  

## ----warning=FALSE, fig.cap = "Link between alpha and spatial inconsistency", out.width="50%", fig.pos="H", fig.align="center"----
ggplot(DFindices_SFCM)+
  geom_smooth(aes(x=alpha,y=spConsistency), color = "black")+
  geom_point(aes(x=alpha,y=spConsistency), color = "red")


## ----warning=FALSE, fig.cap = "Link between alpha and explained inertia", out.width="50%", fig.pos="H", fig.align="center"----
ggplot(DFindices_SFCM)+
  geom_smooth(aes(x=alpha,y=Explained.inertia), color = "black")+
  geom_point(aes(x=alpha,y=Explained.inertia), color = "red")

## ----fig.cap="Link between alpha and silhouette index", message=FALSE, warning=FALSE, fig.pos="H", fig.align="center", out.width="50%"----
ggplot(DFindices_SFCM)+
  geom_smooth(aes(x=alpha,y=Silhouette.index), color = "black")+
  geom_point(aes(x=alpha,y=Silhouette.index), color = "red")

## ----fig.cap="Link between alpha and Xie and Beni index", message=FALSE, warning=FALSE, out.width="50%", fig.pos="H", fig.align="center"----
ggplot(DFindices_SFCM)+
  geom_smooth(aes(x=alpha,y=XieBeni.index), color = "black")+
  geom_point(aes(x=alpha,y=XieBeni.index), color = "red")


## ----warning=FALSE------------------------------------------------------------
SFCM <- SFCMeans(Data, WMat, k = 4, m = 1.5, alpha = 0.7,
                 tol = 0.0001, standardize = FALSE,
                 verbose = FALSE, seed = 456)

## ----warning=FALSE, eval = FALSE----------------------------------------------
#  future::plan(future::multisession(workers=4))
#  DFindices_SFGCM <- selectParameters.mc(algo = "SGFCM", data = Data,
#                                 k = 4, m = 1.5, alpha = seq(0,2,0.05),
#                                 beta = seq(0,0.85,0.05),
#                                 nblistw = WMat, standardize = FALSE, chunk_size = 50,
#                                 tol = 0.0001, verbose = FALSE, seed = 456)

## ----warning=FALSE, fig.cap = "Impact of beta and alpha on silhouette index", out.width = "80%", fig.pos="H", fig.align="center"----
ggplot(DFindices_SFGCM) + 
  geom_raster(aes(x = alpha, y = beta, fill = Silhouette.index),  size = 5) + 
  scale_fill_viridis() +
  coord_fixed(ratio=1)

## ----warning=FALSE, fig.cap = "Impact of beta and alpha on Xie and Beni index", out.width = "80%", fig.pos="H", fig.align="center"----
ggplot(DFindices_SFGCM) + 
  geom_raster(aes(x = alpha, y = beta, fill = XieBeni.index),  size = 5) + 
  scale_fill_viridis() +
  coord_fixed(ratio=1)

## ----warning=FALSE, fig.cap = "Impact of beta and alpha on spatial inconsistency", out.width = "80%", fig.pos="H", fig.align="center"----
ggplot(DFindices_SFGCM) + 
  geom_raster(aes(x = alpha, y = beta, fill = spConsistency),  size = 5) + 
  scale_fill_viridis() +
  coord_fixed(ratio=1)


## ----warning=FALSE------------------------------------------------------------
SGFCM <- SGFCMeans(Data,WMat,k = 4,m=1.5, alpha=0.95, beta = 0.65,
                    tol=0.0001, standardize = FALSE, verbose = FALSE, seed = 456)

## ----warning=FALSE------------------------------------------------------------
r1 <- calcqualityIndexes(Data, SFCM$Belongings,m = 1.5)
r2 <- calcqualityIndexes(Data, SGFCM$Belongings,m = 1.5)

diagSFCM <- spatialDiag(SFCM$Belongings, nblistw = WMat,
                        undecided = 0.45,nrep = 500)
diagSGFCM <- spatialDiag(SGFCM$Belongings, nblistw = WMat,
                         undecided = 0.45,nrep = 500)

df <- cbind(
  c(unlist(r1),diagSFCM$SpConsist),
  c(unlist(r2),diagSGFCM$SpConsist)
)
row.names(df)[length(row.names(df))] <- "sp.consistency"

knitr::kable(df,digits = 3,col.names = c("SFCM","SGFCM"))

## ----warning=FALSE, fig.cap = "Probability of belonging to cluster 1", out.width="80%", fig.pos="H", fig.align="center"----
SFCMMaps <- mapClusters(geodata = LyonIris, object = SFCM$Belongings,undecided = 0.45)
SGFCMMaps <- mapClusters(geodata = LyonIris, object = SGFCM$Belongings,undecided = 0.45)

ggarrange(SFCMMaps$ProbaMaps[[1]],SGFCMMaps$ProbaMaps[[1]], nrow = 1, ncol = 2,
          common.legend = TRUE, legend = "bottom")

## ----warning=FALSE, fig.cap = "Probability of belonging to cluster 2", out.width="80%", fig.pos="H", fig.align="center"----
ggarrange(SFCMMaps$ProbaMaps[[2]],SGFCMMaps$ProbaMaps[[2]], nrow = 1, ncol = 2,
          common.legend = TRUE, legend = "bottom")

## ----warning=FALSE, fig.cap = "Probability of belonging to cluster 3", out.width="80%", fig.pos="H", fig.align="center"----
ggarrange(SFCMMaps$ProbaMaps[[3]],SGFCMMaps$ProbaMaps[[3]], nrow = 1, ncol = 2,
          common.legend = TRUE, legend = "bottom")

## ----warning=FALSE, fig.cap = "Probability of belonging to cluster 4", out.width="80%"----
ggarrange(SFCMMaps$ProbaMaps[[4]],SGFCMMaps$ProbaMaps[[4]], nrow = 1, ncol = 2,
          common.legend = TRUE, legend = "bottom")

## ----warning=FALSE, fig.cap = "Most likely cluster and undecided units", out.width="80%", fig.pos="H", fig.align="center"----
ggarrange(SFCMMaps$ClusterPlot,SGFCMMaps$ClusterPlot, nrow = 1, ncol = 2,
          common.legend = TRUE, legend = "bottom")

## ----warning=FALSE------------------------------------------------------------
spdiag_1 <- spatialDiag(Cmean$Belongings, nblistw = WMat, nrep=250)
spdiag_2 <- spatialDiag(GCmean$Belongings, nblistw = WMat, nrep=250)
spdiag_3 <- spatialDiag(SFCM$Belongings, nblistw = WMat, nrep=250)
spdiag_4 <- spatialDiag(SGFCM$Belongings, nblistw = WMat, nrep=250)

#looking at the moran I values for each group
moran_table <- data.frame(cbind(spdiag_1$MoranValues$MoranI,
                     spdiag_2$MoranValues$MoranI,
                     spdiag_3$MoranValues$MoranI,
                     spdiag_4$MoranValues$MoranI
                     ))
row.names(moran_table) <- paste("cluster ",1:4,sep="")
knitr::kable(moran_table, digits = 3,
             col.names = c("FCM","GFCM","SFCM","SGFCM"),
             caption = "Moran I index for the columns of the membership matrix"
             )

## ----warning=FALSE------------------------------------------------------------
print(c(spdiag_1$SpConsist, spdiag_2$SpConsist,spdiag_3$SpConsist,spdiag_4$SpConsist))

## ----warning=FALSE------------------------------------------------------------
sum(spdiag_4$SpConsist > spdiag_3$SpConsistSamples) / length(spdiag_3$SpConsistSamples)

## ----warning=FALSE------------------------------------------------------------
Undecided <- undecidedUnits(SGFCM$Belongings,0.45)
LyonIris$FinalCluster <- ifelse(Undecided=="Undecided",
                                "Undecided",paste("cluster",Undecided,sep="_"))

# mapping the groups
DFmapping <- merge(FortiData,LyonIris,by.x="id",by.y="OID")

ggplot(data=DFmapping)+
  geom_polygon(aes(x=long,y=lat,group=group,fill=FinalCluster),color=rgb(0,0,0,0))+
  coord_fixed(ratio = 1)+
  scale_fill_manual(name="FinalCluster",values = c("cluster_V1"="palegreen3",
                                 "cluster_V2"="firebrick",
                                 "cluster_V3"="lightyellow2",
                                 "cluster_V4"="steelblue",
                                 "cluster_V5"="pink",
                                 "Undecided"=rgb(0,0,0,0.4)))+
  theme( axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
    )

## ----warning=FALSE------------------------------------------------------------
colors <- c("palegreen3","firebrick","lightyellow2","steelblue","pink")
uncertaintyMap(LyonIris, SGFCM$Belongings, color = colors)

## ----warning=FALSE------------------------------------------------------------
LyonIris$entropyidx  <- calcUncertaintyIndex(SGFCM$Belongings)

# mapping the uncertainty
DFmapping <- merge(FortiData,LyonIris,by.x="id",by.y="OID")

ggplot(data=DFmapping)+
  geom_polygon(aes(x=long,y=lat,group=group,fill=entropyidx),color=rgb(0,0,0,0))+
  coord_fixed(ratio = 1)+
  labs(title = "Uncertainty evaluation", fill = "entropy index") +
  theme(axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
    )

## ----warning=FALSE------------------------------------------------------------
summarizeClusters(LyonIris@data[AnalysisFields],belongmatrix = SGFCM$Belongings,
                  weighted = TRUE, dec = 3)
# equivalent to : 
# summary(SGFCM, LyonIris@data[AnalysisFields])

## ----warning=FALSE, eval = FALSE----------------------------------------------
#  spiderPlots(LyonIris@data[AnalysisFields], SGFCM$Belongings,
#              chartcolors = c("darkorange3","grey4","darkgreen","royalblue"))
#  violinPlots(LyonIris@data[AnalysisFields], SGFCM$Groups)

## ----warning=FALSE, eval = FALSE----------------------------------------------
#  bootvalues <- boot_group_validation(SGFCM, nsim = 1000, maxiter = 1000,
#                                       tol = 0.0001, verbose = FALSE)

## ----warning=FALSE------------------------------------------------------------
melted_df <- reshape2::melt(bootvalues$group_consistency)
melted_df$variable <- as.factor(melted_df$variable)

ggplot() +
  geom_histogram(mapping = aes(x = value), data = melted_df, bins = 30) +
  labs(title = "stability of clusters", subtitle = "for 1000 iterations",
       x = "Jaccard index") +
  facet_wrap(vars(variable), ncol=2)

## ----warning=FALSE------------------------------------------------------------
df_gp3 <- bootvalues$group_centers[["group3"]]

melted_df <- reshape2::melt(df_gp3)
melted_df$variable <- as.factor(melted_df$variable)

ggplot() +
  geom_histogram(mapping = aes(x = value), data = melted_df, bins = 30) +
  labs(title = "stability of group 3 centers", subtitle = "for 1000 iterations") +
  xlim(-3,3)+
  facet_wrap(vars(variable), ncol=3)

## ----warning=FALSE, eval = FALSE----------------------------------------------
#  # create the modified weight matrix
#  WMat_adj <- adjustSpatialWeights(Data, WMat$neighbours, style = "W")
#  
#  # calculate the modified version of FCM with non-local information
#  nl_SFCM <- SFCMeans(Data, WMat_adj, k = 4, m = 1.5, alpha = 0.7,
#                   tol = 0.0001, standardize = FALSE,
#                   verbose = FALSE, seed = 456)
#  

