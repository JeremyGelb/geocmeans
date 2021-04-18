---
title: 'geocmeans: A R package for spatial fuzzy c-means'
tags:
  - R
  - unsupervised classification
  - cmeans
  - spatial analysis
  - fuzzy classification
authors:
  - name: Gelb Jeremy
    orcid: 0000-0002-7114-2714
    affiliation: 1
affiliations:
 - name: Laboratoire d'Équité Environnemental, Institut National de la Recherche Scientifique (INRS)
   index: 1

date: 20 Avril 2021
bibliography: references.bib
---

# Summary

Unsupervised classification methods like *k-means* or the *Hierarchical Ascendant Classification* (HAC) are widely used in geography even though they are not well suited for spatial data. Yet, recent development have been proposed to include the geographical dimension into clustering. As an example, `ClustGeo` [@chavent2018clustgeo] is a spatial extension of the HAC, available in the R package with the same name. We propose in the R package `geocmeans` a spatial extension of the *Fuzzy C-Means* (FCM) algorithm to complete this toolbox with a fuzzy approach. The package provides also several helper functions to assess and compare quality of classifications, select appropriate hyper parameters, and interpret the final groups.

# Statement of need

Traditional clustering algorithms do not account for the geographical dimension of data causing two main concerns:

* First, an important part of the information related to observations' locations and geographical organization is discarded. Often, in social, environmental and economical sciences the geographical dimension is highly structuring, in particular in the context of positive spatial autocorrelation.
* Second, in many applications (like the development of regional policy), it is desirable that close observations have more chances to belong to the same group. In this regard, it is common to observe with traditional unsupervised classification algorithms "holes": observations attributed to a certain group and surrounded by observations attributed to another. Often, the difference of the semantic attributes between those observations does not justify such spatial outliers.

A first approach for spatial clustering explored from the 1990s was to impose to the classification a spatial constraint based on contiguity (also called aggregation based approach). The most probably known methods are the `AZP` [@openshaw1977geographical], `SKATER`[@assunccao2006efficient] and `AMOEBA` [@aldstadt2006using]. However, this approach can only yield spatially continuous groups and can be too strict to obtain meaningful clusters. That is why a new approach has been proposed: including the spatial dimension. The spirit is: space should not act as a constraint in the classification but like supplementary data.

The recent method `ClustGeo` [@chavent2018clustgeo] has received some attention from geographers. However, it yield a hard clustering and may hide situations were observations are located at the border of two groups. There is thus a need for available unsupervised spatial fuzzy clustering methods.

Such methods have been discussed and applied in the brain imagery segmentation field [@cai2007fast; @zhao2013kernel] leading to the development of the Spatial Fuzzy C-means (SFCM) algorithm. The package `geocmeans` is proposed to make this tool available to researchers and professionals working with spatial datasets. A first application to construct a socio-residential and environmental taxonomy in Lyon has outlined the potential of the method [@gelb2021apport].

# Core functionality

## Main algorithm

`geocmeans` provides four fuzzy unsupervised classification algorithms:

* `CMeans`, the original c-means algorithm, requiring two hyper parameters, *m* (fuzzyness degree) and *k* (number of groups).

* `GFCMeans`, the so-called generalized c-means algorithm. It is known to accelerate convergence and yield less fuzzy results by adjusting the membership matrix at each iteration. It requires an extra *β* parameter controlling the strength of the modification. The modification only affects the formula updating the membership matrix.

$$u_{ik} = \frac{(||x_{k} - v{_i}||^{2} - \beta_k) ^{(-1/(m-1))}}{\sum_{j=1}^c(||x_{k} - v{_j}||^2 - \beta_k)^{(-1/(m-1))}}$$

with: \
$\beta_k = min(||x_{k} - v||^2)$ \
$0 \leq \beta \leq 1$ \
$u_{ik}$ the probability of observation *k* to belong to cluster *i* \
$x_k$ the observation *k* in the dataset *x* \
$v_i$ the cluster *i* \
*m* the fuzzyness parameter

* `SFCMeans`, the SFCM algorithm, requiring two more parameters *W* and *α*. *W* is a spatial weight matrix used to calculate a spatially lagged version of the dataset *x*. *α* is used to control the weight of the spatially lagged dataset. If $\alpha = 0$ then SFCM degenerates to a simple FCM. If $\alpha = 1$ the same weight is given to the original and lagged dataset. If $\alpha = 2$ then the spatially lagged dataset has a weight doubled in comparison with the original dataset, and so on... The integration of the spatially lagged dataset modifies the formula updating the membership matrix and the formula updating the centers of clusters.

$$u_{ik} = \frac{(||x_{k} - v{_i}||^2 + \alpha||\bar{x_{k}} - v{_i}||^2)^{(-1/(m-1))}}{\sum_{j=1}^c(||x_{k} - v{_j}||^2 + \alpha||\bar{x_{k}} - v{_j}||^2)^{(-1/(m-1))}}$$

$$v_{i} = \frac{\sum_{k=1}^N u_{ik}^m(x_{k} + \alpha\bar{x_{k}})}{(1 + \alpha)\sum_{k=1}^N u_{ik}^m}$$

with:

$0 \leq \alpha \leq \infty$ \
$v_i$ the cluster *i* \
$\bar{x}$ the spatially lagged version of *x* \

As the formula suggests, the SFCM can be seen as a spatially smoothed version of the FCM and *α* controls the degree of spatial smoothness. This smoothing can be interpreted as an attempt to reduce spatial overfitting of the FCM. 

* `SGFCMeans`, the SGFCM algorithm, combining SFCM and SGFCM and thus requiring the definition of three extra parameters *W*, *α* and *β*. Only the formula to calculate the membership matrix is different from the SFCM.

$$u_{ik} = \frac{(||x_{k} - v{_i}||^2 -\beta_k + \alpha||\bar{x_{k}} - v{_i}||^2)^{(-1/(m-1))}}{\sum_{j=1}^c(||x_{k} - v{_j}||^2 -\beta_k + \alpha||\bar{x_{k}} - v{_j}||^2)^{(-1/(m-1))}}$$

## Selecting parameters

As stated above, up to five hyper parameters have to be selected by the user. Finding the best combination is facilitated by the function `selectParameters` calculating the classifications for all the possible combinations of parameters in specified ranges and returning several metrics of classification quality.

## Interpreting the results

To interpret the results, four functions are provided: 

* `summarizeClusters`: returning summary statistics for each cluster for in depth analysis.
* `spiderPlots`: plotting a spider chart to quickly differentiate  the clusters.
* `violinPlots`: plotting one violin plot split by cluster for each variable in the dataset.
* `mapClusters`: mapping the membership matrix and the most likely cluster for the observations .
* `calcqualityIndexes`: returning several quality indices for a classification.
* `spatialDiag`: performing a complete spatial diagnostic to determine if the inclusion of space in the classification is justified.

# Example

The data used for the socio-residential and environmental taxonomy in Lyon is included in the package. The following example uses this data to demonstrate the basic functionality of the package. More details are given in the vignettes of the package.

```r
library(geocmeans)
library(spdep)

data(LyonIris)

#selecting the columns for the analysis
AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14",
                   "Pct_65","Pct_Img","TxChom1564","Pct_brevet","NivVieMed")
                   
#creating a spatial weights matrix
Neighbours <- poly2nb(LyonIris,queen = TRUE)
WMat <- nb2listw(Neighbours,style="W",zero.policy = TRUE)

#rescaling the columns
Data <- LyonIris@data[AnalysisFields]
for (Col in names(Data)){
  Data[[Col]] <- scale(Data[[Col]])
}

#considering k = 4 and m = 1.5, find an optimal value for alpha
DFindices_SFCM <- selectParameters(algo = "SFCM", data = Data,
                               k = 4, m = 1.5, alpha = seq(0,2,0.05),
                               nblistw = WMat, standardize = FALSE,
                               tol = 0.0001, verbose = FALSE, seed = 456)

#keeping alpha = 0.7
SFCM_results <- SFCMeans(Data, WMat, k = 4, m = 1.5, alpha = 0.7,
                 tol = 0.0001, standardize = FALSE,
                 verbose = FALSE, seed = 456)

#calculating some quality indexes
calcqualityIndexes(Data, SFCM_results$Belongings)

#mapping the results
mapClusters(LyonIris, SFCM_results$Belongings)

```

# Acknowledgements

We are grateful to Professor Philippe Apparicio for his comments and suggestions on the package, its documentation and this article.

# References
