# geocmeans

<img src="man/figures/geocmeans_logo.png" width = 120 alt="geocmeans Logo"/>

A R package to perform spatially constrained cmeans.

## Getting Started

A good start point for this package is the Introduction vignette. It presents the main features


### Installing

you can install this package with the following code in R.
The packages use mainly the following packages in its internal structure :

* sp
* spdep
* fclust
* future
* future.apply
* dplyr
* ggplot2

```{r}
devtools::install_github("JeremyGelb/geocmeans")
```

### Examples

We provide here some short examples of the main features

* calculating a Spatial Fuzzy C-Means (SFCM)

```{r}
library(geocmeans)
data(LyonIris)
AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img","TxChom1564","Pct_brevet","NivVieMed")
dataset <- LyonIris@data[AnalysisFields]
queen <- spdep::poly2nb(LyonIris,queen=TRUE)
Wqueen <- spdep::nb2listw(queen,style="W")
result <- SFCMeans(dataset, Wqueen,k = 5, m = 1.5, alpha = 1.5, standardize = TRUE)
```

* select the best parameters for the classification
```{r}
future::plan(future::multiprocess(workers=2))
values <- select_parameters.mc(dataset, k = 5, m = seq(2,3,0.1),
alpha = seq(0,2,0.1), nblistw = Wqueen)

```

* calcuating some indices to asses the quality of the classification

```{r}
qual_indices <- calcqualityIndexes(result$Data, result$Belongings, m=1.5)
sp_indices <- spatialDiag(result$Belongings, Wqueen, nrep=150)
```

* build some maps and charts to inspect results

```{r}
all_maps <- mapClusters(LyonIris, result$Belongings)
spiderPlots(dataset,result$Belongings)
violinPlots(dataset, result$Groups)
```

* calculate summary statistics for each group

```{r}
summarizeClusters(dataset, result$Belongings)
```

## Authors

* **Jeremy Gelb** - *Creator and maintener*


## License

This project is licensed under the GPL2 License

## Acknowledgments

* Hat tip to Hadley Wickham and its helpfull book *R packages*

