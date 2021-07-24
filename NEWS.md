# geocmeans 0.1.1.9000

## New Features

* Added support to use raster data for clustering (see vignette *rasters*)
* Added a function to predict the membership matrix of a new set of observations (*predict.FCMres*)
* Added a shiny app (function : *sp_clust_explorer*) for result exploration
* The results of the functions *CMeans*, *GCMeans*, *SFCMeans*, *SGFCMeans* are now objects of class *FCMres* and the generic methods *predict*, *summary*, *is* and *print* can be used on them. *FMCres* object can easily be created by hand from results of other classifier if needed, see the new vignette *FMCres*.
* Added some cluster quality indices : Negentropy Increment index, Generalized Dunn’s index (43 and 53), David-Bouldin index, Calinski-Harabasz index
* Added a function to perform clustering validation by bootstrap (see function *bstp_group_validation*)
* Added a function to reorder the results of a classification to match the most similar groups in a second classification (*groups_matching*)
* Added functions to evaluate spatial autocorrelation of a classification results: ELSA and FuzzyELSA (see functions *calcELSA* and *calcFuzzyELSA*)

## corrected bugs

* issue 1 fixed by editing the mapping functions. A bug occurred when the fid of a SpatialDataFrame read from a shapefile was different from 1:nrow(df)

## performance

* an important performance gain can be observed for large dataset, the function to compare matrices between two iterations is now significantly faster.
* Core functions rewritten with Rcpp for massive time gain

# geocmeans 0.1.1

* Accepted CRAN version
