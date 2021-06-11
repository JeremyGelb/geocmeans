# geocmeans 0.1.1.9000

## New Features

* Added a function to predict the membership matrix of a new set of observations (*predict.gcm*)
* Added a shiny app (function : *sp_clust_explorer*) for result exploration
* The results of the functions *CMeans*, *GCMeans*, *SFCMeans*, *SGFCMeans* are now objects of class *FCMres* and the generic functions *predict* and *summary* can be used on them.
* Added some cluster quality indices : Negentropy Increment index, Generalized Dunn’s index (43 and 53), David-Bouldin index, Calinski-Harabasz index
* Added a function to perform clustering validation by boostrap (see function *bstp_group_validation*)

## corrected bugs

* issue 1 fixed by editing the mapping functions. A bug arrised when the fid of a SpatialDataFrame read from a shapefile was different from 1:nrow(df)


# geocmeans 0.1.1

* Accepted CRAN version
