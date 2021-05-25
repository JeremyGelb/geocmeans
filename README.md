
<!-- README.md is generated from README.Rmd. Please edit that file -->

# geocmeans <img src='man/figures/geocmeans_logo.png' align="right" height="138.5" />

## An R package to perform Spatial Fuzzy C-means.

<!-- badges: start -->

[![R-CMD-check](https://github.com/JeremyGelb/geocmeans/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JeremyGelb/geocmeans/actions/workflows/R-CMD-check.yaml)
[![](https://img.shields.io/badge/devel%20version-0.1.1-green.svg)](https://github.com/JeremyGelb/geocmeans)
[![](https://www.r-pkg.org/badges/version/geocmeans?color=blue)](https://cran.r-project.org/package=geocmeans)
[![](http://cranlogs.r-pkg.org/badges/grand-total/geocmeans?color=blue)](https://cran.r-project.org/package=geocmeans)
<!-- badges: end -->

## Installation

The stable version of `geocmeans` is available on CRAN. You can install
it with the command below.

    install.packages("geocmeans")

You can install a development version of the `geocmeans` package using
the command below.

    remotes::install_github(repo = "JeremyGelb/geocmeans", build_vignettes = TRUE, force = TRUE)

## Authors

Jeremy Gelb, Laboratoire d’Équité Environnemental INRS (CANADA), Email:
<jeremy.gelb@ucs.inrs.ca>

## Contributors

Philippe Apparcio, Laboratoire d’Équité Environnemental INRS (CANADA),
Email: <philippe.apparicio@ucs.inrs.ca>

## About the package

Provides functions to apply Spatial Fuzzy c-means Algorithm, visualize
and interpret results. This method is well suited when the user wants to
analyze data with a fuzzy clustering algorithm and to account for the
spatial dimension of the dataset. Indexes for measuring the spatial
consistency and classification quality are proposed in addition. The
algorithms were developed first for brain imagery as described in the
articles of [Cai and
al. 2007](https://doi.org/10.1016/j.patcog.2006.07.011) and [Zaho and
al. 2013](https://doi.org/10.1016/j.dsp.2012.09.016). [Gelb and
Apparicio](https://doi.org/10.4000/cybergeo.36414) proposed to apply the
method to perform a socio-residential and environmental taxonomy in Lyon
(France).

#### Fuzzy classification algorithms

Four Fuzzy classification algorithms are proposed :

-   FCM: Fuzzy C-Means, with the function `CMeans`
-   GFCM: Generalized Fuzzy C-Means, with the function `GFCMeans`
-   SFCM: Spatial Fuzzy C-Means, with the function `SFCMeans`
-   SGFCM: Spatial Generalized Fuzzy C-Means, with the function
    `SGFCMeans`

Each function return a membership matrix, the data used for the
classification (scaled if required) and the centers of the clusters.

#### Parameter selections

The algorithms available require different parameters to be fixed by the
user. The function `selectParameters` is a useful tool to compare the
results of different combinations of parameters. A multicore version,
`selectParameters.mc`, using a plan from the package `future` is also
available to speed up the calculus.

#### Classification quality

Many indices of classification quality can be calculated with the
function `calcqualityIndexes`:

-   *Silhouette.index*: the silhouette index (fclust::SIL.F)
-   *Partition.entropy*: the partition entropy index (fclust::PE)
-   *Partition.coeff*: the partition entropy coefficient (fclust::PC)
-   *Modified.partition.coeff*: the modified partition entropy
    coefficient (fclust::MPC)
-   *XieBeni.index* : the Xie and Beni index (fclust::XB)
-   *FukuyamaSugeno.index* : the Fukuyama and Sugeno index
    (geocmeans::calcFukuyamaSugeno)
-   *Explained.inertia*: the percentage of total inertia explained by
    the solution

#### Interpretation

Several functions are also available to facilitate the interpration of
the classification:

-   summary statistics for each cluster: `summarizeClusters`
-   spider charts: `spiderPlots`
-   violin plots: `violinPlots`
-   maps of the membership matrix: `mapClusters` (support polygon,
    points and polylines)

#### Spatial inconsistency

We proposed an index to quantify the spatial inconsistency of a
classification ([Gelb and
Apparicio](https://doi.org/10.4000/cybergeo.36414)). If in a
classification close observations tend to belong to the same group, then
the value of the index is close to 0. If the index is close to 1, then
the belonging to groups is randomly distributed in space. A value higher
than one can happen in the case of negative spatial autocorrelation. The
index is described in the vignette `adjustinconsistency`. The function
`spatialDiag` does a complete spatial diagnostic of the membership
matrix resulting from a classification.

## Examples

Detailed examples are given in the vignette `introduction`

    vignette("introduction","geocmeans")

## Testing

If you would like to install and run the unit tests interactively,
include `INSTALL_opts = "--install-tests"` in the installation code.

    remotes::install_github(repo = "JeremyGelb/geocmeans", build_vignettes = TRUE, force = TRUE, INSTALL_opts = "--install-tests")
    testthat::test_package("geocmeans", reporter = "stop")

## Contribute

To contribute to `geocmeans`, please follow these
[guidelines](https://github.com/JeremyGelb/geocmeans/blob/master/CONTRIBUTING.md).

Please note that the `geocmeans` project is released with a [Contributor
Code of
Conduct](https://github.com/JeremyGelb/geocmeans/blob/master/CONDUCT.md).
By contributing to this project, you agree to abide by its terms.

## License

`geocmeans` version 0.1.0 is licensed under [GPL2
License](https://github.com/JeremyGelb/geocmeans/blob/master/LICENSE.txt).
