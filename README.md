
<!-- README.md is generated from README.Rmd. Please edit that file -->

# geocmeans <img src='man/figures/geocmeans_logo.png' align="right" height="138.5" />

## An R package to perform Spatial Fuzzy C-means.

<!-- badges: start -->

    [![R-CMD-check](https://github.com/JeremyGelb/geocmeans/workflows/R-CMD-check/badge.svg)](https://github.com/JeremyGelb/geocmeans/actions)

<!-- badges: end -->

## Installation

You can install a development version of the `geocmeans` package using
the command below.

    remotes::install_github(repo = "JeremyGelb/geocmeans", build_vignettes = TRUE, force = TRUE)

## Authors

Jeremy Gelb, Laboratoire d’Équité Environnemental INRS (CANADA), Email:
<jeremy.gelb@ucs.inrs.ca>

Philippe Apparcio, Laboratoire d’Équité Environnemental INRS (CANADA),
Email: <philippe.apparicio@ucs.inrs.ca>

## About the package

Provides functions to apply Spatial Fuzzy c-means Algorithm, visualize
and interpret results. This method is well suited when the user wants to
analyze data with a fuzzy clustering algorithm and to account for the
spatial dimension of the dataset. Tests for estimating the spatial
consistency and classification quality are proposed in addition. The
algorithms were developped first for brain imagery as described in the
articles of [Cai and
al. 2007](https://doi.org/10.1016/j.patcog.2006.07.011) and [Zaho and
al. 2013](https://doi.org/10.1016/j.dsp.2012.09.016). [Gelb and
Apparicio](https://doi.org/10.4000/cybergeo.36414) proposed to apply the
method to perform a socio-residential and environmental taxonomy in Lyon
(France).

Approaches for visualising uncertainty in spatial data are presented in
this package. These include the three approaches developed in [Lucchesi
and Wikle
(2017)](http://faculty.missouri.edu/~wiklec/LucchesiWikle2017Stat) and a
fourth approach presented in [Kuhnert et
al. (2018)](https://publications.csiro.au/publications/#publication/PIcsiro:EP168206).

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
[guidelines](CONTRIBUTING.md).

Please note that the `Vizumap` project is released with a [Contributor
Code of Conduct](CONDUCT.md). By contributing to this project, you agree
to abide by its terms.

## License

`Vizumap` version 1.2.0 is licensed under [GPLv3](LICENSE.md).

## History of Vizumap

Vizumap began as a visualisation project at the University of Missouri
in 2016. Chris Wikle, professor of statistics, posed an interesting
research question to Lydia Lucchesi, a student curious about data
visualisation and R.

How do you include uncertainty on a map displaying areal data estimates?

Over the course of a year, they put together three methods for
visualising uncertainty in spatial statistics: the bivariate choropleth
map, the pixel map, and the glyph map. By mid-2017, there were maps, and
there was a lot of R code, but there was not a tool that others could
use to easily make these types of maps, too. That’s when statistician
Petra Kuhnert recommended developing an R package. Over the course of a
month, Petra and Lydia developed Vizumap (originally named VizU) at
CSIRO Data61 in Canberra, Australia. Since then, the package has been
expanded to include exceedance probability maps, an uncertainty
visualisation method developed by Petra while working on a Great Barrier
Reef (GBR) project.

Vizumap has been used to visualise the uncertainty of American Community
Survey estimates, the prediction errors of sediment estimates in a GBR
catchment, and most recently the [uncertainty of estimated locust
densities in
Australia](https://www.nature.com/articles/s41598-020-73897-1/figures/4).
We would like to assemble a Vizumap gallery that showcases different
applications of the package’s mapping methods. If you use Vizumap to
visualise uncertainty, please feel free to send the map our way. We
would like to see it!

## References

Kuhnert, P.M., Pagendam, D.E., Bartley, R., Gladish, D.W., Lewis, S.E.
and Bainbridge, Z.T. (2018) [Making management decisions in face of
uncertainty: a case study using the Burdekin catchment in the Great
Barrier Reef, Marine and Freshwater
Research](https://publications.csiro.au/publications/#publication/PIcsiro:EP168206),
69, 1187-1200, <https://doi.org/10.1071/MF17237>.

Lucchesi, L.R. and Wikle C.K. (2017) [Visualizing uncertainty in areal
data with bivariate choropleth maps, map pixelation and glyph
rotation](http://faculty.missouri.edu/~wiklec/LucchesiWikle2017Stat),
Stat, <https://doi.org/10.1002/sta4.150>.
