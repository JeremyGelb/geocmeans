#' @details
#' The C-means algorithm (FCM) is a classical method, close to the k-means algorithm
#' but it yields fuzzy results. For each observation, the probability that it
#' belongs to each cluster is calculated and stored in a membership matrix.
#'
#' The generalized version of the FCM (GFCM) can be used to obtain less fuzzy results
#' and to speed up convergence. A parameter (beta) is introduced in the calculation of
#' the membership matrix to reduce the membership value to the furthest cluster of each
#' observation. beta varies between 0 (classical FCM) and 1.
#'
#' The robust version of the FCM is obtained by normalizing for each group the
#' euclidean distances between the observations and their centre at each iteration.
#' It improves the performance of the FCM when clusters do not have the same density
#' or do not have hyperspherical shapes.
