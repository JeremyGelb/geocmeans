#' @param maxiter An integer for the maximum number of iterations
#' @param tol The tolerance criterion used in the evaluateMatrices function for
#'   convergence assessment
#' @param standardize A boolean to specify if the variables must be centred and
#'   reduced (default = True)
#' @param verbose A boolean to specify if the progress should be printed
#' @param init A string indicating how the initial centres must be selected. "random"
#' indicates that random observations are used as centres. "kpp" use a distance-based method
#' resulting in more dispersed centres at the beginning. Both of them are heuristic.
#' @param seed An integer used for random number generation. It ensures that the
#' starting centres will be the same if the same value is selected.
