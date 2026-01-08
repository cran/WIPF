#' Weighted Iterative Proportional Fitting (WIPF) in one dimension
#'
#' @description  Implements WIPF in one dimension. This function updates, using a set of weights,
#'               an initial 1-dimensional array, a vector (referred as the seed), to match a given value
#'               (referred as the margin), in such as way that the weighted sum of the updated values coincide with the margin.
#'
#' @author Jose M. Pavia, \email{pavia@@uv.es}
#'
#' @param seed A vector of non-negative values with the initial values.
#'
#' @param weights A vector of non-negative values with the weights associated to each component of `seed` and with the same length as `seed`.
#'
#' @param margin A non-negative scalar with the (weighted) marginal total to be fitted. Default, `1`.
#'
#' @param normalize `TRUE`/`FALSE` argument indicating whether weights should be normalized
#'                   to sum 1 before building the weighted sum to be compared with the margin value.
#'                   Default, `TRUE`. Normalization is necessary when we want to adjust a set of indexes
#'                   for which the margin correspond to other index that is a theoretical convex
#'                   combination of the inner indexes. This characterizes a typical context where WIPF could be of value.
#'
#' @param tol Stopping criterion. The algorithm stops when the maximum absolute difference between
#'            the solutions obtained in two consecutive iteration is lower than the value specified by `tol`.
#'            Default, `0.000001`.
#'
#' @param maxit Stopping criterion. A positive integer number indicating the maximum number of iterations
#'              allowed. Default, `1000`. The algorithm will stop if the values to be fitted still has not
#'              converged after this many iterations.
#'
#' @param full `TRUE`/`FALSE` argument indicating if either only the solution should be shown or
#'                            a more complete output.
#'
#' @param ... Other arguments to be passed to the function. Not currently used.
#'
#' @return
#' When `full = FALSE` an object similar to `seed` with the solution reached when the algorithm stops.
#' When `full = TRUE` a list with the following components:
#'  \item{sol}{ An object similar to `seed` with the solution reached at convergence (or when the maximum number of iterations is reached).}
#'  \item{iter}{ Number of iterations when the algorithm stops.}
#'  \item{error.margins}{ An object similar to `margin` with the absolute differences between the values in `margin` and
#'                        the weighted sum of the values in `sol`.}
#'  \item{inputs}{ A list containing all the objects with the values used as arguments by the function.}
#'
#' @note Weighted Iterative proportional fitting is an extension of IPF.
#'       WIPF produces the same solutions than IPF with all weights being ones
#'       and when they are not normalized. IPF is also known as RAS in economics,
#'       raking in survey research or matrix scaling in computer science.
#'
#' @export
#'
#' @examples
#' s <- c(1.0595723, 0.9754876, 0.8589494, 0.8589123)
#' w <- c(651301.9, 581185.1, 555610.8, 602595.6)
#' example <- WIPF1(seed = s, weights = w)
#'

WIPF1 <- function(seed, weights,
                  margin = 1,
                  normalize = TRUE,
                  tol = 10^-6,
                  maxit = 1000,
                  full = FALSE,
                  ...)
  {
  # checking if a non-negative values in seed
  if (min(seed) < 0) {
    stop('Error: the object "seed" contains negative values.')
  }
  # checking a non-negative values
  if (min(weights) < 0 | max(weights) == 0) {
    stop('Error: non-negative values with at least a positive value is required in "weights".')
  }
  if(!is.vector(seed) | length(seed) < 2){
    stop('Error: the object "seed" must be a vector of at least two components.')
  }
  if(!is.vector(weights) | length(weights) < 2){
    stop('Error: the object "weights" must be a vector of at least two components.')
  }
  if(length(seed) !=  length(weights)){
    stop('Error: the objects "seed" and weights" must have the same length.')
  }
  if(length(margin) != 1L | min(margin) < 0){
    stop('Error: the object "margin" must be a non-negative constant.')
  }

  inputs <- c(as.list(environment()), list(...))
  if(normalize){
    w <- weights/sum(weights)
  } else {
    w <- weights
  }
  delta <- delta0 <- seed
  suma <- sum(delta*w)
  dif <- abs(margin - suma)
  iter <- 0L
  while(dif > tol & iter < maxit){
    delta <- delta*(margin/suma)
    suma <- sum(delta*w)
    dif <- max(abs(delta - delta0))
    delta0 <- delta
    iter <- iter + 1L
  }
  if(full){
    return(list("sol" = delta, "iter" = iter, "error.margins" = margin - suma, "inputs" = inputs))
  } else {
    return(delta)
  }
}
