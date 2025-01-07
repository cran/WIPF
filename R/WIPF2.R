#' Weighted Iterative Proportional Fitting (WIPF) in two dimensions
#'
#' @description  Implements WIPF in two dimensions. This function updates an initial
#'               2-dimensional array (a matrix, referred to as the seed) using a matrix
#'               of weights to align with a set of two vectors (referred to as the margins),
#'               where one of them can be missing. When `margin1` and `margin2` are compatible
#'               given the weights, the updated values ensure that the weighted sum
#'               across columns matches `margin1` and the weighted sum across rows
#'               matches `margin2`.
#'               If the margins are incompatible given the weights, the function `WIPF1`
#'               is applied to the initial margins to make the margins compatible
#'               with the weights. In those cases, margins are updated (are made compatible)
#'               in increasing order of sub-indices (i.e., `margin2` is adjusted
#'               to be compatible with `margin1`).
#'
#' @author Jose M. Pavia, \email{pavia@@uv.es}
#'
#' @param seed A (RxC) matrix of non-negative values with the initial values.
#'
#' @param weights A (RxC) matrix of non-negative values with the weights associated to each entry of the matrix.
#'
#' @param margin1 A R-length vector of positive values with the target (weighted) marginal sums across columns to be fitted.
#'
#' @param margin2 A C-length vector of positive values with the target (weighted) marginal sums across rows to be fitted.
#'
#' @param normalize Logical (`TRUE`/`FALSE`) argument indicating whether the weights should be normalized
#'                  (across all dimensions, for either row or column weights to sum 1)  before constructing
#'                  the weighted sums for comparison with the margin values. Default, `TRUE`. Normalization is essential when adjusting
#'                  a set of indexes where the margins represent theoretical convex combinations of the inner indexes.
#'                  This characterizes a typical context where WIPF could be of value.
#'
#' @param tol Stopping criterion. The algorithm stops when the maximum difference between the weighted sum(s) of
#'            the values to be fitted and the margin(s) is lower than the value specified by `tol`.
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
#'  \item{dev.margins}{ A list with a set of objects similar to the margins with absolute maximum deviations
#'                        between the values in margins and the corresponding weighted sums of the values in `sol`.}
#'  \item{margin1}{ A R-length vector of positive values with the actual margin1 object used to reach the solution.
#'                    This coincides with `margin1` when all the margins are compatible given the weights.}
#'  \item{margin2}{ A C-length vector of positive values with the actual margin2 object used to reach the solution.
#'                    This coincides with `margin2` when all the margins are compatible given the weights.}
#'  \item{inputs}{ A list containing all the objects with the values used as arguments by the function.}
#'
#' @note Weighted Iterative proportional fitting is an extension of IPF.
#'       WIPF produces the same solutions than IPF with all weights being ones
#'       and when they are not normalized. IPF is also known as RAS in economics,
#'       raking in survey research or matrix scaling in computer science.

#' @export
#'
#' @examples
#' s <- structure(c(1.1279, 1.1304, 1.0304, 0.8554, 1.5606, 1.4171, 1.2862,
#'                  1.2472, 1.0746, 1.0796, 0.9806, 0.928, 1.1607, 1.2436, 1.2191,
#'                  1.0786, 1.0194, 1.1716, 0.9937, 0.8611, 1.0172, 1.2511, 1.1606,
#'                  1.1959), .Dim = c(4L, 6L))
#' w <- structure(c(72161.97, 93725.94, 84408.83, 172774.13, 52875.08,
#'                  31936.92, 14191.44, 12595.46, 291698.94, 231408.32,
#'                  221763.43, 235217.74, 42028.56, 64458.09, 93443.13,
#'                  60348.74, 222482.04, 103695.94, 57066.82, 48657.48,
#'                  9572.75, 75745.02, 83912.38, 94019.92), .Dim = c(4L, 6L))
#' m1 <- c(1.110737, 1.029947, 0.934799, 0.906475)
#' m2 <- c(0.810992, 1.375921, 1.071519, 1.045006, 0.949938, 0.915762)
#' example <- WIPF2(seed = s, weights = w, margin1 = m1, margin2 = m2, full = TRUE)
#'

WIPF2 <- function(seed,
                  weights,
                  margin1,
                  margin2,
                  normalize = TRUE,
                  tol = 10^-6,
                  maxit = 1000,
                  full = FALSE,
                  ...)
  {

  # checking non-negative values
  if (min(seed) < 0) {
    stop('Error: the object "seed" contains negative values.')
  }
  if (min(weights) < 0 | min(apply(weights, 1, max)) == 0 | min(apply(weights, 2, max)) == 0) {
    stop('Error: non-negative values with at least a positive value in each row and column is required in "weights".')
  }
  # checking if margin1 is provided
  if (missing(margin1)) {
    margin1 <- NULL
  }
  # checking if margin2 is provided
  if (missing(margin2)) {
    margin2 <- NULL
  }

  if (!is.null(margin1)) {
    if (min(margin1) <= 0)
      stop('Error: the object "margin1" contains non-positive values.')
  }
  if (!is.null(margin2)) {
    if (min(margin2) <= 0)
      stop('Error: the object "margin2" contains non-positive values.')
  }

  # checking dimensions
  if(length(dim(seed)) != 2){
    stop('Error: the object "seed" must be a matrix.')
  }
  if(length(dim(weights)) != 2){
    stop('Error: the object "weights" must be a matrix.')
  }
  if(!all(dim(seed) == dim(weights))){
    stop('Error: the objects "seed" and weights" must have the same dimensions.')
  }
  if(length(margin1) !=0){
    if(nrow(seed) != length(margin1))
      stop('Error: the dimensions of "seed" and margin1" are incompatible.')
  }
  if(length(margin2) !=0){
    if(ncol(seed) != length(margin2))
      stop('Error: the dimensions of "seed" and margin2" are incompatible.')
  }

  inputs <- c(as.list(environment()), list(...))

  delta <- as.matrix(seed)
  weights <- as.matrix(weights)

  if(normalize){
    W <- weights/rowSums(weights)
    V <- t(t(weights)/rowSums(t(weights)))
    w <- weights/sum(weights)
    pesos1 <- rowSums(w)
    pesos2 <- colSums(w)
  } else {
    V <- W <- weights
    pesos1 <- pesos2 <- 1
  }

  if (is.null(margin1) & is.null(margin2)){
    iter <- 0
    dif <- NULL
  } else if (is.null(margin1) & !is.null(margin2)){
    iter <- 1
    vector.columna <- colSums(delta * V)
    S1 <- diag(as.vector(margin2/vector.columna))
    S1[is.nan(S1)] <- 0L
    S1[is.infinite(S1)] <- 1L
    delta <- delta %*% S1
    dif <- list("dev2" = rep(0, length(margin2)))
  } else if (!is.null(margin1) & is.null(margin2)){
    iter <- 1
    vector.fila <- rowSums(delta * W)
    R1 <- diag(as.vector(margin1/vector.fila))
    R1[is.nan(R1)] <- 0L
    R1[is.infinite(R1)] <- 1L
    delta <- R1 %*% delta
    dif <- list("dev1" = rep(0, length(margin1)))
  } else {
    if (abs(sum(margin1*pesos1) - sum(margin2*pesos2)) > tol){
      warning("The values for arguments 'margin1' and 'margin2' are not compatible given the 'weights', 'margin2' have been updated with WIPF1 using 'margin1' as reference.")
      # margin1 <- WIPF1(seed = margin1, weights = rowSums(W),
      #                    normalize = normalize, tol = tol,
      #                   maxit = maxit)
      margin2 <- WIPF1(seed = margin2, weights = colSums(W),
                       margin = sum(margin1*pesos1), normalize = normalize,
                       tol = tol, maxit = maxit)
    }

    vector.fila <- rowSums(delta*W)
    vector.columna <- colSums(delta*V)

    dif <- max(c(abs(vector.fila - margin1), abs(vector.columna - margin2)))
    iter <- 0L

    while(dif > tol & iter < maxit){

      R1 <- diag(as.vector(margin1/vector.fila))
      R1[is.nan(R1)] <- 0L
      R1[is.infinite(R1)] <- 1L
      delta <- R1 %*% delta

      vector.columna <- colSums(delta*V)
      S1 <- diag(as.vector(margin2/vector.columna))
      S1[is.nan(S1)] <- 0L
      S1[is.infinite(S1)] <- 1L
      delta <- delta %*% S1

      vector.fila <- rowSums(delta*W)
      vector.columna <- colSums(delta*V)

      dif <- max(c(abs(vector.fila - margin1), abs(vector.columna - margin2)))
      iter <- iter + 1L
    }
    delta[is.nan(delta)] <- 0L
    delta[is.infinite(delta)] <- 1L
    dif <- list("dev1" = abs(vector.fila - margin1),
                "dev2" = abs(vector.columna - margin2))
  }

  if(full){
    return(list("sol" = delta, "iter" = iter, "dev.margins" = dif,
                "margin1" = margin1, "margin2" = margin2, "inputs" = inputs))
  } else {
    return(delta)
  }

}
