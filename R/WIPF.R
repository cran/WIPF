#' Weighted Iterative Proportional Fitting (WIPF) in N (greater than 1) dimensions
#'
#' @description Implements the Weighted Iterative Proportional Fitting (WIPF)
#'              algorithm to adjust an N-dimensional array subject to weighted
#'              marginal constraints.
#'
#' @details The function updates an initial N-dimensional array (referred to as
#'          the seed) using an N-dimensional array of weights to align with a
#'          collection of lower-dimensional margins of different orders, some of
#'          which may be missing.
#'          When all provided margins are compatible given the weights, the
#'          updated array ensures that the corresponding weighted sums over all
#'          specified index subsets coincide with the supplied margins.
#'          If the provided margins are incompatible given the weights, the
#'          functions `WIPF1` and `WIPF` are applied to the provided margins to
#'          guarantee their compatibility with the weights. In these cases, the margins
#'          are updated in increasing order of sub-indices, with higher-order
#'          indices running faster.
#'
#' @author Jose M. Pavia, \email{pavia@@uv.es}
#'
#' @param seed An N-dimensional array of non-negative values with the initial values.
#'
#' @param weights An N-dimensional array of non-negative values with the weights associated to each entry of the `seed` array.
#'
#' @param margins A list of arrays of non-negative values of order less than N with the
#'                target (weighted) marginal sums over the specified dimensions.
#'
#' @param indices A list of vectors with the same length as `margins`, indicating the indices
#'                corresponding to each array in `margins`.
#'
#' @param normalize `TRUE`/`FALSE` argument indicating whether the weights should be normalized
#'                  across all dimensions before constructing the weighted sums for comparison
#'                  with the margin values. Default is `TRUE`.
#'                  Normalization is essential when adjusting a set of indices
#'                  for which the margins represent theoretical convex combinations of the inner
#'                  indices. This characterizes a typical context in which WIPF may be of value.
#'
#' @param tol Stopping criterion. The algorithm stops when the maximum difference between the weighted sums of
#'            the values to be fitted and the margins is lower than the value specified by `tol`.
#'            Default, `0.000001`.
#'
#' @param maxit Stopping criterion. A positive integer number indicating the maximum number of iterations
#'              allowed. Default, `1000`. The algorithm will stop if the values to be fitted still has not
#'              converged after these many iterations.
#'
#' @param full `TRUE`/`FALSE` argument indicating if either only the solution should be saved or
#'              a more complete output. Default `FALSE`.
#'
#' @param ... Other arguments to be passed to the function. Not currently used.
#'
#' @return
#' When `full = FALSE` an object similar to `seed` with the solution reached when the algorithm stops.
#' When `full = TRUE` a list with the following components:
#'  \item{sol}{ An object similar to `seed` with the solution reached at convergence (or when the maximum number of iterations is reached).}
#'  \item{iter}{ Number of iterations when the algorithm stops.}
#'  \item{dev.margins}{ A list of arrays similar to the `margins` output containing the absolute maximum deviations
#'                      between the values in `margins` and the corresponding weighted sums of the values in `sol`.}
#'  \item{margins}{ A list of arrays similar `margins` with the actual margins used to reach the solution.
#'                    Each array whose margins are compatible given the weights coincides with the original array.}
#'  \item{inputs}{ A list containing all the objects with the values used as arguments by the function.}
#'
#' @note Weighted Iterative Proportional Fitting is an extension of IPF.
#'       WIPF produces the same solutions than IPF with all weights being ones
#'       and when they are not normalized. IPF is also known as RAS in economics,
#'       raking in survey research or matrix scaling in computer science.
#'
#' @export
#'
#' @examples
#' s <- structure(c(0.9297, 0.9446, 0.8763, 0.92, 0.8655, 0.8583, 0.8132,
#'                  0.8679, 0.7968, 0.7834, 0.721, 0.7859, 0.7747, 0.7851, 0.8632,
#'                  1.041, 1.5617, 1.5642, 1.4847, 1.5176, 1.4157, 1.3851, 1.3456,
#'                  1.4012, 1.3017, 1.2626, 1.1904, 1.2668, 1.3203, 1.3181, 1.1965,
#'                  1.1654, 1.2219, 1.3863, 1.306, 1.1963, 1.1376, 1.35, 1.2595,
#'                  1.1289, 1.0456, 1.2863, 1.1274, 1.0208, 1.0542, 1.1272, 1.1594,
#'                  1.1668, 1.1931, 1.1328, 1.1221, 1.1011, 1.1298, 1.0454, 1.0573,
#'                  1.0557, 1.0599, 0.973, 0.9545, 0.9721, 1.0489, 0.9934, 0.9382,
#'                  0.876, 1.339, 1.1939, 1.0229, 1.0378, 1.0402, 0.9554, 0.9794,
#'                  1.0089, 0.9422, 0.8584, 0.8563, 0.9013, 0.9252, 0.8706, 0.8354,
#'                  0.8071, 0.9737, 1.0008, 0.9593, 0.9257, 0.9556, 0.9534, 0.9313,
#'                  0.9151, 0.883, 0.8731, 0.8285, 0.8309, 0.9131, 0.9258, 0.8467,
#'                  0.7785), .Dim = c(4L, 4L, 6L))
#' w <- structure(c(18520.3, 11776.3, 19479.5, 22497.6, 18968.7, 17263.7,
#'                  36494.7, 21707, 13406.3, 13570.4, 37746.1, 20593.2, 6595.6, 25444.6,
#'                  59868.2, 81777.2, 3380.4, 20610.7, 22247.3, 6800.9, 5236.3, 14877.8,
#'                  7205, 5028.4, 1130.7, 6603.2, 4007.4, 2620.5, 374.8, 1624.3,
#'                  4963.7, 9551.3, 31806, 93615.9, 121986.6, 44640.3, 32110.6, 95814.4,
#'                  72827.9, 30922.5, 43197.3, 72050.8, 66673.4, 40370.1, 31488.2,
#'                  55014.9, 69457.2, 80021.2, 17701.7, 8765.2, 11790.9, 3872.8,
#'                  30544.5, 12141.2, 12415.2, 9471.9, 36138.6, 19198.1, 23120.1,
#'                  15597.9, 12140.2, 8058.3, 20948.3, 19380.2, 78543.9, 86503.6,
#'                  28727.8, 29208.7, 26300.6, 42363, 20786.6, 14380.3, 9493.5, 17816.2,
#'                  19844.1, 10898.2, 1419, 4211.5, 20615, 22748.2, 3365.8, 2639.8,
#'                  2433.3, 930.5, 22119.6, 31022.7, 12748.5, 10161.4, 15450.2, 32747.1,
#'                  22596.4, 13228.1, 17289.2, 30189.2, 31476.6, 15338.7),
#'                 .Dim = c(4L, 4L, 6L))
#' m1 <- c(1.025527, 1.018229, 0.969744, 0.994998)
#' m2 <- c(1.111023, 1.030213, 0.935041, 0.906709)
#' m3 <- c(0.810568, 1.375203, 1.07096, 1.044461, 0.949441, 0.915284)
#' m12 <- structure(c(1.061059, 1.120345, 1.097519, 1.188501, 1.017091,
#'                    0.967245, 1.03447, 1.18867, 0.9797, 0.900885, 0.85575, 1.070772,
#'                    1.041953, 1.074653, 0.887316, 0.791906), .Dim = c(4L, 4L))
#' m13 <- structure(c(0.779029, 0.865343, 0.757887, 0.852708, 1.351367,
#'                    1.409585, 1.350907, 1.361528, 1.091867, 1.107661, 0.99364, 1.127478,
#'                    1.13439, 0.948428, 1.075919, 0.916096, 1.031958, 0.835103, 1.006321,
#'                    0.982888, 0.86109, 0.976673, 0.961731, 0.764211), .Dim = c(4L, 6L))
#' m23 <- structure(c(0.962955, 0.880973, 0.798545, 0.714783, 1.547556,
#'                    1.277098, 1.149491, 1.210108, 1.186342, 1.084436, 0.976822, 1.003611,
#'                    1.092564, 1.066306, 1.038601, 0.996779, 0.971751, 1.016173, 0.867197,
#'                    0.803929, 0.831913, 0.933863, 0.857392, 0.960169), .Dim = c(4L, 6L))
#'
#' m <- list(m1, m3, m12, m23)
#' ind <- list(1, 3, c(1, 2), c(2, 3))
#'
#'
#' example <- WIPF(seed = s, weights = w, margins = m, indices = ind)
#'

WIPF <- function(seed,
                 weights,
                 margins,
                 indices,
                 normalize = TRUE,
                 tol = 10^-6,
                 maxit = 1000,
                 full = FALSE,
                 ...)
{

  # Inputs
  inputs <- c(as.list(environment()), list(...))

  # Initial checks
  if (!is.list(indices)) {
    stop("The argument 'indices' must be a list.")
  }
  if (!is.list(margins)) {
    stop("The argument 'margins' must be a list.")
  }
  if ( !(is.array(seed)) )
    stop("The argument 'seed' must be an array.")
  N <- length(dim(seed))
  test_indices(indices, N)
  mg.in <- test_indices_margins(margins, indices, N)
  margins <- mg.in$margins
  indices <- mg.in$indices

  # Sort and complete margins
  mg.in <- reorder_and_fill_margins(margins, indices, N)
  margins <- mg.in$margins
  indices <- mg.in$indices

  # Compatible margins
  margins <- test_WIPF(seed = seed,
                       weights = weights,
                       margins = margins,
                       indices = indices,
                       normalize = normalize,
                       tol = tol,
                       maxit = maxit)

  errors <- updating_errors_WIPF(seed = seed,
                                 weights = weights,
                                 margins = margins,
                                 indices = indices,
                                 normalize = normalize)

  dif <- max(errors$errors)
  iter <- 0L
  delta <- seed

  while(dif > tol & iter < maxit){
    delta <- updating_delta_WIPF(delta = delta,
                                 weights = weights,
                                 margins = margins,
                                 indices = indices,
                                 normalize = normalize)
    errors <- updating_errors_WIPF(seed = delta,
                                   weights = weights,
                                   margins = margins,
                                   indices = indices,
                                   normalize = normalize)
    dif <- max(errors$errors)
    iter <- iter + 1L
  }

  margin_names <- paste0("margin_", vapply(indices, paste0, character(1), collapse = "_"))

  names(margins) <- names(errors$dif) <- margin_names

  if(full){
    return(list("sol" = delta, "iter" = iter, "dev.margins" = errors$dif,
                "margins" = margins, "inputs" = inputs))
  } else {
    return(delta)
  }

}


