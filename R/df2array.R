#' Long-format data frame to array
#'
#' @description Organizes the information in a long-format data frame with N factor
#'              columns and one value column into a K-dimensional array (K \eqn{\le} N),
#'              where the size of each dimension corresponds to the number of levels
#'              of the respective factor.
#'
#' @author Jose M. Pavia, \email{pavia@@uv.es}
#'
#' @param df A long-format data frame.
#'
#' @param margins A vector of integers indicating which columns of `df` should be mapped to the
#'                dimensions of the array, and in what order. By default, `1:(ncol(df)-1)`,
#'                meaning all columns of `df` except the last are mapped to dimensions in the
#'                order they appear in the data frame.
#'
#' @param values An integer specifying which column of `df` is mapped to the cell values
#'               of the array. By default, `ncol(df)`, meaning that the values of the last column of `df`
#'               are used as the cell values.
#'
#' @param NA2zeros A `TRUE`/`FALSE` argument indicating whether intersections of levels
#'                 not present in `df` should be imputed with zero. Default, `TRUE`.
#'
#' @param names A `TRUE`/`FALSE` argument indicating whether the level labels of the factors
#'              should be used as dimension names (`dimnames`) for the array. Default is `TRUE`.
#'
#' @param ... Other arguments to be passed to the function. Not currently used.
#'
#' @return
#' An array with `length(margins)` dimensions.
#'
#' @export
#'
#' @examples
#'  x <- structure(list(LF1 = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
#'                              1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L,
#'                              2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
#'                              2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L,
#'                              3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L,
#'                              3L, 3L, 3L, 3L, 3L, 3L, 3L, 4L, 4L, 4L, 4L, 4L, 4L,
#'                              4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L,
#'                              4L, 4L, 4L, 4L, 4L),
#'                      LF2 = c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 3L,
#'                              3L, 3L, 3L, 3L, 3L, 4L, 4L, 4L, 4L, 4L, 4L, 1L, 1L,
#'                              1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L,
#'                              3L, 3L, 3L, 4L, 4L, 4L, 4L, 4L, 4L, 1L, 1L, 1L, 1L,
#'                              1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L,
#'                              3L, 4L, 4L, 4L, 4L, 4L, 4L, 1L, 1L, 1L, 1L, 1L, 1L,
#'                              2L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 3L, 4L,
#'                              4L, 4L, 4L, 4L, 4L),
#'                      LF3 = c(1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L,
#'                              2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L,
#'                              3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L,
#'                              4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L,
#'                              5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L,
#'                              6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L,
#'                              1L, 2L, 3L, 4L, 5L, 6L, 1L, 2L, 3L, 4L, 5L, 6L, 1L,
#'                              2L, 3L, 4L, 5L, 6L),
#'                     ix = c(0.8812, 0.8887, 1.0035, 0.8782, 1.1580, 1.3894, 0.7986,
#'                            1.0170, 1.0875, 0.9499, 0.9524, 1.1707, 0.4907, 1.4251,
#'                            0.8045, 0.7830, 0.7144, 0.9673, 0.5705, 6.8399, 0.6700,
#'                            0.6110, 2.1088, 0.7673, 0.8206, 1.0989, 1.0824, 0.7626,
#'                            1.1863, 1.6287, 0.8107, 0.8689, 1.0907, 0.9404, 0.9957,
#'                            1.2035, 0.5604, 0.9439, 0.8367, 0.7845, 0.8614, 1.0996,
#'                            0.3270, 1.1892, 0.6776, 0.5313, 0.7801, 0.9651, 1.2576,
#'                            1.1939, 1.2554, 1.1225, 1.5741, 1.5718, 0.8092, 0.8460,
#'                            1.0899, 1.0742, 1.0668, 1.0680, 0.8204, 0.8988, 1.0015,
#'                            1.0354, 0.9541, 1.0639, 0.5223, 0.6963, 0.6749, 0.7230,
#'                            0.6616, 0.9579, 1.1752, 1.3359, 1.2824, 1.6836, 1.5313,
#'                            2.5715, 1.0579, 0.8304, 1.0632, 1.0016, 0.9370, 1.1711,
#'                            0.7874, 1.4360, 1.0949, 0.8646, 0.8430, 1.4736, 0.7795,
#'                            0.9362, 0.8489, 0.8246, 0.8449, 0.6331)),
#'                            row.names = c(NA, -96L), class = "data.frame")
#'
#' example <- df2array(df = x, margins = c(2, 1, 3))
#'

df2array <- function(df,
                     margins = 1:(ncol(df)-1),
                     values = ncol(df),
                     NA2zeros = TRUE,
                     names = TRUE,
                     ...) {

  if ( !(is.data.frame(df)) )
    stop("The argument 'df' must be a data.frame.")
  if (max(abs(margins - round(margins))) > 0)
    stop("The argument 'margins' must be a vector of integers.")
  if (length(margins) < 2)
    stop("The argument 'margins' must be a vector of at least length 2.")
  if (abs(values - round(values)) > 0)
    stop("The argument 'values' must an integer.")
  if (ncol(df) < length(margins))
    stop("The length of argument 'margins' cannot be larger than the number of columns in 'df'.")
  if (max(margins, values) > ncol(df))
    stop("Not enough columns in 'df' given the values declared in 'margins' and/or 'values'.")
  if(any(duplicated(df)))
    stop("There are duplicated combinations of level of factors in the columns selected in 'df' to build the array.")

  niveles <- apply(df[margins], 2, unique)

  out <- array(NA, dim = sapply(niveles, length))
  if (names) dimnames(out) <- niveles

  combinaciones <- expand.grid(niveles)
  for (ii in 1L:nrow(combinaciones)) {
    fila <- combinaciones[ii, ]
    fila_coincidente <- which(apply(df[margins], 1L,
                                    function(x) all(x == fila)))
    if (length(fila_coincidente) > 0)
       out[as.matrix(fila)] <- df[values][fila_coincidente, ]
  }

  if (NA2zeros) out[is.na(out)] <- 0

  return(out)
}
