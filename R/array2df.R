#' Array to long-format data frame
#'
#' @description Organizes the information in an array with N dimensions into a
#'              long-format data frame with N + 1 columns, where the last column
#'              contains the values of the array.
#'
#' @author Jose M. Pavia, \email{pavia@@uv.es}
#'
#' @param arr An array of N-dimensiones
#'
#' @param value_name Name for the column containing the values of the array.
#'

array2df <- function(arr, value_name = "value") {
  if ( !(is.array(arr)) )
    stop("The argument 'arr' must be an array.")

  dims <- dim(arr)
  dim_names <- dimnames(arr)
  n_dims <- length(dims)

  # Si no hay nombres de dimensiones, crear índices
  if (is.null(dim_names)) {
    dim_names <- lapply(dims, seq_len)
    names(dim_names) <- paste0("dim", seq_len(n_dims))
  } else {
    # Rellenar nombres faltantes
    for (i in seq_len(n_dims)) {
      if (is.null(dim_names[[i]])) {
        dim_names[[i]] <- seq_len(dims[i])
      }
      if (is.null(names(dim_names)[i]) || names(dim_names)[i] == "") {
        names(dim_names)[i] <- paste0("dim", i)
      }
    }
  }

  # Expandir combinaciones
  out <- expand.grid(dim_names, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

  # Añadir valores del array
  out[[value_name]] <- as.vector(arr)

  return(out)
}
