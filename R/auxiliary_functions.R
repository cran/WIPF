### Auxiliary functions

## test_WIPF3
test_WIPF3 <- function(seed,
                       weights,
                       margin1,
                       margin2,
                       margin3,
                       margin12,
                       margin13,
                       margin23,
                       normalize,
                       tol,
                       maxit,
                       full,
                       ...) {

  if (!identical(dim(seed), dim(weights)))
    stop('Arguments "seed" and "weights" must have the same dimensions')

  # checking non-negative values
  if (min(seed) < 0) {
    stop('Error: the object "seed" contains negative values.')
  }
  if (min(weights) < 0 | min(apply(weights, 1, max)) == 0 |
      min(apply(weights, 2, max)) == 0 | min(apply(weights, 3, max)) == 0 ) {
    stop('Error: non-negative values with at least a positive value in each row, column and layer is required in "weights".')
  }

  # Checking non-positive values in margins
  if (!is.null(margin1)) {
    if (min(margin1) <= 0)
      stop('Error: the object "margin1" contains non-positive values.')
  }
  if (!is.null(margin2)) {
    if (min(margin2) <= 0)
      stop('Error: the object "margin2" contains non-positive values.')
  }
  if (!is.null(margin3)) {
    if (min(margin3) <= 0)
      stop('Error: the object "margin3" contains non-positive values.')
  }
  if (!is.null(margin12)) {
    if (min(margin12) < 0)
      stop('Error: the object "margin12" contains negative values.')
  }
  if (!is.null(margin13)) {
    if (min(margin13) < 0)
      stop('Error: the object "margin13" contains negative values.')
  }
  if (!is.null(margin23)) {
    if (min(margin23) < 0)
      stop('Error: the object "margin23" contains negative values.')
  }

  # checking dimensions
  if(length(dim(seed)) != 3){
    stop('Error: the object "seed" must be a 3D array.')
  }
  if(length(dim(weights)) != 3){
    stop('Error: the object "weights" must be a 3D array.')
  }
  if(!all(dim(seed) == dim(weights))){
    stop('Error: the objects "seed" and weights" must have the same dimensions.')
  }
  if(length(margin1) !=0){
    if(dim(seed)[1L] != length(margin1))
      stop('Error: the dimensions of "seed" and margin1" are incompatible.')
  }
  if(length(margin2) !=0){
    if(dim(seed)[2L] != length(margin2))
      stop('Error: the dimensions of "seed" and margin2" are incompatible.')
  }
  if(length(margin3) !=0){
    if(dim(seed)[3L] != length(margin3))
      stop('Error: the dimensions of "seed" and margin3" are incompatible.')
  }
  if(length(margin12) !=0){
    if(!all(dim(seed)[c(1L, 2L)] == dim(margin12)))
      stop('Error: the dimensions of "seed" and margin12" are incompatible.')
  }
  if(length(margin13) !=0){
    if(!all(dim(seed)[c(1L, 3L)] == dim(margin13)))
      stop('Error: the dimensions of "seed" and margin13" are incompatible.')
  }
  if(length(margin23) !=0){
    if(!all(dim(seed)[c(2L, 3L)] == dim(margin23)))
      stop('Error: the dimensions of "seed" and margin23" are incompatible.')
  }

  # Compatibility of margins
  # Compatibility 1D
  if(normalize){
    weights12 <- apply(weights, c(1L, 2L), sum)
    w12 <- weights12/sum(weights12)
    pesos1 <- rowSums(w12)
    pesos2 <- colSums(w12)
    weights13 <- apply(weights, c(1L, 3L), sum)
    w13 <- weights13/sum(weights13)
    pesos3 <- colSums(w13)
  } else {
    pesos1 <- pesos2 <- pesos3 <- 1
  }

  if(!is.null(margin1) & !is.null(margin2) |
     !is.null(margin1) & !is.null(margin3) |
     !is.null(margin2) & !is.null(margin3) ){
    if (!is.null(margin1) & !is.null(margin2)){
      d12 <- abs(sum(margin1*pesos1) - sum(margin2*pesos2))
    } else {
      d12 <- 0
    }
    if (!is.null(margin1) & !is.null(margin3)){
      d13 <- abs(sum(margin1*pesos1) - sum(margin3*pesos3))
    } else {
      d13 <- 0
    }
    if (!is.null(margin2) & !is.null(margin3)){
      d23 <- abs(sum(margin2*pesos2) - sum(margin3*pesos3))
    } else {
      d23 <- 0
    }
  } else {
    d12 <- d13 <- d23 <- 0
  }

  if (max(c(d12, d13, d23)) > tol){
    warning("The values of at least a couple of 1D margins are not compatible given the 'weights',
            they have been updated with WIPF1 in ascending order of sub-indices.")
    if(!is.null(margin1) & !is.null(margin2))
      margin2 <- WIPF1(seed = margin2, weights = apply(weights, 2L, sum),
                       margin = sum(margin1*pesos1), normalize = normalize,
                       tol = tol, maxit = maxit)
    if(!is.null(margin1) & !is.null(margin3))
      margin3 <- WIPF1(seed = margin3, weights = apply(weights, 3L, sum),
                       margin = sum(margin1*pesos1), normalize = normalize,
                       tol = tol, maxit = maxit)
    if(!is.null(margin2) & !is.null(margin3))
      margin3 <- WIPF1(seed = margin3, weights = apply(weights, 3L, sum),
                       margin = sum(margin2*pesos2), normalize = normalize,
                       tol = tol, maxit = maxit)
  }

  # Compatibility 1D and 2D
  if(normalize){
    weights12 <- apply(weights, c(1L, 2L), sum)
    weights13 <- apply(weights, c(1L, 3L), sum)
    weights23 <- apply(weights, c(2L, 3L), sum)
    w12 <- weights12/sum(weights12)
    w13 <- weights13/sum(weights13)
    w23 <- weights23/sum(weights23)
    pesos1 <- rowSums(w12)
    pesos2 <- colSums(w12)
    pesos3 <- colSums(w13)
  } else {
    w12 <- array(1, dim = dim(weights)[c(1L, 2L)])
    w13 <- array(1, dim = dim(weights)[c(1L, 3L)])
    w23 <- array(1, dim = dim(weights)[c(2L, 3L)])
    pesos1 <- pesos2 <- pesos3 <- 1
  }

  if(!is.null(margin12)){
    if (!is.null(margin1)){
      d121 <- max(abs(rowSums(margin12 * w12) - margin1 * pesos1))
    } else {
      d121 <- 0
    }
    if (!is.null(margin2)){
      d122 <- max(abs(colSums(margin12 * w12) - margin2 * pesos2))
    } else {
      d122 <- 0
    }
  } else {
    d122 <- d121 <- 0
  }
  if(max(d121, d122) > tol){
      margin12 <- WIPF2(seed = margin12, weights = apply(weights, c(1L,2L), sum),
                        margin1 = margin1, margin2 = margin2,
                        normalize = normalize, tol = tol, maxit = maxit)
  }

  # Dimensions 1 and 3
  if(!is.null(margin13)){
    if (!is.null(margin1)){
      d131 <- max(abs(rowSums(margin13 * w13) - margin1 * pesos1))
    } else {
      d131 <- 0
    }
    if (!is.null(margin3)){
      d133 <- max(abs(colSums(margin13 * w13) - margin3 * pesos3))
    } else {
      d133 <- 0
    }
  } else {
    d131 <- d133 <- 0
  }
  if(max(d131, d133) > tol){
      margin13 <- WIPF2(seed = margin13, weights = apply(weights, c(1L,3L), sum),
                        margin1 = margin1, margin2 = margin3,
                        normalize = normalize, tol = tol, maxit = maxit)
  }

  # Dimensions 2 and 3
  if(!is.null(margin23)){
    if (!is.null(margin2)){
      d232 <- max(abs(rowSums(margin23 * w23) - margin2 * pesos2))
    } else {
      d232 <- 0
    }
    if (!is.null(margin3)){
      d233 <- max(abs(colSums(margin23 * w23) - margin3 * pesos3))
    } else {
      d233 <- 0
    }
  } else {
    d232 <- d233 <- 0
  }
  if(max(d232, d233) > tol){
    margin23 <- WIPF2(seed = margin23, weights = apply(weights, c(2L,3L), sum),
                      margin1 = margin2, margin2 = margin3,
                      normalize = normalize, tol = tol, maxit = maxit)
  }

  if (max(d121, d122, d131, d133, d232, d233) > tol)
    warning("The values of at least a 2D margin and a 1D margin are not compatible given the 'weights', they have been updated using WIPF2 in ascending order of sub-indices.")

  return(list("margin1" = margin1, "margin2" = margin2, "margin3" = margin3,
              "margin12" = margin12, "margin13" = margin13, "margin23" = margin23))

}
## -----------------------------------------------------------------------
## updating errors
updating_errors <- function(seed, weights, margins, normalize){
  margin1 <- margins$margin1
  margin2 <- margins$margin2
  margin3 <- margins$margin3
  margin12 <- margins$margin12
  margin13 <- margins$margin13
  margin23 <- margins$margin23

  if(normalize){
    weights1 <- apply(weights, 1L, sum)
    weights2 <- apply(weights, 2L, sum)
    weights3 <- apply(weights, 3L, sum)
    weights12 <- apply(weights, c(1L, 2L), sum)
    weights13 <- apply(weights, c(1L, 3L), sum)
    weights23 <- apply(weights, c(2L, 3L), sum)
  } else {
    weights12 <- weights13 <- weights23 <- weights1 <- weights2 <- weights3 <- 1
  }

  # errors
  delta <- seed
  if (!is.null(margin1)){
    error1 <- abs(apply(weights * delta, 1L, sum) - weights1 * margin1)/weights1
  } else {
    error1 <- 0
  }
  if (!is.null(margin2)){
    error2 <- abs(apply(weights * delta, 2L, sum) - weights2 * margin2)/weights2
  } else {
    error2 <- 0
  }
  if (!is.null(margin3)){
    error3 <- abs(apply(weights * delta, 3L, sum) - weights3 * margin3)/weights3
  } else {
    error3 <- 0
  }
  if (!is.null(margin12)){
    error12 <- abs(apply(weights * delta, c(1L, 2L), sum) - weights12 * margin12)/weights12
  } else {
    error12 <- 0
  }
  if (!is.null(margin13)){
    error13 <- abs(apply(weights * delta, c(1L, 3L), sum) - weights13 * margin13)/weights13
  } else {
    error13 <- 0
  }
  if (!is.null(margin23)){
    error23 <- abs(apply(weights * delta, c(2L, 3L), sum) - weights23 * margin23)/weights23
  } else {
    error23 <- 0
  }

  dif <- list("dev1" = error1, "dev2" = error2, "dev3" = error3,
              "dev12" = error12, "dev13" = error13, "dev23" = error23)
  dif <- dif[sapply(dif, length) != 1]

  errors <- c(max(error1), max(error2), max(error3), max(error12),
              max(error13), max(error23))
  return(list("errors" = errors, "dif" = dif))
}
## -----------------------------------------------------------------------
## updating delta
updating_delta <- function(delta, weights, margins, normalize){
  margin1 <- margins$margin1
  margin2 <- margins$margin2
  margin3 <- margins$margin3
  margin12 <- margins$margin12
  margin13 <- margins$margin13
  margin23 <- margins$margin23

  if(normalize){
    weights1 <- apply(weights, 1L, sum)
    weights2 <- apply(weights, 2L, sum)
    weights3 <- apply(weights, 3L, sum)
    weights12 <- apply(weights, c(1L, 2L), sum)
    weights13 <- apply(weights, c(1L, 3L), sum)
    weights23 <- apply(weights, c(2L, 3L), sum)
  } else {
    weights12 <- weights13 <- weights23 <- weights1 <- weights2 <- weights3 <- 1
  }

  tamanyo <- dim(weights)
  # updating_process
  if (!is.null(margin1)){
    factor1 <- (weights1 * margin1)/apply(weights * delta, 1L, sum)
    factor1[is.nan(factor1)] <- 0
    factor1[is.infinite(factor1)] <- 1
    factor1 <- replicate(tamanyo[3L],
                         matrix(factor1, nrow = tamanyo[1L], ncol = tamanyo[2L]),
                         simplify = "array")
    delta <- delta * factor1
  }
  if (!is.null(margin2)){
    factor2 <- (weights2 * margin2)/apply(weights * delta, 2L, sum)
    factor2[is.nan(factor2)] <- 0
    factor2[is.infinite(factor2)] <- 1
    factor2 <- aperm(replicate(tamanyo[3L],
                               matrix(factor2, nrow = tamanyo[1L], ncol = tamanyo[2L]),
                               simplify = "array"), c(2L, 1L, 3L))
    delta <- delta * factor2
  }
  if (!is.null(margin3)){
    factor3 <- (weights3 * margin3)/apply(weights * delta, 3L, sum)
    factor3[is.nan(factor3)] <- 0
    factor3[is.infinite(factor3)] <- 1
    factor3 <- aperm(replicate(tamanyo[2L],
                               matrix(factor3, nrow = tamanyo[3L], ncol = tamanyo[1L]),
                               simplify = "array"), c(2L, 3L, 1L))
    delta <- delta * factor3
  }
  if (!is.null(margin12)){
    factor12 <- (weights12 * margin12)/apply(weights * delta, c(1L, 2L), sum)
    factor12[is.nan(factor12)] <- 0
    factor12[is.infinite(factor12)] <- 1
    factor12 <- replicate(tamanyo[3L], factor12, simplify = "array")
    delta <- delta * factor12
  }
  if (!is.null(margin13)){
    factor13 <- (weights13 * margin13)/apply(weights * delta, c(1L, 3L), sum)
    factor13[is.nan(factor13)] <- 0
    factor13[is.infinite(factor13)] <- 1
    factor13 <- aperm(replicate(tamanyo[2L], factor13, simplify = "array"), c(1L, 3L, 2L))
    delta <- delta * factor13
  }
  if (!is.null(margin23)){
    factor23 <- (weights23 * margin23)/apply(weights * delta, c(2L, 3L), sum)
    factor23[is.nan(factor23)] <- 0
    factor23[is.infinite(factor23)] <- 1
    factor23 <- aperm(replicate(tamanyo[1L], factor23, simplify = "array"), c(3L, 1L, 2L))
    delta <- delta * factor23
  }

  return(delta)
}

###.......................................................................................
############
##  WIPF  ##
############

## test indices
test_indices <- function(indices, N) {

  for (i in seq_along(indices)) {
    idx <- indices[[i]]

    if (!is.vector(idx)) {
      stop("An element %d of argument `indices` is not a vector.")
    }

    if (!is.numeric(idx)) {
      stop("An element of argument `indices` is not numeric.")
    }

    if (any(idx %% 1 != 0)) {
      stop("An element of argument `indices` contains non-integer values.")
    }

    if (any(idx < 1 | idx > N)) {
      stop("An element of argument `indices` contains values outside [1, Dimension of seed].")
    }

    if (any(duplicated(idx))) {
      stop("An element of argument `indices` contains duplicated indices.")
    }
  }

  # Check de márgenes no repetidos, incluso tras permutación
  indices_norm <- vapply(
    indices,
    function(x) paste(sort(x), collapse = ","),
    character(1)
  )

  if (any(duplicated(indices_norm))) {
    dup <- unique(indices_norm[duplicated(indices_norm)])
    stop("Repeated index sets detected in argument 'indices'.")
  }

  invisible(TRUE)
}

#-------------------------------------------------------------------------------
## test indices and margins and order of margins
test_indices_margins <- function(margins, indices, N) {

  if (length(margins) != length(indices)) {
    stop("Length of arguments 'margins' and 'indices' must be the same.")
  }

  if (!all(vapply(margins,
                  function(x) is.null(x) || is.vector(x) || is.array(x),
                  logical(1)))) {
    stop("All elements of argument 'margins' must be vectors, matrices, or arrays.")
  }

  # Lista de márgenes reordenados
  margins_perm <- vector("list", length(margins))

  for (i in seq_along(margins)) {

    idx <- indices[[i]]
    margin <- margins[[i]]

    # si es vector, convertir a array 1D
    if (is.vector(margin) && is.null(dim(margin))) {
      margin <- array(margin, dim = length(margin))
    }

    # Chequear que la dimensión del margen coincide con la longitud del índice
    if (length(dim(margin)) != length(idx) & !is.null(margin)) {
      stop("Dimension mismatch between elements of arguments 'margings' and 'indices'.")
    }

    # Normalizamos el orden de los índices (de menor a mayor)
    ord <- order(idx)
    idx_sorted <- idx[ord]

    # Permutamos el array si es necesario
    if (!all(ord == seq_along(ord))) {
      margin <- aperm(margin, perm = ord)
    }

    margins_perm[[i]] <- margin
    indices[[i]] <- idx_sorted
  }

  # Márgenes reordenados y sus índices correspondientes
  return( list(margins = margins_perm, indices = indices) )
}

#-------------------------------------------------------------------------------
# Auxiliary function to generate a list of all possible combinations of indices
canonical_indices <- function(N) {
  out <- list()
  for (k in seq_len(N)) {
    combs <- utils::combn(N, k, simplify = FALSE)
    out <- c(out, combs)
  }
  return(out[-length(out)])
}

#-------------------------------------------------------------------------------
## Reordenamos márgenes (e índices) y completamos con márgenes NULL
reorder_and_fill_margins <- function(margins, indices, N) {

  # 1. Generamos todas las combinaciones canónicas de índices
  canon <- canonical_indices(N)

  # 2. Normalizamos los índices de entrada (orden ascendente)
  indices_norm <- lapply(indices, function(x) sort(as.integer(x)))

  # 3. Crear las listas finales de márgenes e índices
  margins_out <- vector("list", length(canon))
  indices_out <- vector("list", length(canon))

  for (i in seq_along(canon)) {
    idx <- canon[[i]]

    # Buscar si existe un margen correspondiente
    pos <- which(vapply(indices_norm, function(x) identical(x, idx), logical(1)))

    if (length(pos) == 1) {
      margins_out[[i]] <- margins[[pos]]
      indices_out[[i]] <- idx
    } else {
      # margins_out[i]] <- NULL
      indices_out[[i]] <- idx
    }
  }

  return(list( margins = margins_out, indices = indices_out) )
}

#-------------------------------------------------------------------------------
## Auxiliar function of check_and_update_margins_WIPF
### Function to reindex a list of indices
reindex_indices <- function(indices, idx) {
  lapply(indices, function(v) {
    match(v, idx)
  })
}

#-------------------------------------------------------------------------------
## Function to check the compatibility of margins and updated them if necessary
check_and_update_margins_WIPF <- function(seed, weights, margins, indices,
                                           tol = 1e-8,
                                           normalize = TRUE, maxit = 1000) {

  N <- length(dim(seed))
  canon <- canonical_indices(N)
  pesos <- margins
  for (i in 1L:length(canon)){
    if (normalize){
      w <- apply(weights, canon[[i]], sum)
      pesos[[i]] <- w/sum(w)
    } else {
      pesos[[i]] <- 1
    }
  }

  # Dimensiones 1
  if (N > 1L){
    sel <- sapply(canon, length) == 1L & !sapply(margins, is.null)
    if (sum(sel) > 1){
      pos <- which(sel)
      ref <- min(pos) # Referencia
      disc <- 0 # discrepancies
      for (ii in pos[-1L]){
        error <- abs(sum(margins[[ref]] * pesos[[ref]]) - sum(margins[[ii]] * pesos[[ii]]))
        disc <- c(disc, error)
        if (error > tol)
          margins[[ii]] <- WIPF1(seed = as.vector(margins[[ii]]),
                                 weights = apply(weights, ii, sum),
                                 margin = sum(margins[[ref]] * pesos[[ref]]),
                                 normalize = normalize,
                                 tol = tol, maxit = maxit)
      }
    }
  }

  # Dimensiones > 1
  if (N > 2L){
    for (d in 2L:(N - 1L)){ # iteramos dimensiones de orden 1 a N-1
      sel <- sapply(canon, length) == d & !sapply(margins, is.null)
      pos <- which(sel)
      for (ii in pos){
        disc2 <- NULL
        idx <- indices[[ii]] # Margin objetivo
        subidx <- unlist(lapply(1:(length(idx) - 1),
                                function(k) utils::combn(idx, k, simplify = FALSE)),
                         recursive = FALSE) # Submárgenes
        wsel <- NULL
        for (sii in 1L:length(subidx)){
          sel <- which(vapply(indices, identical, logical(1), subidx[[sii]]))
          wsel <- c(wsel, sel)
          if (!is.null(margins[[sel]])){
            error <- abs(sum(margins[[ii]] * pesos[[ii]]) - sum(margins[[sel]] * pesos[[sel]]))
            disc2 <- c(disc2, error)
          }
        }
        if (max(disc2) > tol){
          margins_t <- margins[wsel]
          indices_t <- reindex_indices(indices[wsel], idx)
          indices_t <- indices_t[!sapply(margins_t, is.null)]
          margins_t <- margins_t[!sapply(margins_t, is.null)]
          margins[[ii]] <- WIPF(seed = margins[[ii]],
                                weights = apply(weights, idx, sum),
                                margins = margins_t,
                                indices = indices_t, ##
                                tol = tol,
                                normalize = normalize,
                                maxit = maxit)
        } # if max(disc) > tol
        disc <- c(disc, disc2)
      } # for (ii in pos)
    } # for (d in 2:(N - 1))
  } # if (N > 2)

  if (max(disc) > tol)
    warning("The values of at least a couple of margins are not compatible given the 'weights',
              they have been updated in ascending order of sub-indices.")

  return(margins)
}

#-------------------------------------------------------------------------------
## test_WIPF
test_WIPF <- function(seed,
                      weights,
                      margins,
                      indices,
                      normalize,
                      tol,
                      maxit) {

  if (!identical(dim(seed), dim(weights)))
    stop('Error: arguments "seed" and "weights" must have the same dimensions')

  # checking non-negative values
  if (min(seed) < 0) {
    stop('Error: the argument "seed" contains negative values.')
  }

  if (min(weights) < 0) {
    stop('Error: the argument "weights" contains negative values.')
  }

  # At least a positive value in each dimension
  for (d in seq_along(dim(weights))) {
    margin_max <- apply(weights, MARGIN = d, max)
    if (min(margin_max) <= 0) {
      stop(sprintf(
        'Error: Argument "weights" must have at least one positive value in each slice along dimension %d.',
        d
      ))
    }
  }

  # Checking non-positive values in margins
  for (i in 1L:length(margins)){
    margin <- margins[[i]]
    if (!is.null(margin)) {
      if (min(margin) < 0)
        stop('Error: At least a component of "margins" contains negative values.')
    }
  }

  # Checking for compatibility between margins and seed
  for (i in seq_along(margins)) {
    margin <- margins[[i]]
    idx <- indices[[i]]
    if (!is.null(margin)) {
      if (!identical(dim(seed)[idx], dim(margin))) {
        stop('Error: Dimension mismatch between "seed" and one of the elements of "margins".')
      }
    }
  }

  # Compatibility of margins
  margins <- check_and_update_margins_WIPF(seed = seed,
                                           weights = weights,
                                           margins = margins,
                                           indices = indices,
                                           tol = tol,
                                           normalize = normalize,
                                           maxit = maxit)

  return(margins)
}

## -----------------------------------------------------------------------
## updating errors WIPF
updating_errors_WIPF <- function(seed, weights, margins, indices, normalize){

  # Pesos no ponderados
  N <- length(dim(seed))
  canon <- canonical_indices(N)
  pesos <- margins
  for (i in 1L:length(canon)){
    if (normalize){
      pesos[[i]] <- apply(weights, canon[[i]], sum)
    } else {
      pesos[[i]] <- 1
    }
  }

  # errors
  dif <- margins
  errors <- NULL
  delta <- seed

  for (ii in 1:length(margins)){
    margin <- margins[[ii]]
    peso <- pesos[[ii]]
    if (!is.null(margin)){
      dif[[ii]] <- abs(apply(weights * delta, indices[[ii]], sum) - peso * margin)/peso
      error <- max(dif[[ii]])
    } else {
      error <- 0
    }
    errors <- c(errors, error)
  }

  return(list("errors" = errors, "dif" = dif))
}

## ----------------------------------------------------------------------------
## updating delta WIPF
updating_delta_WIPF <- function(delta,
                                weights,
                                margins,
                                indices,
                                normalize) {

  # Pesos no ponderados
  N <- length(dim(delta))
  canon <- canonical_indices(N)
  pesos <- margins
  for (i in 1L:length(canon)){
    if (normalize){
      pesos[[i]] <- apply(weights, canon[[i]], sum)
    } else {
      pesos[[i]] <- 1
    }
  }

  tamanyo <- dim(weights)

  for (ii in seq_along(margins)) {
    margin <- margins[[ii]]
    if (!is.null(margin)){
      dims <- indices[[ii]]
      peso <- pesos[[ii]]
      denom <- apply(weights * delta, dims, sum)
      factor <- (peso * margin) / denom
      factor[is.nan(factor)] <- 0
      factor[is.infinite(factor)] <- 1

      delta <- sweep(delta, dims, factor, `*`)
    }
  }

  return(delta)
}
