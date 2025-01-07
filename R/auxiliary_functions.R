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
    warning("The values of at least a couple the 1D margins are not compatible given the 'weights', they have been updated with WIPF1 in ascending order of sub-indices.")
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

