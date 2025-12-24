library(MASS)
library(expm)
library(rTensor)
library(tensr)
library(glasso)

Matrix.simulation.summary <- function(Omega.hat.list, Omega.true.list, 
                               norm.off.diag = FALSE,
                               crit = 1e-3) {
  
  if (!is.list(Omega.hat.list)) {
    stop("argument Omega.hat.list should be a list")
  }
  else if (!is.list(Omega.true.list)) {
    stop("argument Omega.true.list should be a list")
  }
  else if (any(!sapply(Omega.hat.list, is.matrix))) {
    stop("argument Omega.hat.list should be a list of precision matrices")
  }
  else if (any(!sapply(Omega.true.list, is.matrix))) {
    stop("argument Omega.true.list should be a list of precision matrices")
  }
  else if (length(Omega.hat.list) != length(Omega.true.list)) {
    stop("arguments Omega.hat.list and Omega.true.list should share the same length")
  }
  else if (any(!(sapply(Omega.hat.list, dim)[1, ] == sapply(
    Omega.true.list,
    dim
  )[1, ]))) {
    stop("dimension of elements in argument Omega.hat.list should match argument Omega.true.list")
  }
  else if (!is.logical(norm.off.diag)) {
    stop("argument norm.off.diag should be a logical TRUE or FALSE ")
  }
  else if(length(Omega.hat.list)!=2 | length(Omega.true.list)!=2){
    stop("This function needs to be run in Matrix distribution.")
  }
  
  K <- dim(as.array(Omega.hat.list))
  Est.Sigma <- solve(Omega.hat.list[[1]]);C <- diag(sqrt(diag(Est.Sigma)))
  scale.2nd <- norm(Omega.true.list[[2]],"F") / norm(Omega.hat.list[[2]],"F")
  
  # estimation error in Frobenius norm and Maximum norm
  error.f <- rep(0, K) 
  error.max <- rep(0, K)
  # true positive rate and true negative rate
  tpr <- rep(0, K)
  fpr <- rep(0, K)
  tnr <- rep(0, K)
  fnr <- rep(0, K)
  nonzero <- rep(0,K)
  
  scaled.precision <- array(list(),K)
  scaled.precision[[1]] <- C%*%Omega.hat.list[[1]]%*%C;scaled.precision[[2]] <- Omega.hat.list[[2]] * scale.2nd
  
  if(norm.off.diag == FALSE){
    for(mode in 1:K){
      error.f[mode] <- norm(scaled.precision[[mode]] - Omega.true.list[[mode]],"F") / norm(Omega.true.list[[mode]],"F")
      error.max[mode] <- norm(scaled.precision[[mode]] - Omega.true.list[[mode]],"M") / norm(Omega.true.list[[mode]],"M")
    }
  }
  if(norm.off.diag == TRUE){
    Omega.hat.list.off <- scaled.precision
    Omega.true.list.off <- Omega.true.list
    for(mode in 1:K){
      diag(Omega.hat.list.off[[mode]]) <- 0;diag(Omega.true.list.off[[mode]]) <- 0
      error.f[mode] <- norm(Omega.hat.list.off[[mode]] - Omega.true.list.off[[mode]],"F") / norm(Omega.true.list.off[[mode]],"F")
      error.max[mode] <- norm(Omega.hat.list[[mode]] - Omega.true.list.off[[mode]],"M") / norm(Omega.true.list.off[[mode]],"M")
    }
  }
  Omega.hat.list.off <- scaled.precision
  Omega.true.list.off <- Omega.true.list
  for(mode in 1:K){
    diag(Omega.hat.list.off[[mode]]) <- NA;diag(Omega.true.list.off[[mode]]) <- NA
    tpr[mode] <- length(intersect(which(abs(Omega.hat.list.off[[mode]]) > crit), 
                                  which(abs(Omega.true.list.off[[mode]]) > crit)))/ length(which(abs(Omega.true.list.off[[mode]]) > crit))
    fpr[mode] <- length(intersect(which(abs(Omega.hat.list.off[[mode]]) > crit),
                                  which(abs(Omega.true.list.off[[mode]]) <= crit))) / length(which(abs(Omega.true.list.off[[mode]]) <= crit))                                                                                                 
    tnr[mode] <- length(intersect(which(abs(Omega.hat.list.off[[mode]]) <= crit), 
                                  which(abs(Omega.true.list.off[[mode]]) <= crit))) / length(which(abs(Omega.true.list.off[[mode]]) <= crit))
    fnr[mode] <- length(intersect(which(abs(Omega.hat.list.off[[mode]]) <= crit), 
                                  which(abs(Omega.true.list.off[[mode]]) > crit))) / length(which(abs(Omega.true.list.off[[mode]]) > crit))
    B <- Omega.hat.list.off[[mode]];diag(B) <- 0
    nonzero[mode] <- sum(abs(B)> crit)
  }
  # output
  Out <- list()
  Out$error.f <- error.f
  Out$error.max <- error.max
  Out$tpr <- tpr
  Out$fpr <- fpr
  Out$tnr <- tnr
  Out$fnr <- fnr
  Out$nonzero <- nonzero
  return(Out)
}
