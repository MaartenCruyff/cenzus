softmax <- function(x)exp(c(0, x))/sum(exp(c(0, x)))

complete_LLM <- function(beta, D, n, dimp, f0){

  bx <- c(-log(sum(exp(D[, -1] %*% beta))), beta)

  p  <- exp(D %*% bx)/n

  pr <- Matrix::crossprod(p, dimp)

  -as.vector(Matrix::tcrossprod(f0, log(pr)))
}

complete_LL <- function(beta, D, n, dimp, f0){

  bx <- c(-log(sum(exp(D[, -1] %*% beta))), beta)

  p  <- exp(D %*% bx)/n

  pr <- crossprod(p, dimp)

  -tcrossprod(f0, log(pr))

}


loglin <- function(b0, f, D) {

  b   <- c(-log(sum(exp(D[, -1] %*% b0))) + log(sum(f)), b0)

  tol <- 1

  while (tol > 1e-2) {

    db   <- crossprod(D, f) - crossprod(D, exp(D %*% b))

    mu   <- diag(c(exp(D %*% b)))

    I    <- crossprod(D, mu) %*% D

    cvc  <- chol2inv(chol(I))

    bt   <- b + cvc %*% db

    tol  <- sum(abs(bt - b))

    b    <- bt
  }

  b
}


boot_fit <- function(d0, dx, n, N, D, n0, ds, fx, dimp, betas, model){


  bs <- data.frame(ds[, -ncol(ds)], Freq = rmultinom(1:nrow(ds),
                                                     size  = N,
                                                     prob  = ds$Freq/N))[-nrow(ds), ]


  f0   <- bs$Freq

  fxc  <- sum(f0) * fx / n


  ########################################################
  ### Dense matrix procedure if cells in dimp < sparse ###
  ########################################################

  if (class(dimp) == "matrix") {

    b      <- loglin(b0 = betas, f = fxc, D = D)      # M-step: parameter estimation

    fxc    <- exp(D %*% b)                                          # update observed frequencies

    p      <- c(fxc / sum(fxc))                                     # probabilities observed profiles

    num    <- t(p * dimp) * f0

    denom  <- crossprod(p, dimp)

    fxc    <- c(denom^-1 %*% num)                                   # E-step: update

    LL     <- log(denom) %*% f0                                     # loglikehood

    old  <- LL

    tol  <- 1

    while(tol > 1e-7){

      b      <- loglin(b0 = b[-1], f = fxc, D = D)

      fxc    <- exp(crossprod(t(D), b))

      p      <- c(fxc / sum(fxc))

      num    <- t(p * dimp) * f0

      denom  <- crossprod(p, dimp)

      fxc    <- c(denom^-1 %*% num )

      LL     <- log(denom) %*% f0

      tol <- LL - old

      old <- LL

    }

  }

  #######################################################
  ### Sparse matrix procedure if cells dimp > 10^5   ####
  #######################################################


  if (class(dimp) != "matrix") {

    b      <- loglin(b0 = rep(0, ncol(D) - 1), f = fxc, D = D)

    fxc    <- exp(D %*% b)

    p      <- c(fxc / sum(fxc))

    num    <- Matrix::t(p * dimp) * f0

    denom  <- Matrix::crossprod(p, dimp)

    fxc    <- as.vector(denom^-1 %*% num)

    LL     <- Matrix::tcrossprod(f0, log(denom))

    old  <- LL

    tol  <- 1

    while(tol > 1e-7){

      b      <- loglin(b0 = b[-1], f = fxc, D = D)

      fxc    <- exp(D %*% b)

      p      <- c(fxc / sum(fxc))

      num    <- Matrix::t(p * dimp) * f0

      denom  <- Matrix::crossprod(p, dimp)

      fxc    <- as.vector(denom^-1 %*% num)

      LL     <- Matrix::tcrossprod(f0, log(denom))

      tol <- as.vector(LL - old)

      old <- LL

    }


  }

  colnames(dx)[which(colnames(dx) == "fitted")] <- "Freq"

  c(exp(model.matrix(model, dx) %*% b))


}

