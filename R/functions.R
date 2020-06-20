#' Multiple Systems Estimation
#'
#' @description \code{mse} fits a loglinear model with the EM algorithm.
#' @param model formula of the model to be fitted. A register should be denoted by
#' a single upper case letter, a covariate by single lower case letters. Latent variables
#' are denoted by single upper case letters, as specified in \code{lat}.
#' @param dframe a data frame with the data, either as individual
#' records or as a contingency table. The registers should be coded 0 if not observed,
#' and 1 if observed. It is also possible to use a subset of the variables in the
#' data frame.
#' @param lat a character (vector) with the letter(s) denoting the latent variable(s) in
#' \code{model}. Defaults to \code{NULL} for the absence of latent variables.
#' @param nclass a numeric vector with the number of classes of the latent variables
#' specified in \code{lat}, if any.
#' @param crit the convergence criterion of the EM algorithm computed as the difference
#' of the complete data log-likelihood at ultimate and penultimate iteration. Defaults to 1e-7.
#' @param maxiter the maximum number of EM iterations. Defaults to 5000.
#' @param seed random seed determining the starting values of the EM algorithm. Defaults to 1.
#' Changing the seed may be helpful when the Hessian is non-invertible. See 'Details'.
#' @return A list with fitted objects (the first four are also printed to the screen, the others are
#' required by the \code{\link{boot_mse}} function).
#' \item{loglike}{the maximized incomplete data log-likelihood.}
#' \item{coefs}{a matrix with the log-linear parameter estimates, SE's, and t- and p-values.}
#' \item{probs}{a matrix with the fitted population probabilities of the variable levels.}
#' \item{fitted}{a matrix with variables and fitted frequencies.}
#' \item{hist}{iteration history.}
#' \item{model}{the model formula.}
#' \item{nclass}{a with the number(s) of classes of the latent variable(s), if any.}
#' \item{dobsmis}{contingency table excluding structural zeros and including missings.}
#' \item{dimp}{a matrix used for imputing the missings.}
#' \item{D}{the design matrix of the log-linear model.}
#' \item{fxc}{starting values for the bootstrap function.}
#' @details The Hessian may not be invertible when the model is overparametrized or when the
#' starting values for the EM algorithm are poorly chosen. In the latter case, changing the seed
#' may result in better starting values and consequently in an invertible Hessian.
#' @examples
#' ## Fit the maximal model for registers A and B and covariates a and b
#' AB <- mse(Freq ~ A*b + B*a + a*b, dframe = NewZealand[, c(1:2, 5:6, 9)])
#'
#' ## Fit a model with two latent classes for a, b, c and d.
#' mse(Freq ~ A + B + C + D + X*(a + b + c + d), dframe = NewZealand, lat = "X", nclass = 2)
#' @importFrom stats complete.cases model.matrix na.omit optim pnorm xtabs rnorm rmultinom
#' @export


mse <- function(model, dframe, lat = NULL, nclass = 1, crit = 1e-7,
                maxiter = 5000, seed = 1){


  set.seed(seed)

  mse <- ifelse(any(colnames(dframe) %in% LETTERS), TRUE, FALSE)      # if registers in data, mse is TRUE

  if ("Freq" %in% colnames(dframe)){                                  # input is contingency table

    d0 <- as.data.frame(xtabs(Freq ~ ., dframe, addNA = TRUE))   # contingency table population including missing labels

  } else {                                                            # input are indivdual records

    d0 <- as.data.frame(xtabs(     ~ ., dframe, addNA = TRUE))   # contingency table population including missing labels

  }

  g  <- function(x) eval(substitute(x))                           # help function for creating logical expression str0

  if (mse){

    lists <- colnames(d0)[colnames(d0) %in% LETTERS]

    if (!all(sapply(d0[, lists], FUN = levels) %in% 0:1)) stop("only values 0 and 1 allowed for lists")

    str0 <- NULL

    for(i in 1:length(lists)){

      if(i == length(lists)){

        str0 <- paste(str0, g(lists[i]), "== 0 ")

      } else {

        str0 <- paste(str0, g(lists[i]), "== 0 &")                # logical test for structural zeros
      }
    }

    dobsmis  <- subset(d0, !eval(parse(text = str0)))             # remove rows with structural zeros

  } else {

    dobsmis <- d0

  }


  ########################
  ## initialize objects ##
  ########################

  f0   <- dobsmis$Freq

  n    <- sum(f0)

  dobs <- dobsmis[complete.cases(dobsmis),
                  colnames(dobsmis) != "Freq"]                     # contingency table excluding missing (and structural zeros)

  ########################################
  ## add latent class variables to data ##
  ########################################

  dx   <- dobsmis                                                  # contingency table including latent variables and missings

  if (!is.null(lat)){

    dtmp <- NULL

    for(j in length(lat):1){

      for(i in 1:nclass[j]){

        dtmp <- rbind(dtmp, dx)
      }

      dtmp <- data.frame(factor(rep(1:nclass[j], each = nrow(dx))), dtmp)

      names(dtmp)[1] <- g(lat[j])

      dx <- dtmp

      dtmp <- NULL
    }
  }

  dxc  <- dx[complete.cases(dx), ]                                  # profiles with latent variables excluding missings

  D    <- model.matrix(model, dxc)                                  # design matrix


  imp  <- matrix(1, nrow(dobs), nrow(dobsmis))                       # matrix for imputation missings without latent variables


  for(k in 1:nrow(dobsmis)){

    for(j in 1:ncol(dobs)){

      if(!is.na(dobsmis[k, j])){

        imp[, k] <- imp[, k] * (dobs[, j] == dobsmis[k, j])

      }

    }

  }


  dimp <- NULL                                                       # matrix for imputation missing with latent variables

  if(is.null(lat)){

    dimp <- imp

  } else {

    for(i in 1:prod(nclass)){

      dimp <- rbind(dimp, imp)

    }
  }

  ########################################################
  ### Dense matrix procedure if cells in dimp < sparse ###
  ########################################################

  if (prod(dim(dimp)) <= 1e5) {

    fxc    <- n * softmax(rep(0, nrow(D) - 1))                      # Uniform starting frequencies

    b      <- loglin(b0 = rep(0, ncol(D) - 1), f = fxc, D = D)       # M-step: parameter estimation

    fxc    <- exp(D %*% b)                                          # update observed frequencies

    p      <- c(fxc / sum(fxc))                                     # probabilities observed profiles

    num    <- t(p * dimp) * f0

    denom  <- crossprod(p, dimp)

    fxc    <- c(denom^-1 %*% num)                                   # E-step: update

    LL     <- log(denom) %*% f0                                     # loglikehood

    hist   <- data.frame(loglike   = rep(NA, maxiter),
                         tolerance = rep(NA, maxiter))              # iteration history

    hist[1, ] <- c(LL, -LL)

    old  <- LL

    tol  <- 1

    iter <- 1

    while(tol > crit){

      b      <- loglin(b0 = b[-1], f = fxc, D = D)

      fxc    <- exp(D %*% b)

      p      <- c(fxc/sum(fxc))

      num    <- t(p * dimp) * f0

      denom  <- crossprod(p, dimp)

      fxc    <- c(denom^-1 %*% num )

      LL     <- log(denom) %*% f0

      hist[iter + 1, ] <- c(LL, round(LL - old, 8))

      tol <- LL - old

      old <- LL

      iter <- iter + 1

      if (iter == maxiter) break

    }

    clgl <- optim(par  = b[-1], hessian = T,
                  fn   = complete_LL,
                  D    = D,
                  n    = n,
                  dimp = dimp,
                  f0   = f0)
  }

  #######################################################
  ### Sparse matrix procedure if cells dimp > sparse ####
  #######################################################


  if (prod(dim(dimp)) > 1e5) {

    dimp   <- Matrix::Matrix(dimp)

    fxc    <- n * softmax(rnorm(nrow(dxc) - 1))

    b      <- loglin(b0 = rep(0, ncol(D) - 1), f = fxc, D = D)

    fxc    <- exp(D %*% b)

    p      <- c(fxc / sum(fxc))

    num    <- Matrix::t(p * dimp) * f0

    denom  <- Matrix::crossprod(p, dimp)

    fxc    <- as.vector(denom^-1 %*% num)

    LL     <- Matrix::tcrossprod(f0, log(denom))

    hist   <- data.frame(loglike = rep(NA, maxiter), tol = rep(NA, maxiter))

    hist[1, ] <- c(LL, -LL)

    old  <- LL

    tol  <- 1

    iter <- 1

    while(tol > crit){

      b      <- loglin(b0 = b[-1], f = fxc, D = D)

      fxc    <- exp(D %*% b)

      p      <- c(fxc / sum(fxc))

      num    <- Matrix::t(p * dimp) * f0

      denom  <- Matrix::crossprod(p, dimp)

      fxc    <- as.vector(denom^-1 %*% num)

      LL     <- Matrix::tcrossprod(f0, log(denom))

      hist[iter + 1, ] <- c(LL, round(LL - old, 8))

      tol <- as.vector(LL - old)

      old <- LL

      iter <- iter + 1

      if (iter == maxiter) break

    }

    clgl <- optim(par  = b[-1, 1, drop = F], hessian = T, method =  "BFGS",
                  fn   = complete_LLM,
                  D    = D,
                  n    = n,
                  dimp = dimp,
                  f0   = f0)
  }

  hist  <- na.omit(hist)

  se    <- c(NA, sqrt(diag(solve(clgl$hessian))))

  eig   <- eigen(clgl$hessian, only.values = TRUE)$values

  coefs <- data.frame(beta = b, se,
                      zval = b / se,
                      pval = 2 * pnorm(-abs(b / se)), row.names = colnames(D))


  dpop   <-  as.data.frame(xtabs(Freq ~ ., dxc))

  fitted <- exp(model.matrix(model, dpop) %*% b)

  Nhat   <- sum(fitted)

  phat   <- fitted / Nhat

  dpop   <- dpop[, colnames(dpop) != "Freq"]

  fitted <- cbind(dpop, fitted)

  dpop$phat  <- phat


  if (!is.null(lat)){

    probs <- list()

    for(j in 1:length(lat)){

      probs[[paste0("P(", colnames(dpop[j]), ")")]] <- round(xtabs(paste("phat ~", colnames(dpop[j])), data = dpop), 4)

      q     <- NULL

      names <- NULL

      for (k in (length(lat) + 1):(ncol(dpop) - 1)){

        form  <- paste("phat ~ ", colnames(dpop)[k], "+", colnames(dpop)[j])

        q     <- round(rbind(q, xtabs(form, dpop)/matrix(colSums(xtabs(form, dpop)),
                                                         nlevels(dpop[, k]),
                                                         nlevels(dpop[, j]),
                                                         byrow = T)), 4)

        names <- c(names, rep(colnames(dpop)[k], nlevels(dpop[, k])))
      }

      rownames(q) <- paste(names, "=",  rownames(q))

      colnames(q) <- paste(colnames(dpop)[j], "=", 1:nlevels(dpop[, j]))

      probs[[paste0("P(.|", colnames(dpop)[j],")")]] <- q

    }

  } else if (is.null(lat)) {

    probs <- list()

    q     <- NULL

    names <- NULL

    for (k in 1:(ncol(dpop) - 1)) {

      form  <- paste("phat ~ ", colnames(dpop)[k])

      q     <- round(c(q, xtabs(form, dpop)), 4)

      names <- c(names, rep(colnames(dpop)[k], nlevels(dpop[, k])))
    }


    probs[["Probabilities"]] <- matrix(q, ncol = 1, dimnames = list(paste(names, "=",  names(q)), "est"))


  }



  print(matrix(as.vector(LL), 1, 1, dimnames = list(NULL, "loglikelihood")), digits = 12)
  cat("\n")
  print(matrix(c(n, Nhat, round(Nhat - n)), 1, 3, dimnames = list(NULL, c("nobs", "Nhat", "n0"))))
  cat("\n")
  print(probs)

  cat("\n")
  print(as.matrix(round(coefs, 4), 3))


  invisible(list(
    loglike = LL,
    coefs   = coefs,
    probs   = probs,
    fitted  = fitted,
    hist    = format(hist, scientific = F),
    model   = model,
    nclass  = nclass,
    dobsmis = dobsmis,
    dimp    = dimp,
    D       = D,
    fxc     = fxc))


}


################
## bootstrap ####
################

#' Perform a semi-parametric bootstrap
#'
#' @description \code{boot_mse} Draws random multinomial bootstrap samples of
#' size N from the observed frequencies, including those of cells with missing values,
#' and the fitted frequencies of the cells with structural zeros, and fits the model
#' to the bootstrap samples excluding the cells with structural zeros.
#' @param object an object created with the function \code{\link{mse}}.
#' @param B, the number of bootstrap replications. Defaults to 2000.
#' @param seed random seed for reproducibility. Defaults to 1.
#' @return A matrix with the observed (and latent) variables
#' and the fitted frequencies of the \code{B} bootstrap samples in the columns.
#' @examples
#' # Bootstrap the object:
#' AB <- mse(Freq ~ A*b + B*a + a*b, dframe = NewZealand[, c(1:2, 5:6, 9)])
#' \dontrun{
#' boot_AB <- boot_mse(AB)
#' }
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster detectCores clusterSetRNGStream stopCluster
#' @importFrom iterators icount
#' @importFrom foreach foreach %dopar%
#' @export



boot_mse <- function(object, B = 2000, seed = 1){

  d0    <- object$dobsmis

  dx    <- object$fitted

  n     <- sum(d0$Freq)

  N     <- sum(dx$fitted)

  D     <- object$D

  n0    <- N - n

  ds    <- rbind(d0, c(rep(NA, ncol(d0)  - 1), n0))

  fx    <- object$fxc

  dimp  <- object$dimp

  betas <- object$coefs[-1, 1]

  model <- object$model


  cl <- makeCluster(detectCores()/2)

  registerDoParallel(cl)

  RNGkind("L'Ecuyer-CMRG")

  clusterSetRNGStream(cl = cl, iseed = seed)

  out <- foreach(icount(B), .combine = cbind) %dopar% {

    boot_fit(d0 = d0, dx = dx, n = n, N = N, D = D, n0 = n0, ds = ds, fx = fx, dimp = dimp, betas = betas, model = model)

  }

  stopCluster(cl)

  return(cbind(dx, out))

}






