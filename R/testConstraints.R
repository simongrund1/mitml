testConstraints <- function(model, qhat, uhat, constraints, method = c("D1", "D2"), ariv = c("default", "positive"), df.com = NULL){
# test constraints with multiply imputed data

  # ***
  # check input
  #

  if(missing(model) == (missing(qhat) || missing(uhat))){
    stop("Either 'model' or both 'qhat' and 'uhat' must be supplied.")
  }

  # match methods
  method <- match.arg(method)
  ariv <- match.arg(ariv)

  # warnings for ignored arguments
  if(!is.null(df.com) && method == "D2") warning("Complete-data degrees of freedom are not available for use with 'D2', and thus were ignored.")
  if(ariv == "positive" && method == "D1") warning("The 'positive' option is not available with method 'D1' and was ignored.")

  # clean names in constraints
  constraints <- gsub("\\(Intercept\\)", "Intercept", constraints)
  k <- length(constraints)

  # ***
  # process matrix, array or list arguments
  #

  if(!missing(qhat)){

    # check input
    if(missing(uhat)) stop("Either 'model' or both 'qhat' and 'uhat' must be supplied.")
    if(!is.matrix(qhat) && is.array(qhat)) stop("The 'qhat' argument must be either a matrix or a list.")
    if(is.matrix(uhat)) stop("The 'uhat' argument must be either an array or a list.")

    # convert point estimates
    if(is.matrix(qhat)){
      qhat <- lapply(seq_len(ncol(qhat)), function(i, Q) Q[,i], Q = qhat)
    }
 
    # convert variance estimates
    if(is.array(uhat)){
      uhat <- lapply(seq_len(dim(uhat)[3]), function(i, U) as.matrix(U[,,i]), U = uhat)
    }

    # ensure proper format
    m <- length(qhat)
    p <- length(qhat[[1]])
    if(m != length(uhat) || !is.matrix(uhat[[1]]) || p != ncol(uhat[[1]]) || p != nrow(uhat[[1]])) stop("Dimensions of 'qhat' and 'uhat' do not match.")

    Qhat <- matrix(unlist(qhat), ncol = m)
    Uhat <- array(unlist(uhat), dim = c(p, p, m))
    if(any(!is.finite(Uhat))) stop("Missing values in 'uhat' are not allowed.")

    nms <- names(qhat[[1]])

  }

  # ***
  # process fitted models
  #

  if(!missing(model)){

    if(!is.list(model)) stop("The 'model' argument must be a list of fitted statistical models.")

    # get class (check required packages)
    cls <- class(model[[1]])
    .checkNamespace(cls)

    # extract parameter estimates
    est <- .extractParameters(model)

    Qhat <- est$Qhat
    Uhat <- est$Uhat
    nms <- est$nms

    m <- length(model)
    p <- nrow(Qhat)

  }

  # ***
  # delta method
  #

  # prepare parameter names
  if(is.null(nms)) stop("Could not determine parameter names.")
  nms <- gsub("\\(Intercept\\)", "Intercept", nms)
  rownames(Qhat) <- nms
  dimnames(Uhat) <- list(nms, nms, NULL)

  newQhat <- array(NA, dim = c(k, m))
  newUhat <- array(NA, dim = c(k, k, m))

  for(ii in 1:m){

    theta <- Qhat[,ii]
    Sigma <- as.matrix(Uhat[,,ii])

    g <- parse(text = constraints)
    env.g <- new.env()
    for(pp in 1:p) assign(names(theta)[pp], theta[pp], pos = env.g)

    # new parameter estimates
    newtheta <- numeric(k)
    for(kk in seq_len(k)) newtheta[kk] <- eval(g[kk], envir = env.g)

    # derivative, new covariance matrix
    gdash.theta <- matrix(NA, k, p)
    for(kk in seq_len(k)){
      tmp <- numericDeriv(g[[kk]], names(theta), env.g)
      gdash.theta[kk,] <- attr(tmp, "gradient")
    }
    newSigma <- gdash.theta %*% Sigma %*% t(gdash.theta)

    newQhat[,ii] <- newtheta
    newUhat[,,ii] <- newSigma

  }

  # ***
  # pool results
  #

  # calculate pooled estimates and covariance matrix (for output)
  Qbar <- apply(newQhat, 1, mean)
  Ubar <- apply(newUhat, c(1, 2), mean)

  B <- cov(t(newQhat))
  r <- (1+m^(-1)) * sum(diag(B%*%solve(Ubar))) / k

  Ttilde <- (1+r) * Ubar

  # D1 (Li et al., 1991)
  if(method == "D1"){

    D <- .D1(Qhat = newQhat, Uhat = newUhat, df.com = df.com)

    r <- D$r
    val <- D$F
    v <- D$v

  }

  # D2 (Li, Meng et al., 1991)
  if(method == "D2"){

    dW <- numeric(m)
    for(ii in seq_len(m)) dW[ii] <- t(newQhat[,ii]) %*% solve(newUhat[,,ii]) %*% newQhat[,ii]

    D <- .D2(d = dW, k = k)

    r <- D$r
    if(ariv == "positive") r <- max(0, r)
    val <- D$F
    v <- D$v

  }

  # create output
  pval <- pf(val, k, v, lower.tail = FALSE)
  out <- matrix(c(val, k, v, pval, r), ncol = 5)
  colnames(out) <- c("F.value", "df1", "df2", "P(>F)", "RIV") # new label for p-value, SiG 2017-02-09

  out <- list(
    call = match.call(),
    constraints = constraints,
    test = out,
    Qbar = Qbar,
    T = Ttilde,
    m = m,
    adj.df = !is.null(df.com),
    df.com = df.com,
    method = method
  )

  class(out) <- "mitml.testConstraints"
  return(out)

}

