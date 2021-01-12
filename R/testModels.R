testModels <- function(model, null.model, method = c("D1", "D2", "D3", "D4"), use = c("wald", "likelihood"), ariv = c("default", "positive", "robust"), df.com = NULL, data = NULL){
# model comparison and hypothesis tests for k-dimensional estimands

  # ***
  # check input
  #

  # check model specification
  m <- length(model)
  if(!(is.list(model) && is.list(null.model))) stop("The 'model' and 'null.model' arguments must be lists of fitted statistical models.")
  if(length(null.model) != m) stop("The 'model' and 'null.model' arguments must be lists with the same length.")

  # match methods
  method <- match.arg(method)
  use <- match.arg(use)
  ariv <- match.arg(ariv)

  # check for incompatible arguments
  if(!is.null(df.com) && method != "D1") warning("Complete-data degrees of freedom are not available for use with '", method, "' and were ignored.")
  if(use == "likelihood" && method != "D2") warning("The 'likelihood' option is not available with method '", method ,"' and was ignored.")
  if(!is.null(data) && method != "D4") warning("The 'data' argument is not used with method '", method ,"' and was ignored.")
  if(ariv == "positive" && method == "D1") warning("The 'positive' option is not available with method 'D1' and was ignored.")
  if(ariv == "robust" && method != "D4") warning("The 'robust' option is not available with method '", method ,"' and was ignored.")

  # check model classes
  cls <- class(model[[1]])
  cls.null <- class(null.model[[1]])

  if(cls[1] != cls.null[1]) warning("The 'model' and 'null.model' arguments appear to include objects of different classes. Results may not be trustworthy.")

  .checkNamespace(union(cls, cls.null))

  # check for REML and refit (if needed)
  reml.model <- sapply(model, .checkREML)
  reml.null.model <- sapply(null.model, .checkREML)
  reml <- any(reml.model, reml.null.model)

  if(reml){
    model[reml.model] <- lapply(model[reml.model], .updateML)
    null.model[reml.null.model] <- lapply(null.model[reml.null.model], .updateML)
  }

  # ***
  # D1
  #

  if(method == "D1"){

    # FIXME: better way to handle this?
    if(inherits(model[[1]], "lavaan")) stop("The 'D1' method is currently not supported for objects of class 'lavaan'. Please see '?testModels' for a list of supported model types.")

    est <- .extractParameters(model, diagonal = FALSE)
    est.null <- .extractParameters(null.model, diagonal = FALSE)

    par.diff <- est$nms[!(est$nms %in% est.null$nms)]
    par.ind <- match(par.diff, est$nms)
    if(length(par.diff) == 0L) stop("The 'model' and 'null.model' appear not to be nested or include the same set of parameters.")

    k <- length(par.diff)
    Qhat <- est$Qhat[par.ind,, drop = FALSE]
    Uhat <- est$Uhat[par.ind, par.ind,, drop = FALSE]

    Qbar <- apply(Qhat, 1, mean)
    Ubar <- apply(Uhat, c(1,2), mean)
    B <- cov(t(Qhat))
    r <- (1+m^(-1))*sum(diag(B%*%solve(Ubar)))/k
    Ttilde <- (1 + r)*Ubar

    # D1 (Li, Raghunathan and Rubin, 1991)
    val <- t(Qbar) %*% solve(Ttilde) %*% Qbar / k

    t <- k*(m-1)
    if(!is.null(df.com)){
      a <- r*t/(t-2)
      vstar <- ( (df.com+1) / (df.com+3) ) * df.com
      v <- 4 + ( (vstar-4*(1+a))^(-1) + (t-4)^(-1) * ((a^2*(vstar-2*(1+a))) /
           ((1+a)^2*(vstar-4*(1+a)))) )^(-1)
    }else{
      if (t>4){
        v <- 4 + (t-4) * (1 + (1 - 2*t^(-1)) * (r^(-1)))^2
      }else{
        v <- t * (1 + k^(-1)) * ((1 + r^(-1))^2) / 2
      }
    }

  }

  # ***
  # D2
  #

  if(method == "D2"){

    if(use == "wald"){

      # FIXME: better way to handle this?
      if(inherits(model[[1]], "lavaan")) stop("The 'D2' method currently only supports likelihood-based comparisons for objects of class 'lavaan'. Please see '?testModels' for a list of supported model types.")

      reml <- FALSE

      # extract parameter estimates
      est <- .extractParameters(model, diagonal = FALSE)
      est.null <- .extractParameters(null.model, diagonal = FALSE)

      par.diff <- est$nms[!(est$nms %in% est.null$nms)]
      par.ind <- match(par.diff, est$nms)
      if(length(par.diff) == 0L) stop("The 'model' and 'null.model' appear not to be nested or include the same set of parameters.")

      # Wald tests
      k <- length(par.diff)
      Qhat <- est$Qhat[par.ind,, drop = FALSE]
      Uhat <- est$Uhat[par.ind, par.ind,, drop = FALSE]

      dW <- sapply(seq_len(m), function(z) t(Qhat[,z]) %*% solve(Uhat[,,z]) %*% Qhat[,z])

    }

    if(use == "likelihood"){

      # extract logLik
      ll <- .evaluateLogLik(model)
      ll.null <- .evaluateLogLik(null.model)
      ll.diff <- ll.null$LL - ll$LL

      if(is.null(ll$df) || is.null(ll.null$df)) stop("Degrees of freedom for the model comparison could not be detected.")
      k <- ll$df - ll.null$df

      # account for numerical imprecision
      isEqual <- mapply(function(x, y) isTRUE(all.equal(x, y)), x = ll$LL, y = ll.null$LL)
      ll.diff[isEqual] <- 0L

      # LR tests
      dW <- -2 * (ll.diff)

    }

    # D2 (Li, Meng et al., 1991)
    dWbar <- mean(dW)
    r <- (1+m^(-1)) * var(sqrt(dW))
    if(ariv == "positive") r <- max(0, r)
    val <- (dWbar/k - (m+1)/(m-1) * r) / (1+r)

    v <- k^(-3/m) * (m-1) * (1+r^(-1))^2

  }

  # ***
  # D3
  #

  if(method == "D3"){

    # evaluate log-likelihood at estimated and pooled values of model parameters
    ll <- .evaluateUserLogLik(model)
    ll.null <- .evaluateUserLogLik(null.model)

    k <- ll$df - ll.null$df

    # D3 (Meng & Rubin, 1992)
    dL.bar <- mean(-2 * (ll.null$LL - ll$LL))
    dL.tilde <- mean(-2 * (ll.null$LL.pooled - ll$LL.pooled))

    r <- (m+1) * (k*(m-1))^(-1) * (dL.bar - dL.tilde)
    if(ariv == "positive") r <- max(0, r)
    val <- dL.tilde / (k*(1+r))

    t <- k*(m-1)
    if( t > 4 ){
      v <- 4 + (t-4) * (1 + (1-2*t^(-1)) * r^(-1))^2
    }else{
      v <- t * (1+k^(-1)) * (1+r^(-1))^2 / 2
    }

    use <- "likelihood"

  }

  # ***
  # D4
  #

  if(method == "D4"){

    # evaluate log-likelihood at estimated and pooled values of model parameters
    ll <- .evaluateStackedLogLik(model, datalist = data)
    ll.null <- .evaluateStackedLogLik(null.model, datalist = data)
    ll.diff <- ll.null$LL - ll$LL

    ll.stacked.diff <- ll.null$LL.stacked - ll$LL.stacked

    k <- ll$df - ll.null$df
    h <- ll$df

    # account for numerical imprecision
    if(isTRUE(all.equal(ll.stacked.diff[1], 0))) ll.stacked.diff <- 0L
    isEqual <- mapply(function(x, y) isTRUE(all.equal(x, y)), x = ll$LL, y = ll.null$LL)
    ll.diff[isEqual] <- 0L

    # D4 (Chan & Meng, 2019)
    dbar <- mean(-2 * ll.diff)
    dhat <- -2 * ll.stacked.diff

    if(ariv == "robust"){

      deltabar <- 2 * mean(ll$LL)
      deltahat <- 2 * ll$LL.stacked
      r <- (m+1) / (h*(m-1)) * (deltabar - deltahat)
      v <- (h*(m-1)) * (1 + 1/r)^2

    }else{

      r <- (m+1) / (k*(m-1)) * (dbar - dhat)
      if(ariv == "positive") r <- max(0, r)
      v <- (k*(m-1)) * (1 + r^(-1))^2

    }

    val <- dhat / (k*(1+r))
    use <- "likelihood"

  }

  # create output
  pval <- pf(val, k, v, lower.tail = FALSE)
  out <- matrix(c(val, k, v, pval, r), ncol = 5)
  colnames(out) <- c("F.value", "df1", "df2", "P(>F)", "RIV") # new label for p-value, SiG 2017-02-09

  out <- list(
    call = match.call(),
    test = out,
    m = m,
    method = method,
    adj.df = !is.null(df.com),
    df.com = df.com,
    use = use,
    ariv = ariv,
    data = !is.null(data),
    reml = reml
  )

  class(out) <- "mitml.testModels"
  return(out)

}

