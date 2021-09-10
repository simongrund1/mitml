# ***
# global functions
#

# * log-likelihood

.evaluateLogLik <- function(model){

  cls <- class(model[[1]])

  # evaluate log-likelihood
  ll <- lapply(model, .getLL)
  if(any(sapply(ll, is.null))) stop("Could not evaluate likelihood for objects of class '", paste0(cls, collapse = "|"), "'. Please see '?testModels' for a list of supported model types.")
  df <- attr(ll[[1]], "df")

  # ensure proper dimensions
  ll <- unlist(ll)

  out <- list(LL = ll, df = df)
  return(out)

}

# * log-likelihood evaluated at user-defined values

.evaluateUserLogLik <- function(model){

  m <- length(model)
  cls <- class(model[[1]])

  # extract arguments and function for likelihood evaluation
  ll.args <- lapply(model, .getArgsLL)
  if(any(sapply(ll.args, is.null))) stop("Could not evaluate likelihood for objects of class '", paste0(cls, collapse = "|"), "'. Please see '?testModels' for a list of supported model types.")
  narg <- length(ll.args[[1]][["parameters"]])
  nms <- names(ll.args[[1]][["parameters"]])

  # evaluate log-likelihood at imputation-specific parameter values
  ll <- lapply(lapply(ll.args, c, force.update = FALSE), do.call, what = .getUserLL)
  df <- attr(ll[[1]], "df")

  ll <- unlist(ll)

  # pool parameter estimates
  psi.bar <- vector("list", narg)
  names(psi.bar) <- nms

  for(i in seq_along(psi.bar)){

    psi <- lapply(ll.args, function(x, .i) x$parameters[[.i]], .i = i)
    isMatrix <- is.matrix(psi[[1]])

    if(isMatrix){
      q <- nrow(psi[[1]])
      pp <- array(unlist(psi), dim = c(q, q, m))
      pp <- apply(pp, c(1, 2), mean)
      rownames(pp) <- colnames(pp) <- names(psi[[1]])
      psi.bar[[i]] <- pp
    }else{
      q <- length(psi[[1]])
      pp <- matrix(unlist(psi), nrow = q, ncol = m)
      pp <- apply(pp, 1, mean)
      names(pp) <- names(psi[[1]])
      psi.bar[[i]] <- pp
    }

  }

  # evaluate log-likelihood at pooled parameter estimates
  for(i in seq_len(m)) ll.args[[i]]$parameters <- psi.bar
  ll.pooled <- sapply(ll.args, do.call, what = .getUserLL)

  out <- list(LL = ll, LL.pooled = ll.pooled, df = df)
  return(out)

}

# * log-likelihood evaluated with stacked data sets

.evaluateStackedLogLik <- function(model, datalist = NULL){

  m <- length(model)
  cls <- class(model[[1]])

  # evaluate log-likelihood
  ll <- lapply(model, .getLL)
  if(any(sapply(ll, is.null))) stop("Could not evaluate likelihood for objects of class '", paste0(cls, collapse = "|"), "'. Please see '?testModels' for a list of supported model types.")
  df <- attr(ll[[1]], "df")

  # ensure proper dimensions
  ll <- unlist(ll)

  # extract data for stacking
  nullData <- is.null(datalist)
  if(nullData) datalist <- lapply(model, .getDataLL)

  # check data
  if(!is.list(datalist) || length(datalist) != m || !all(sapply(datalist, is.data.frame))){
    if(nullData){
      stop("Could not extract data from fitted model objects. Please specify 'data' and see '?testModels' for details.")
    }else{
      stop("The 'data' argument must be a list of imputed data sets that correspond to the fitted model objects. Please see '?testModels' for details.")
    }
  }

  # evaluate log-likelihood with stacked data
  model.stacked <- .updateStackedLL(model[[1]], datalist = datalist)
  ll.stacked <- .getLL(model.stacked) / m

  out <- list(LL = ll, LL.stacked = ll.stacked, df = df)
  return(out)

}

# ***
# generic functions
#

.getLL <- function(object, ...) UseMethod(".getLL", object)
.getArgsLL <- function(object, ...) UseMethod(".getArgsLL", object)
.getUserLL <- function(object, ...) UseMethod(".getUserLL", object)
.getDataLL <- function(object, ...) UseMethod(".getDataLL", object)
.updateStackedLL <- function(object, ...) UseMethod(".updateStackedLL", object)

# ***
# default methods
#

.getLL.default <- function(object) return(logLik(object))
.getArgsLL.default <- function(object) return(NULL)
.getUserLL.default <- function(object, ...) return(NULL)
.getDataLL.default <- function(object) return(model.frame(object))
.updateStackedLL.default <- function(object, datalist) return(update(object, data = do.call(rbind, datalist)))

# ***
# class-specific methods
#

# * stats::lm

.getArgsLL.lm <- function(object){

  # extract arguments to evaluate LL
  n <- nrow(object$model)
  beta <- coef(object)
  sigma2 <- sum(resid(object)^2) / n

  out <- list(object = object, parameters = list(beta = beta, sigma2 = sigma2))
  return(out)

}

.getUserLL.lm <- function(object, parameters, ...){

  n <- nrow(object$model)
  df <- object$rank + 1
  trm <- attributes(object$terms)

  # model matrices
  y <- eval(trm$variables, envir = object$model)[[trm$response]]
  X <- model.matrix(object)

  # parameters
  beta <- parameters[["beta"]]
  sigma2 <- parameters[["sigma2"]]

  ll <- .logLik_lm(y = y, X = X, beta = beta, sigma2 = sigma2)
  attr(ll, "df") <- df

  return(ll)

}

.logLik_lm <- function(y, X, beta, sigma2){
  n <- length(y)
  - (n/2) * log(2*pi*sigma2) - (1/(2*sigma2)) * sum((y - X %*% beta)^2)
}

# * stats::glm

.getArgsLL.glm <- function(object) return(NULL)

# * geepack::geeglm

.getLL.geeglm <- function(object) return(NULL)
.getArgsLL.geeglm <- function(object) return(NULL)

# * lme4::lmer (for only LMMs)

.getArgsLL.lmerMod <- function(object){

  beta <- lme4::getME(object, "fixef")
  theta <- lme4::getME(object, "theta")
  sig <- sigma(object)

  # split theta by clustering variables
  cl <- lme4::getME(object, "cnms")
  ncl <- length(cl)
  nvc <- lengths(cl)
  theta.cl <- split(theta, rep.int(seq_along(cl), (nvc * (nvc + 1))/2))

  # transform theta from scaled Cholesky factors into variance-covariance matrices (for pooling)
  Tau <- vector("list", ncl)
  names(Tau) <- paste0("Tau", seq_len(ncl))

  for(i in seq_len(ncl)){
    q <- sqrt(2*length(theta.cl[[i]]) + 0.25) - 0.5
    m <- matrix(0, nrow = q, ncol = q)
    m[lower.tri(m, diag = TRUE)] <- theta.cl[[i]] * sig
    Tau[[i]] <- m %*% t(m)
  }

  out <- list(object = object, parameters = c(list(beta = beta), Tau, list(sigma2 = sig^2)))
  return(out)

}

.getUserLL.lmerMod <- function(object, parameters, force.update = TRUE, ...){

  if(any(abs(lme4::getME(object, "offset") - 0) > .Machine$double.eps)) stop("The 'D3' method cannot be used for 'lmerMod' objects fitted with an offset.")

  cl <- lme4::getME(object, "cnms")
  ncl <- length(cl)

  # evaluate standard logLik
  ll0 <- logLik(object)
  df <- attr(ll0, "df")

  if(force.update){

    # get fixed-effects linear predictor
    X <- lme4::getME(object, "X")
    beta <- parameters$beta
    linpred <- X %*% beta

    # update formula
    fml <- as.formula(sub("~", "~ 0 +", deparse(formula(object, random.only = TRUE))))

    # update model with fixed contribution of fixed effects
    newobj <- .localUpdate(object, formula = fml, data = model.frame(object), offset = linpred)

    # get variance components
    Tau <- parameters[grep("^Tau", names(parameters))]
    sig <- sqrt(parameters$sigma2)

    # transform variance-covariance matrices into correlations and SDs (for devfun)
    theta.cl <- vector("list", ncl)
    for(i in seq_len(ncl)){
      v <- Tau[[i]]
      r <- lme4::cov2sdcor(v)
      theta.cl[[i]] <- r[lower.tri(r, diag = TRUE)]
    }
    theta <- c(do.call(c, theta.cl), sig)

    # evaluate (profiled) deviance with fixed theta
    dev.fun <- lme4::devfun2(newobj)
    ll <- -dev.fun(pars = theta) / 2
    attr(ll, "df") <- df

  }else{

    ll <- ll0[1]
    attr(ll, "df") <- df

  }

  return(ll)

}

# * lme4::(g)lmer (for both LMMs and GLMMs)

.updateStackedLL.merMod <- function(object, datalist){

  # create imputation-specific levels for clustering variables
  cl <- lme4::getME(object, "cnms")
  for(ii in seq_along(datalist)){
    for(cc in names(cl)){
      datalist[[ii]][,cc] <- paste0("imp", ii, "_", datalist[[ii]][,cc])
    }
  }

  # stack data
  stackdat <- do.call(rbind, datalist)
  for(cc in names(cl)) stackdat[,cc] <- as.integer(as.factor(stackdat[,cc]))

  # update model with stacked data
  # NOTE: update.merMod will find global objects of the same name before local ones (very bad),
  # so we need to update in a separate environment
  env <- new.env()
  assign("stackdat", value = stackdat, envir = env)
  newobj <- .localUpdate(object, envir = env, data = stackdat)

  return(newobj)

}

# * nlme::lme

.getArgsLL.lme <- function(object){

  beta <- nlme::fixef(object)
  Tau <- .listVC_lme(object)
  names(Tau) <- paste0("Tau", seq_along(Tau))
  sigma2 <- sigma(object)^2

  out <- list(object = object, parameters = c(list(beta = beta), Tau, list(sigma2 = sigma2)))
  return(out)

}

.getUserLL.lme <- function(object, parameters, ...){

  ncl <- object$dims$Q # see nlme:::print.summary.lme

  if(ncl > 1) stop("The 'D3' method is only supported for models of class 'lme' with a single cluster variable. Please see '?testModels' for a list of supported model types.")

  # evaluate standard logLik
  p <- object$dims$ncol[[object$dims$Q + 1]] # see nlme:::logLik.lme
  fixedSigma <- attr(object[["modelStruct"]], "fixedSigma")
  df <- p + length(coef(object[["modelStruct"]])) + as.integer(!fixedSigma)

  # response and cluster variables
  y <- nlme::getResponse(object)
  clus <- nlme::getGroups(object)

  # fixed and random effects formulas
  fe.fml <- eval(eval(object$call$fixed)[-2]) # see nlme:::predict.lme
  re.str <- object$modelStruct$reStruct

  # fixed effects and design matrix
  X <- model.matrix(fe.fml, object$data)
  beta <- parameters$beta

  # random effects variance components and design matrix
  Z <- model.matrix(re.str, object$data)
  Tau <- parameters[[grep("^Tau", names(parameters))]]
  sigma2 <- parameters$sigma2

  # evaluate log-likelihood
  ll <- .logLik_lmm(y = y, X = X, Z = Z, cluster = clus, beta = beta, 
                    Tau = Tau, sigma2 = sigma2)
  attr(ll, "df") <- df

  return(ll)

}

.getDataLL.lme <- function(object){

  out <- nlme::getData(object)
  return(out)

}

.updateStackedLL.lme <- function(object, datalist){

  # add levels to clustering variables
  re <- rev(object$modelStruct$reStruct) # see nlme:::VarCorr.lme
  cl <- names(re)
  for(ii in seq_along(datalist)){
    for(cc in cl){
      datalist[[ii]][,cc] <- paste0("imp", ii, "_", datalist[[ii]][,cc])
    }
  }

  # update model with stacked data
  stackdat <- do.call(rbind, datalist)
  for(cc in names(cl)) stackdat[,cc] <- as.integer(as.factor(stackdat[,cc]))
  newobj <- update(object, data = stackdat)

  return(newobj)

}

.logLik_lmm <- function(y, X, Z, cluster, beta, Tau, sigma2){

  p <- length(beta)
  q <- dim(Tau)[1]
  y <- split(y, cluster)
  X <- split(X, cluster)
  Z <- split(Z, cluster)

  lvls <- unique(cluster)
  L <- numeric(length(lvls))
  for(i in seq_along(lvls)){

    yi <- y[[i]]
    ni <- length(yi)
    Xi <- matrix(X[[i]], nrow = ni, ncol = p)
    Ri <- yi - Xi%*%beta
    Zi <- matrix(Z[[i]], nrow = ni, ncol = q)

    V <- diag(sigma2, ni) + Zi %*% Tau %*% t(Zi)
    Vinv <- chol2inv(chol(V))

    dV <- determinant(V, logarithm = TRUE)
    dV <- dV$modulus * dV$sign
    L[i] <- dV + t(Ri) %*% Vinv %*% (Ri)
  }

  -sum(L)/2

}

# * lavaan::lavaan

.getLL.lavaan <- function(object){

  # FIXME: catch scaled LRT statistics (currently not supported)
  # see lavaan::lavTestLRT
  tests <- unlist(sapply(slot(object, "test"), "[", "test"))
  isScaled <- c("satorra.bentler", "yuan.bentler", "yuan.bentler.mplus", "mean.var.adjusted",
                "scaled.shifted") %in% tests
  if(any(isScaled)){
    return(NULL)
  }

  ll <- lavaan::logLik(object)
  return(ll)

}

.getArgsLL.lavaan <- function(object){

  # FIXME: catch scaled LRT statistics (currently not supported)
  # see lavaan::lavTestLRT
  tests <- unlist(sapply(slot(object, "test"), "[", "test"))
  isScaled <- c("satorra.bentler", "yuan.bentler", "yuan.bentler.mplus", "mean.var.adjusted",
                "scaled.shifted") %in% tests
  if(any(isScaled)){
    return(NULL)
  }

  # get parameter table
  pt <- lavaan::parTable(object)
  isFree <- pt[["free"]] > 0

  out <- list(object = object, parameters = list(free = pt[isFree, "est"]))
  return(out)

}

.getUserLL.lavaan <- function(object, parameters, force.update = TRUE, ...){

  ll0 <- lavaan::logLik(object)
  df <- attr(ll0, "df")

  if(force.update){

    # get parameter table
    pt <- lavaan::parTable(object)
    isFree <- pt[["free"]] > 0
    isConstraint <- pt[["op"]] %in% c(":=", "==", "<", ">")

    # fix free parameters to user-defined values
    pt[isFree, c("est", "se", "start")] <- NA
    pt[isFree, "ustart"] <- parameters$free[pt[isFree, "free"]]
    pt[["free"]] <- 0
    pt[["user"]] <- 1

    # remove defined parameters
    pt <- pt[!isConstraint,]

    # extract data
    data <- .restoreData_lavaan(object)

    # update model with fixed parameters
    newobj <- .localUpdate(object, model = pt, data = data)
    ll <- lavaan::logLik(newobj)[1]

  }else{

    ll <- ll0[1]

  }

  attr(ll, "df") <- df
  return(ll)

}

.getDataLL.lavaan <- function(object){

  out <- .restoreData_lavaan(object)
  return(out)

}

.updateStackedLL.lavaan <- function(object, datalist){

  # create imputation-specific levels for clustering variables
  cl <- lavaan::lavInspect(object, "cluster")
  hasLevels <- length(cl) > 0

  if(hasLevels){

    # add levels to clustering variables
    for(ii in seq_along(datalist)) datalist[[ii]][,cl] <- paste0("imp", ii, "_", datalist[[ii]][,cl])

    # stack data
    stackdat <- do.call(rbind, datalist)
    stackdat[,cl] <- as.integer(as.factor(stackdat[,cl]))

  }else{
    stackdat <- do.call(rbind, datalist)
  }

  # update model with stacked data
  newobj <- .localUpdate(object, data = stackdat)

  return(newobj)

}

.restoreData_lavaan <- function(object){

    # extract data
    data <- lavaan::lavInspect(object, "data")
    grp <- lavaan::lavInspect(object, "group")
    cl <- lavaan::lavInspect(object, "cluster")
    
    # re-add group and cluster indicators
    hasGroups <- length(grp) > 0
    hasLevels <- length(cl) > 0
    data <- if(hasGroups) lapply(data, as.data.frame) else as.data.frame(data)
    if(hasGroups){
      grp.nms <- lavaan::lavInspect(object, "group.label")
      for(ii in seq_along(grp.nms)) data[[ii]][,grp] <- grp.nms[ii]
      data <- do.call(rbind, data)
    }
    if(hasLevels){
      cc <- lavaan::lavInspect(object, "cluster.label")
      if(hasGroups) cc <- do.call(c, cc)
      data[,cl] <- cc
    }

    return(data)

}

