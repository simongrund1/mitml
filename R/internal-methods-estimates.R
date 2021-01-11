# ***
# global functions
#

# * parameter estimates and variance-covariance matrix

.extractParameters <- function(model, diagonal = FALSE, include.extra.pars = FALSE){

  # number of imputations
  m <- length(model)

  # extract parameter estimates and variance-covariance matrices
  Qhat <- lapply(model, .getCoef, include.extra.pars = include.extra.pars)
  Uhat <- lapply(model, .getVcov, include.extra.pars = include.extra.pars)
  p <- length(Qhat[[1]])
  nms <- names(Qhat[[1]])

  # ensure proper dimensions
  stopifnot(all(p == dim(Uhat[[1]])))
  Qhat <- matrix(unlist(Qhat), nrow = p, ncol = m)
  Uhat <- array(unlist(Uhat), dim = c(p, p, m))

  # extract diagonal
  if(diagonal){
    Uhat <- apply(Uhat, 3, diag)
    if(is.null(dim(Uhat))) dim(Uhat) <- dim(Qhat)
  }

  out <- list(Qhat = Qhat, Uhat = Uhat, nms = nms)
  return(out)

}

# * misc. parameter estimates (e.g., variance components)

.extractMiscParameters <- function(model){

  # number of imputations
  m <- length(model)

  # extract parameter estimates and variance-covariance matrices
  Qhat <- lapply(model, .getMisc)
  p <- length(Qhat[[1]])
  nms <- names(Qhat[[1]])

  # ensure proper dimensions
  if(is.null(Qhat[[1]])){
    Qhat <- NULL
  }else{
    Qhat <- matrix(unlist(Qhat), nrow = p, ncol = m)
  }

  out <- list(Qhat = Qhat, nms = nms)
  return(out)

}

# ***
# generic functions
#

.getCoef <- function(object, ...) UseMethod(".getCoef", object)
.getVcov <- function(object, ...) UseMethod(".getVcov", object)
.getMisc <- function(object, ...) UseMethod(".getMisc", object)

# ***
# default methods
#

.getCoef.default <- function(object, ...) return(coef(object))
.getVcov.default <- function(object, ...) return(as.matrix(vcov(object)))
.getMisc.default <- function(object) return(NULL)

# ***
# class-specific methods
#

# * stats::lm

.getMisc.lm <- function(object){

  # residual variance
  res <- resid(object)
  rv <- sum(res^2) / df.residual(object)
  names(rv) <- "Residual~~Residual"

  return(rv)

}

# * stats::glm

.getMisc.glm <- function(object){

  fam <- tolower(object$family$family)
  if(fam == "gaussian") .getMisc.lm(object)

  return(NULL)

}

# * lme4::(g)lmer

.getCoef.merMod <- function(object, ...) return(lme4::fixef(object))

.getMisc.merMod <- function(object){

  # check if model uses scale
  useSc <- getME(object, "devcomp")$dims["useSc"] == 1

  # variance components by cluster variable
  vc <- lme4::VarCorr(object)
  clus <- names(vc)

  # loop over cluster variables
  out.list <- list()
  for(cc in clus){

    vc.cc <- vc[[cc]]
    if(is.null(dim(vc.cc))) dim(vc.cc) <- c(1, 1)
    nms <- sub("^[(]Intercept[)]$", "Intercept", rownames(vc.cc))

    vc.out <- diag(vc.cc)
    names(vc.out) <- paste0(nms, "~~", nms, "|", cc)

    vc.ind <- which(upper.tri(vc.cc), arr.ind = TRUE)
    for(ii in seq_len(nrow(vc.ind))){
      vc.ii <- vc.cc[vc.ind[ii, , drop=FALSE]]
      names(vc.ii) <- paste0(nms[vc.ind[ii,1]], "~~", nms[vc.ind[ii,2]], "|", cc)
      vc.out <- c(vc.out, vc.ii)
    }

    out.list[[cc]] <- vc.out

  }

  # residual variance (if model uses scale)
  if(useSc){

    rv <- attr(vc, "sc")^2
    names(rv) <- "Residual~~Residual"

    out.list[["Residual"]] <- rv

  }

  # get additional parameters (ICC; only for single clustering)
  if(useSc && length(clus) == 1){

    hasIntercept <- "(Intercept)" %in% colnames(vc[[clus]])
    if(hasIntercept){
      iv <- vc[[clus]]["(Intercept)", "(Intercept)"]
      icc <- iv / (iv + rv)
      names(icc) <- paste("ICC|", clus, sep = "")
    }

    out.list[["ICC"]] <- icc

  }

  out <- do.call(c, unname(out.list))
  return(out)

}

# * nlme::lme
 
.getCoef.lme <- function(object, ...) return(nlme::fixef(object))

.getMisc.lme <- function(object){

  # check if model uses fixed sigma (no scale)
  fixedSigma <- attr(object$modelStruct, "fixedSigma")

  # variance components by cluster variable
  vc.list <- .listVC_lme(object)
  out.list <- list()
  cl <- names(vc.list)
  for(cc in names(vc.list)){

    vc.cc <- vc.list[[cc]]
    nms <- sub("^[(]Intercept[)]$", "Intercept", attr(vc.cc, "nms"))

    vc.out <- diag(vc.cc)
    names(vc.out) <- paste0(nms, "~~", nms, "|", cc)

    vc.ind <- which(upper.tri(vc.cc), arr.ind = TRUE)
    for(ii in seq_len(nrow(vc.ind))){
      vc.ii <- vc.cc[vc.ind[ii, , drop=FALSE]]
      names(vc.ii) <- paste0(nms[vc.ind[ii,1]], "~~", nms[vc.ind[ii,2]], "|", cc)
      vc.out <- c(vc.out, vc.ii)
    }

    out.list[[cc]] <- vc.out

  }

  # residual variance (if model does not use fixed sigma)
  if(!fixedSigma){

    rv <- object$sigma^2
    names(rv) <- "Residual~~Residual"

    out.list[["Residual"]] <- rv

  }

  # get additional parameters (ICC; only for single clustering)
  if(!fixedSigma && length(cl) == 1){

    vc <- vc.list[[cl]]
    rownames(vc) <- colnames(vc) <- attr(vc.list[[cl]], "nms")

    hasIntercept <- "(Intercept)" %in% rownames(vc)
    if(hasIntercept){
      iv <- vc["(Intercept)", "(Intercept)"]
      icc <- iv / (iv + rv)
      names(icc) <- paste("ICC|", cl, sep = "")
    }

    out.list[["ICC"]] <- icc

  }

  out <- do.call(c, unname(out.list))
  return(out)

}

.listVC_lme <- function(object){

  # read random effects structure
  re <- rev(object$modelStruct$reStruct) # see nlme:::VarCorr.lme
  vc <- lapply(re, nlme::VarCorr, rdig = 10^6, sigma = object$sigma)
  cl <- names(vc)
 
  # loop over cluster variables
  vc.list <- list()
  for(cc in cl){

    vc.cc <- vc[[cc]]
    if(is.null(dim(vc.cc))) dim(vc.cc) <- c(1, 1)

    # standard deviation of random effects
    vc.sd <- vc.cc[,"StdDev"]

    # correlation and covariance matrix of random effects
    if(length(vc.sd) == 1){
      vc.cov <- as.matrix(vc.sd^2)
      attr(vc.cov, "nms") <- rownames(vc.cc)[1]
    }else{
      vc.cor <- cbind(attr(vc.cc, "corr"), "")
      vc.cor[upper.tri(vc.cor, diag = TRUE)] <- ""
      storage.mode(vc.cor) <- "numeric"
      diag(vc.cor) <- 1
      vc.cor[upper.tri(vc.cor)] <- t(vc.cor)[upper.tri(t(vc.cor))]

      # calculate covariance matrix
      vc.sd <- diag(vc.cc[,"StdDev"])
      vc.cov <- vc.sd %*% vc.cor %*% vc.sd
      attr(vc.cov, "nms") <- rownames(vc.cc)
    }

    vc.list[[cc]] <- vc.cov

  }

  return(vc.list)

}

# * geepack::geeglm

.getMisc.geeglm <- function(object){

  fixedScale <- length(object$geese$gamma) == 0
  fixedCor <- length(object$geese$alpha) == 0

  # scale parameter (gamma)
  out.list <- list()
  if(!fixedScale){

    gam <- object$geese$gamma
    nms <- gsub("^[(]Intercept[)]$", "Intercept", names(gam))
    names(gam) <- paste0("Scale:", nms)

    out.list[["gamma"]] <- gam

  }

  # correlation parameters (alpha)
  if(!fixedCor){

    alpha <- object$geese$alpha
    names(alpha) <- paste0("Correlation:", names(alpha))

    out.list[["alpha"]] <- alpha

  }

  out <- do.call(c, unname(out.list))
  return(out)

}

# * lavaan::lavaan

.getCoef.lavaan <- function(object, include.extra.pars = FALSE, ...){

  # extract parameter estimates
  pt <- lavaan::parTable(object)
  isFree <- pt[["free"]] > 0 & !duplicated(pt[["free"]]) # see lavaan:::lav_object_inspect_coef
  isDefined <- pt[["op"]] == ":="
  isCoef <- if(include.extra.pars) isFree | isDefined else isFree
  
  hasLevels <- length(lavaan::lavInspect(object, "cluster")) > 0
  hasGroups <- length(lavaan::lavInspect(object, "group")) > 0

  # parameter names
  # NOTE: replaces names in coef() with names that are independent of labels
  # assigned to parameters (can be inconsistent across models)
  nms <- pt[, c("lhs", "op", "rhs")]
  if(hasLevels) nms[["level"]] <- paste0(".l", pt[, "level"])
  if(hasGroups) nms[["group"]] <- paste0(".g", pt[, "group"])
  nms <- do.call(mapply, c(as.list(nms), list(FUN = paste0)))

  out <- pt[isCoef, "est"]
  names(out) <- nms[isCoef]
  return(out)

}

.getVcov.lavaan <- function(object, include.extra.pars = FALSE, ...){

  pt <- lavaan::parTable(object)
  hasDefined <- any(pt[["op"]] == ":=")

  if(hasDefined && include.extra.pars){
    out <- lavaan::lavInspect(object, "vcov.def.joint")
  }else{
    out <- lavaan::lavInspect(object, "vcov")
  }
  rownames(out) <- colnames(out) <- NULL

  return(out)

}

.getMisc.lavaan <- function(object){

  # extract (nonfree) parameter estimates
  pt <- lavaan::parTable(object)
  isFree <- pt[["free"]] > 0 & !duplicated(pt[["free"]]) # see lavaan:::lav_object_inspect_coef
  isDefined <- pt[["op"]] == ":="
  isCoef <- isFree | isDefined
  
  hasLevels <- length(lavaan::lavInspect(object, "cluster")) > 0
  hasGroups <- length(lavaan::lavInspect(object, "group")) > 0

  # parameter names
  nms <- pt[, c("lhs", "op", "rhs")]
  if(hasLevels) nms[["level"]] <- paste0(".l", pt[, "level"])
  if(hasGroups) nms[["group"]] <- paste0(".g", pt[, "group"])
  nms <- do.call(mapply, c(as.list(nms), list(FUN = paste0)))

  out <- pt[!isCoef, "est"]
  names(out) <- nms[!isCoef]
  return(out)

}


