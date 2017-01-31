# ***
# Functions to extract variance components from supported classes
# of statistical models
#

# *** lmer method
.getVC.lmer <- function(model){

  if(!requireNamespace("lme4", quietly=TRUE)) stop("The 'lme4' package must be installed in order to use this function.")
  m <- length(model)
  vlist <- addp <- NULL

  # variance components
  vlist <- list()
  vc <- lapply(model,lme4::VarCorr)
  clus <- names(vc[[1]])
  for(vv in clus){
    q <- dim(vc[[1]][[vv]])[1]
    v.cl <- vapply(vc, function(z) z[[vv]], FUN.VALUE=matrix(0,q,q))
    if(is.null(dim(v.cl))) dim(v.cl) <- c(1,1,m)
    dimnames(v.cl)[1:2] <- lapply(dimnames(vc[[1]][[vv]]), function(z) sub("^[(]Intercept[)]$","Intercept",z))
    vlist[[paste("|",vv,sep="")]] <- v.cl
  }
  # residual variance (if model uses scale)
  usesc <- attr(vc[[1]], "useSc")
  if(usesc){
    rv <- sapply(vc, function(z) attr(z,"sc")^2)
    dim(rv) <- c(1,1,m)
    dimnames(rv) <- list("Residual","Residual",NULL)
    vlist <- c(vlist, list(rv))
  }

  # additional parameters
  # 1. ICC (only single clustering)
  if(usesc & length(clus)==1){
    if("(Intercept)"%in%colnames(vc[[1]][[clus]])){
      iv <- sapply(vc, function(z) z[[clus]]["(Intercept)","(Intercept)"])
      icc <- iv / (iv + rv[1,1,])
      addp <- c(addp, mean(icc))
      names(addp) <- paste("ICC|",clus,sep="")
    }
  }

  list(vlist=vlist,addp=addp)
}

# *** nlme method
.getVC.nlme <- function(model){

  if(!requireNamespace("nlme", quietly=TRUE)) stop("The 'nlme' package must be installed in order to use this function.")
  m <- length(model)
  vlist <- addp <- NULL

  # variance components (single clustering, limited for multiple clustering)
  cls <- class(model[[1]])[1]
  clus <- names(model[[1]]$coefficients$random)
  vlist <- list()
  # single cluster variable
  if(cls=="lme" & length(clus)==1){
    vc <- lapply(model,nlme::getVarCov)
    clus <- attr(vc[[1]],"group.levels")
    q <- dim(vc[[1]])[1]
    v.cl <- vapply(vc, identity, FUN.VALUE=matrix(0,q,q))
    if(is.null(dim(v.cl))) dim(v.cl) <- c(1,1,m)
    dimnames(v.cl)[1:2] <- lapply(dimnames(vc[[1]]), function(z) sub("^[(]Intercept[)]$","Intercept",z))
    vlist[[paste("|",clus,sep="")]] <- v.cl
  }else{
    vc <- lapply(model,nlme::VarCorr)
    # by variable
    for(vv in clus){
      q <- dim(model[[1]]$coefficients$random[[vv]])[2]
      if(length(clus)==1){
        rind <- 1:q
      }else{
        rind <- grep( paste0("^",vv," =$"), rownames(vc[[1]]))
        rind <- (rind+1):(rind+q)
      }
      # ... by term
      for(qq in rind){
        v.cl <- sapply(vc, function(x) as.numeric(x[qq,1]))
        if(is.null(dim(v.cl))) dim(v.cl) <- c(1,1,m)
        dimnames(v.cl)[1:2] <- list(sub("^[(]Intercept[)]$", "Intercept", rownames(vc[[1]])[qq]))
        vlist[[paste("|",vv,sep="")]] <- v.cl
      }
    }
  }
  # residual variance (if estimated)
  fixsigma <- attr(model[[1]]$modelStruct,"fixedSigma")
  if(!fixsigma){
    rv <- sapply(model, function(z) z$sigma^2)
    dim(rv) <- c(1,1,m)
    dimnames(rv) <- list("Residual","Residual",NULL)
    vlist <- c(vlist, list(rv))
  }

  # additional parameters
  # 1. ICC (only lme, single clustering)
  if(!fixsigma & cls=="lme" & length(clus)==1){
    if("(Intercept)"%in%colnames(vc[[1]])){
      iv <- sapply(vc, function(z) z["(Intercept)","(Intercept)"])
      icc <- iv / (iv + rv[1,1,])
      addp <- c(addp, mean(icc))
      names(addp) <- paste("ICC|",clus,sep="")
    }
  }

  list(vlist=vlist,addp=addp)
}

# *** geeglm method
.getVC.geeglm <- function(model){

  if(!requireNamespace("geepack", quietly=TRUE)) stop("The 'geepack' package must be installed in order to use this function.")
  m <- length(model)
  vlist <- addp <- NULL

  # variance components (currently not used)
  # vlist <- list()

  # additional parameters
  # 1. scale parameter (gamma)
  isfix <- model[[1]]$geese$model$scale.fix
  if(!isfix){
    gamma <- sapply(model, function(x) x$geese$gamma)
    if(is.null(dim(gamma))){
      dim(gamma) <- c(1,m)
      rownames(gamma) <- names(model[[1]]$geese$gamma)
    }
    addp <- c(addp,rowMeans(gamma))
    nms <- gsub("^[(]Intercept[)]$", "Intercept", names(addp))
    names(addp) <- paste0("Scale:",nms)
  }
  # 2. correlation parameters (alpha)
  corstr <- model[[1]]$geese$model$corstr
  isfix <- corstr%in%c("fixed","userdefined")
  if(!isfix){
    alpha <- sapply(model, function(x) x$geese$alpha)
    if(is.null(dim(alpha))){
      dim(alpha) <- c(1,m)
      rownames(alpha) <- names(model[[1]]$geese$alpha)
    }
    rownames(alpha) <- paste0("Correlation:",rownames(alpha))
    addp <- c(addp,rowMeans(alpha))
  }

  list(vlist=vlist,addp=addp)
}

# *** lm method
.getVC.lm <- function(model,ML=FALSE){

  m <- length(model)
  vlist <- addp <- NULL

  if(ML){  # SiG 16-04-2016
    rv <- sapply(model, function(z) sum(resid(z)^2)/length(resid(z)) )
  }else{
    rv <- sapply(model, function(z) sum(resid(z)^2)/df.residual(z) )
  }
  dim(rv) <- c(1,1,m)
  dimnames(rv) <- list("Residual","Residual",NULL)
  vlist <- c(vlist, list(rv))

  list(vlist=vlist,addp=addp)
}

