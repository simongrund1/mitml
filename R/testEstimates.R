testEstimates <- function(model, qhat, uhat = NULL, extra.pars = FALSE, df.com = NULL, ...){
# combine scalar estimates from the analysis of multiply imputed data

  # ***
  # check input
  #

  # handle deprecated arguments
  dots <- list(...)
  extra.pars <- .checkDeprecated(extra.pars, arg.list = dots, name = "var.comp")

  # check model specification
  if(missing(model) == missing(qhat)){
    stop("Either 'model' or 'qhat' must be supplied.")
  }

  # check extra parameters
  if(!extra.pars) pooled.ep.est <- NULL

  # ***
  # process matrix, array or list input
  #

  if(!missing(qhat)){

    # check arguments
    if(extra.pars) warning("The 'extra.pars' argument is ignored when 'qhat' is used.")

    # check qhat
    if(!is.matrix(qhat) && is.array(qhat)) stop("The 'qhat' argument must be either a matrix or a list.")
    if(is.matrix(qhat)){
      qhat <- lapply(seq_len(ncol(qhat)), function(i, Q) Q[,i], Q = qhat)
    }
 
    # check uhat
    if(!is.null(uhat)){

      # convert formats
      if(is.matrix(uhat)){
        uhat <- lapply(seq_len(ncol(uhat)), function(i, U) U[,i], U = uhat)
      }else
      if(is.array(uhat)){
        uhat <- lapply(seq_len(dim(uhat)[3]), function(i, U) diag(as.matrix(U[,,i])), U = uhat)
      }

    }

    # convert to common format
    m <- length(qhat)
    Qhat <- matrix(unlist(qhat), ncol = m)
    if(is.null(uhat)){
      Uhat <- NULL
    }else{
      Uhat <- matrix(unlist(uhat), ncol = m)
      if(any(!is.finite(Uhat))) stop("Missing values in 'uhat' are not allowed.")
    }

  # get parameter names
    nms <- names(qhat[[1]])
    if(is.null(nms)) nms <- paste0("Parameter.", 1:nrow(Qhat))

    cls <- NULL
    pooled.ep.est <- NULL

  }

  # ***
  # process model input
  #

  if(!missing(model)){

    if(!is.list(model)) stop("The 'model' argument must be a list of fitted statistical models.")
    m <- length(model)

    # get class (check required packages)
    cls <- class(model[[1]])
    .checkNamespace(cls)

    # extract parameter estimates
    est <- .extractParameters(model, diagonal = TRUE, include.extra.pars = TRUE)

    Qhat <- est$Qhat
    Uhat <- est$Uhat
    nms <- est$nms

    if(extra.pars){

      ep.est <- .extractMiscParameters(model)
      ep.Qhat <- ep.est$Qhat
      ep.nms <- ep.est$nms

    }else{

      ep.Qhat <- ep.nms <- NULL

    }

  }

  # ***
  # pool results
  #

  # pool estimates
  pooled.est <- .pool.estimates(Qhat = Qhat, Uhat = Uhat, m = m, df.com = df.com, par.labels = nms)

  # pool estimates of extra parameters
  if(extra.pars && !missing(model)){

    if(is.null(ep.Qhat)){
      pooled.ep.est <- NULL
      warning("Computation of variance components not supported for objects of class '", paste(cls, collapse = "|"), "' (see ?with.mitml.list for manual calculation).")
    }else if(length(ep.Qhat) == 0){
      pooled.ep.est <- NULL
    }else{
      pooled.ep.est <- .pool.estimates(Qhat = ep.Qhat, Uhat = NULL, par.labels = ep.nms)
    }

  }

  out <- list(
    call = match.call(),
    estimates = pooled.est,
    extra.pars = pooled.ep.est,
    m = m,
    adj.df = !is.null(df.com),
    df.com = df.com,
    cls.method = cls
  )

  class(out) <- "mitml.testEstimates"
  return(out)

}

