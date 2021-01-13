testEstimates <- function(model, qhat, uhat, var.comp = FALSE, df.com = NULL, ...){
# combine scalar estimates from the analysis of multiply imputed data

  # ***
  # check input
  #

  # check model specification
  if(missing(model) == (missing(qhat) || missing(uhat))){
    stop("Either 'model' or both 'qhat' and 'uhat' must be supplied.")
  }

  # check misc. parameters
  if(!var.comp) vc.out <- NULL

  # ***
  # process matrix, array or list arguments
  #

  if(!missing(qhat)){

    # check input
    if(missing(uhat)) stop("Either 'model' or both 'qhat' and 'uhat' must be supplied.")
    if(!is.matrix(qhat) && is.array(qhat)) stop("The 'qhat' argument must be either a matrix or a list.")

    # convert point estimates
    if(is.matrix(qhat)){
      qhat <- lapply(seq_len(ncol(qhat)), function(i, Q) Q[,i], Q = qhat)
    }
 
    # convert variance estimates
    if(is.matrix(uhat)){
      uhat <- lapply(seq_len(ncol(uhat)), function(i, U) U[,i], U = uhat)
    }else
    if(is.array(uhat)){
      uhat <- lapply(seq_len(dim(uhat)[3]), function(i, U) diag(as.matrix(U[,,i])), U = uhat)
    }

    # ensure proper format
    m <- length(qhat)
    if(m != length(uhat)) stop("Dimensions of 'qhat' and 'uhat' do not match.")

    Qhat <- matrix(unlist(qhat), ncol = m)
    Uhat <- matrix(unlist(uhat), ncol = m)
    if(any(!is.finite(Uhat))) stop("Missing values in 'uhat' are not allowed.")

    nms <- names(qhat[[1]])
    if(is.null(nms)) nms <- paste0("Parameter.", 1:nrow(Qhat))

    cls <- NULL
    vc.out <- NULL
    if(var.comp) warning("The 'var.comp' argument is ignored when 'qhat' and 'uhat' are used.")

  }

  # ***
  # process fitted models
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

    if(var.comp){

      vc.est <- .extractMiscParameters(model)
      vc.Qhat <- vc.est$Qhat
      vc.nms <- vc.est$nms

    }else{

      vc.Qhat <- vc.nms <- NULL

    }

  }

  # ***
  # pool results
  #

  Qbar <- apply(Qhat, 1, mean)
  Ubar <- apply(Uhat, 1, mean)
  B <- apply(Qhat, 1, var)
  T <- Ubar + (1+m^(-1)) * B

  se <- sqrt(T)
  t <- Qbar/se

  r <- (1+m^(-1))*B/Ubar

  # compute degrees of freedom
  v <- vm <- (m-1)*(1+r^(-1))^2
  if(!is.null(df.com)){
    lam <- r/(r+1)
    vobs <- (1-lam)*((df.com+1)/(df.com+3))*df.com
    v <- (vm^(-1)+vobs^(-1))^(-1)
  }

  fmi <- (r+2/(v+3))/(r+1)

  # create output for parameter estimates
  pval <- 2 * (1 - pt(abs(t), df = v))   # two-tailed p-value, SiG 2017-02-09
  out <- matrix(c(Qbar, se, t, v, pval, r, fmi), ncol = 7)
  colnames(out) <- c("Estimate", "Std.Error", "t.value", "df", "P(>|t|)", "RIV", "FMI") # two-tailed p-value, SiG 2017-02-09
  rownames(out) <- nms

  # create output for other parameter estimates
  if(var.comp && !missing(model)){

    if(is.null(vc.Qhat)){
      vc.out <- NULL
      warning("Computation of variance components not supported for objects of class '", paste(cls, collapse = "|"), "' (see ?with.mitml.list for manual calculation).")
    }else{
      vc.Qbar <- apply(vc.Qhat, 1, mean)
      vc.out <- matrix(vc.Qbar, ncol = 1)
      colnames(vc.out) <- "Estimate"
      rownames(vc.out) <- vc.nms
    }

  }

  out <- list(
    call = match.call(),
    estimates = out,
    var.comp = vc.out,
    m = m,
    adj.df = !is.null(df.com),
    df.com = df.com,
    cls.method = cls
  )

  class(out) <- "mitml.testEstimates"
  return(out)

}

