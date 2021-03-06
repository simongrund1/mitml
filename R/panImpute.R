panImpute <- function(data, type, formula, n.burn = 5000, n.iter = 100, m = 10,
                      group = NULL, prior = NULL, seed = NULL, save.pred = FALSE,
                      keep.chains = c("full", "diagonal"), silent = FALSE){

# wrapper function for the Gibbs sampler in the pan package

  # *** checks
  if(!missing(type) && !missing(formula)) stop("Only one of 'type' or 'formula' may be specified.")
  if(save.pred && !missing(type)){
    warning("Option 'save.pred' is ignored if 'type' is specified")
    save.pred = FALSE
  }
  keep.chains <- match.arg(keep.chains)

  # convert type
  if(!missing(type)){
    formula <- .type2formula(data, type)
    group <- attr(formula, "group")
  }

  # empty objects to assign to
  clname <- yvrs <- y <- ycat <- zcol <- xcol <- pred <- clus <- psave <-
    pvrs <- qvrs <- pnames <- qnames <- NULL

  # preserve original order
  if(!is.data.frame(data)) as.data.frame(data)
  data <- cbind(data, original.order = 1:nrow(data))

  # address additional grouping
  grname <- group
  if(is.null(group)){
    group <- rep(1, nrow(data))
  }else{
    group <- data[,group]
    if(length(group) != nrow(data)) stop("Argument 'group' is not correctly specified.")
  }
  group.original <- group
  group <- as.numeric(factor(group, levels = unique(group)))

  # ***
  # model input

  # populate local frame
  .model.byFormula(data, formula, group, group.original, method = "pan")

  # check model input
  if(any(is.na(group)))
    stop("Grouping variable must not contain missing data.")
  if(any(is.na(pred)))
    stop("Predictor variables must not contain missing data.")
  if(sum(is.na(y)) == 0)
    stop("Target variables do not contain any missing data.")
  if(any(!sapply(y, is.numeric)))
    stop("Target variables must be numeric. You may either convert them or use jomoImpute() instead.")
  if(any(duplicated(yvrs)))
    stop("Found duplicate target variables.")

  # reorder colums
  cc <- which(colnames(data) %in% c(clname, grname, yvrs))
  data.ord <- cbind(data[c(clname, grname, yvrs)], data[-cc])

  # ***
  # pan setup

  if(is.null(prior)){
    prior <- list( a = ncol(y), Binv = diag(1, ncol(y)),
      c = ncol(y)*length(zcol), Dinv = diag(1, ncol(y)*length(zcol)) )
  }

  if(is.null(seed)){
    set.seed(as.integer(runif(1, 0, 10^6)))
  }else{
    set.seed(as.integer(seed))
  }
  rns <- sapply(unique(group), function(x, m) as.integer(runif(m+1, 0, 10^6)), m = m)

  # prepare output
  ind <- which(is.na(data.ord), arr.ind = TRUE, useNames = FALSE)
  ind <- ind[ ind[,2] %in% which(colnames(data.ord) %in% colnames(y)), ,drop = FALSE ]
  rpm <- matrix(NA, nrow(ind), m)

  # standard dimensions
  ng <- length(unique(group))
  np <- length(xcol)
  nq <- length(zcol)
  nr <- ncol(y)

  # reduced dimensions
  dpsi <- nr*nq
  dsig <- nr
  if(keep.chains == "diagonal"){
    dpsi <- dsig <- 1
  }

  bpar <- list(beta = array( NA, c(np, nr, n.burn, ng) ),
               psi = array( NA, c(nr*nq, dpsi, n.burn, ng) ),
               sigma = array( NA, c(nr, dsig, n.burn, ng) ))
  ipar <- list(beta = array( NA, c(np, nr, n.iter*m, ng) ),
               psi = array( NA, c(nr*nq, dpsi, n.iter*m, ng) ),
               sigma = array( NA, c(nr, dsig, n.iter*m, ng) ))

  # burn-in
  if(!silent){
    cat("Running burn-in phase ...\n")
    flush.console()
  }
  glast <- as.list(unique(group))
  for(gg in unique(group)){

    gi <- group == gg
    gy <- y[gi,]
    gpred <- pred[gi,]
    gclus <- clus[gi]
    # sort 1, ..., k
    gclus <- match(gclus, unique(gclus))

    cur <- pan::pan(gy, subj = gclus, gpred, xcol, zcol, prior, seed = rns[1, gg], iter = n.burn)
    glast[[gg]] <- cur$last

    # save parameter chains
    bpar[["beta"]][,,,gg] <- cur$beta
    if(keep.chains == "diagonal"){
      bpar[["psi"]][,,,gg] <- .adiag( cur$psi )
      bpar[["sigma"]][,,,gg] <-.adiag( cur$sigma )
    }else{
      bpar[["psi"]][,,,gg] <- cur$psi
      bpar[["sigma"]][,,,gg] <- cur$sigma
    }

  }

  # imputation
  for(ii in 1:m){
    if(!silent){
      cat("Creating imputed data set (", ii, "/", m,") ...\n")
      flush.console()
    }

    gy.imp <- as.list(unique(group))
    for(gg in unique(group)){

      gi <- group == gg
      gy <- y[gi,]
      gpred <- pred[gi,]
      gclus <- clus[gi]
      # sort 1, ..., k
      gclus <- match(gclus, unique(gclus))

      cur <- pan::pan(gy, subj = gclus, gpred, xcol, zcol, prior, seed = rns[ii+1, gg], iter = n.iter,
        start = glast[[gg]])
      glast[[gg]] <- cur$last

      # save imputations
      gy.imp[[gg]] <- cur$y

      # save parameter chains
      i0 <- seq.int(n.iter*(ii-1)+1, n.iter*ii)
      ipar[["beta"]][,,i0, gg] <- cur$beta
      if(keep.chains == "diagonal"){
        ipar[["psi"]][,,i0, gg] <- .adiag( cur$psi ) 
        ipar[["sigma"]][,,i0, gg] <- .adiag( cur$sigma )
      }else{
        ipar[["psi"]][,,i0, gg] <- cur$psi
        ipar[["sigma"]][,,i0, gg] <- cur$sigma
      }

    }
    y.imp <- do.call(rbind, gy.imp)
    rpm[,ii] <- y.imp[is.na(y)]

  }

  if(!silent){
    cat("Done!\n")
  }

  # clean up
  srt <- data.ord[,ncol(data.ord)]
  data.ord <- data.ord[,-ncol(data.ord)]

  # prepare output data
  if( save.pred && !missing(formula) ) data.ord <- cbind(data.ord, pred[, psave, drop = F])
  # ordering
  attr(data.ord, "sort") <- srt
  attr(data.ord, "group") <- group.original
  # model summary
  model <- list(clus = clname, yvrs = yvrs, pvrs = pvrs, qvrs = qvrs)
  attr(model, "is.ML") <- TRUE
  attr(model, "is.L2") <- FALSE
  attr(model, "full.names") <- list(pvrs = pnames, qvrs = qnames)

  out <- list(
    data = data.ord,
    replacement.mat = rpm,
    index.mat = ind,
    call = match.call(),
    model = model,
    random.L1 = "none",
    prior = prior,
    iter = list(burn = n.burn, iter = n.iter, m = m),
    keep.chains = keep.chains,
    par.burnin = bpar,
    par.imputation = ipar
  )
  class(out) <- c("mitml", "pan")
  return(out)

}

