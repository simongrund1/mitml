jomoImpute <- function(data, type, formula, random.L1=c("none","mean","full"),
                       n.burn=5000, n.iter=100, m=10, group=NULL, prior=NULL,
                       seed=NULL, save.pred=FALSE,
                       keep.chains=c("full","diagonal"), silent=FALSE){

# wrapper function for the different samplers of the jomo package

  # checks arguments
  if(!missing(type) & !missing(formula)) stop("Only one of 'type' or 'formula' may be specified.")
  if(save.pred & !missing(type)){
    warning("Option 'save.pred' is ignored if 'type' is specified")
    save.pred=FALSE
  }
  random.L1 <- match.arg(random.L1)
  keep.chains <- match.arg(keep.chains)

  # convert type
  if(!missing(type)){
    formula <- .type2formula(data,type)
    group <- attr(formula, "group")
  }
  # check for number of model equations
  formula  <- .check.modelL2( formula )
  isL2 <- attr(formula,"is.L2")

  # objects to assign to
  clname <- yvrs <- y <- ycat <- zcol <- xcol <- pred <- clus <- psave <- pvrs <-
  qvrs <- pnames <- qnames <- yvrs.L2 <- y.L2 <- ycat.L2 <- xcol.L2 <- pred.L2 <-
  pvrs.L2 <- pnames.L2 <- NULL

  # preserve original order
  if(!is.data.frame(data)) as.data.frame(data)
  data <- cbind(data, original.order=1:nrow(data))

  # address additional grouping
  grname <- group
  if(is.null(group)){
    group <- rep(1,nrow(data))
  }else{
    if(!group%in%colnames(data)) stop("Argument 'group' is not correctly specified.")
    group <- data[,group]
  }
  group.original <- group
  group <- as.numeric(factor(group,levels=unique(group)))

  # ***
  # model input

  # populate local frame
  .model.byFormula(data, formula, group, group.original, method="jomo.matrix")

  # check model input
  if(any(is.na(group)))
    stop("Grouping variable must not contain missing data.")
  if(any(is.na(pred)))
    stop("Predictor variables must not contain missing data.")
  if(any(!sapply(data[yvrs], function(a) is.factor(a) || is.numeric(a))))
    stop("Target variables must either be numeric or factors.")
  if((sum(is.na(y)) + sum(is.na(ycat)) + ifelse(isL2, sum(is.na(y.L2))+sum(is.na(ycat.L2)), 0))==0)
    stop("Target variables do not contain any missing data.")
  if(any(duplicated(c(yvrs,yvrs.L2))))
    stop("Found duplicate target variables.")
  if(isL2){
    if(any(is.na(pred.L2)))
      stop("Predictor variables must not contain missing data.")
    if(any(!sapply(data[yvrs.L2], function(a) is.factor(a) || is.numeric(a))))
      stop("Target variables must either be numeric or factors.")
  }

  # check for L1 variables in L2 models
  if(isL2){
    y.L1 <- !.check.variablesL2(y.L2,clus)
    x.L1 <- !.check.variablesL2(pred.L2,clus)

    if(any(y.L1)) stop("Target variables at level 1 are not allowed in level-2 equation.")
    if(any(x.L1)){
      for(i in which(x.L1)) pred.L2[,i] <- clusterMeans(pred.L2[,i],clus)
      message("NOTE: Predictor variables at level 1 were found in level-2 equation and were replaced with cluster means (", paste0(pvrs.L2[x.L1], collapse=", "), ").")
    }
  }

  # reorder colums
  cc <- which(colnames(data) %in% c(clname,grname,yvrs,yvrs.L2))
  data.ord <- cbind(data[c(clname,grname,yvrs,yvrs.L2)],data[-cc])

  # *** jomo setup
  #

  ycat.labels <- lapply(data[,c(colnames(ycat),colnames(ycat.L2)),drop=F], levels)
  # select function
  func <- if(ncol(ycat)==0) "con" else if(ncol(y)==0) "cat" else "mix"
  func <- get( paste0(ifelse(!isL2,"jomo1ran","jomo2"), if(!isL2) func,
                      if(isL2 & random.L1=="none") "com",
                      if(random.L1!="none") "hr", ".MCMCchain"), asNamespace("jomo"))

  # standard dimensions and data properties
  ng <- length(unique(group))
  np <- length(xcol)
  nq <- length(zcol)
  ncon <- ncol(y)
  ncat <- ncol(ycat)
  nr <- ncon + ncat   # combined con + cat (variables)
  ynumcat <- matrix(0,ng,ncat)
  nc <- nr2 <- integer(ng)
  if(isL2){
    np.L2 <- length(xcol.L2)
    ncon.L2 <- ncol(y.L2)
    ncat.L2 <- ncol(ycat.L2)
    nr.L2 <- ncon.L2 + ncat.L2   # combined con + cat (variables)
    ynumcat.L2 <- matrix(0,ng,ncat.L2)
    nc.L2 <- nr2.L2 <- integer(ng)
  }else{
    nr2.L2 <- integer(ng)    # zero counts for compatibility
    ncon.L2 <- ncat.L2 <- 0  # of shared code
  }

  # ... manage categories groupwise
  for(gg in unique(group)){

    ynumcat[gg,] <- apply(ycat[group==gg,,drop=F], 2,
                          FUN=function(x) length(unique(x[!is.na(x)])))
    nc[gg] <- length(unique(clus[group==gg]))
    nr2[gg] <- ncon+sum(ynumcat[gg,])-length(ynumcat[gg,]) # combined con + cat (indicators)

    if(isL2){
      ynumcat.L2[gg,] <- apply(ycat.L2[group==gg,,drop=F], 2,
                               FUN=function(x) length(unique(x[!is.na(x)])))
      nc.L2[gg] <- length(unique(clus[group==gg]))
      nr2.L2[gg] <- ncon.L2+sum(ynumcat.L2[gg,])-length(ynumcat.L2[gg,])
    }

  }

  # reduced dimensions
  dpsi <- max(nr2)*nq+max(nr2.L2)
  dsig1 <- ifelse(random.L1=="full", max(nr2)*max(nc), max(nr2))
  dsig2 <- max(nr2)
  if(keep.chains=="diagonal"){
    dpsi <- dsig2 <- 1
  }

  # * * * * * * * * * * * * * * * * * * * *

  # save original seed (if seed is provided)
  original.seed <- NULL
  if(!is.null(seed)){
    if(exists(".Random.seed", .GlobalEnv)) original.seed <- .Random.seed
    set.seed(seed)
  }

  # priors
  if(is.null(prior)){
    prior <- as.list(unique(group))
    for(gg in unique(group)){
      prior[[gg]] <- list( Binv=diag(1,nr2[gg]),
                           Dinv=diag(1,nq*nr2[gg]+nr2.L2[gg]) )
      if(random.L1!="none") prior[[gg]]$a <- nr2[gg]
    }
  }else{ # check if prior is given as simple list
    if(!is.list(prior[[1]])) prior <- rep(list(prior),ng)
  }

  # prepare output
  ind <- which(is.na(data.ord), arr.ind=TRUE, useNames=FALSE)
  ind <- ind[ ind[,2] %in% which(colnames(data.ord)%in%c(yvrs,yvrs.L2)),,drop=FALSE ]
  rpm <- matrix(NA, nrow(ind), m)

  bpar <- c(list(beta=array( NA, c(np,max(nr2),n.burn,ng) )),
            if(isL2) list(beta2=array( NA, c(np.L2,max(nr2.L2),n.burn,ng) )),
            list(psi=array( NA, c(max(nr2)*nq+max(nr2.L2),dpsi,n.burn,ng) ),
                 sigma=array( NA, c(dsig1,dsig2,n.burn,ng) )))
  ipar <- c(list(beta=array( NA, c(np,max(nr2),n.iter*m,ng) )),
            if(isL2) list(beta2=array( NA, c(np.L2,max(nr2.L2),n.iter*m,ng) )),
            list(psi=array( NA, c(max(nr2)*nq+max(nr2.L2),dpsi,n.iter*m,ng) ),
                 sigma=array( NA, c(dsig1,dsig2,n.iter*m,ng) )))

  # burn-in
  if(!silent){
    cat("Running burn-in phase ...\n")
    flush.console()
  }
  glast <- as.list(unique(group))
  for(gg in unique(group)){

    gi <- group==gg
    gprior <- prior[[gg]]

    # function arguments (group specific)
    gclus <- clus[gi]
    gclus <- matrix( match(gclus, unique(gclus))-1, ncol=1 )
    func.args <- list( Y=if(ncon>0 & ncat==0 & !isL2) y[gi,,drop=F] else NULL,
                       Y.con=if(ncon>0 & (ncat>0 | isL2)) y[gi,,drop=F] else NULL,
                       Y.cat=if(ncat>0) ycat[gi,,drop=F] else NULL,
                       Y.numcat=if(ncat>0) ynumcat[gg,] else NULL,
                       Y2.con=if(ncon.L2>0) y.L2[gi,,drop=F] else NULL,
                       Y2.cat=if(ncat.L2>0) ycat.L2[gi,,drop=F] else NULL,
                       Y2.numcat=if(ncat.L2>0) ynumcat.L2[gg,] else NULL,
                       X=pred[gi,xcol,drop=F],
                       X2=if(isL2) pred.L2[gi,xcol.L2,drop=F] else NULL,
                       Z=pred[gi,zcol,drop=F],
                       clus=gclus,
                       beta.start=matrix(0,np,nr2[gg]),
                       l2.beta.start=if(isL2) matrix(0,np.L2,nr2.L2[gg]) else NULL,
                       u.start=matrix(0,nc[gg],nq*nr2[gg]+nr2.L2[gg]),
                       l1cov.start=if(random.L1!="none"){
                         matrix(diag(1,nr2[gg]),nr2[gg]*nc[gg],nr2[gg],byrow=T)
                       }else{
                         diag(1,nr2[gg])
                       },
                       l2cov.start=diag(1,nq*nr2[gg]+nr2.L2[gg]),
                       start.imp=NULL,
                       l2.start.imp=NULL,
                       l1cov.prior=gprior$Binv,
                       l2cov.prior=gprior$Dinv,
                       a=gprior$a,
                       meth=if(random.L1!="none") "random" else NULL,
                       nburn=n.burn,
                       output=0
    )
    func.args <- func.args[!sapply(func.args,is.null)]

    cur <- do.call( func, func.args )
    glast[[gg]] <- cur

    # current parameter dimensions (group-specific)
    bdim <- dim(cur$collectbeta)[1:2]
    pdim <- dim(cur$collectcovu)[1:2]
    sdim <- dim(cur$collectomega)[1:2]

    # save parameter chains
    bpar[["beta"]][1:bdim[1],1:bdim[2],,gg] <- cur$collectbeta
    if(keep.chains=="diagonal"){
      bpar[["psi"]][1:pdim[1],1,,gg] <- .adiag(cur$collectcovu)
    }else{
      bpar[["psi"]][1:pdim[1],1:pdim[2],,gg] <- cur$collectcovu
    }
    # ... random covariance matrices at L1
    if(random.L1=="mean"){
      tmp <- cur$collectomega
      dim(tmp) <- c(nr2[gg],nc[gg],nr2[gg],n.burn)
      if(keep.chains=="diagonal"){
        bpar[["sigma"]][1:sdim[2],1,,gg] <- .adiag(apply(tmp,c(1,3,4),mean))
      }else{
        bpar[["sigma"]][1:sdim[2],1:sdim[2],,gg] <- apply(tmp,c(1,3,4),mean)
      }
    }else{
      if(keep.chains=="diagonal"){
        bpar[["sigma"]][1:sdim[1],1,,gg] <- .adiag(cur$collectomega,
          stacked=(random.L1=="full"))
      }else{
        bpar[["sigma"]][1:sdim[1],1:sdim[2],,gg] <- cur$collectomega
      }
    }
    # ... L2 model
    if(isL2){
      bdim2 <- dim(cur$collect.l2.beta)[1:2]
      bpar[["beta2"]][1:bdim2[1],1:bdim2[2],,gg] <- cur$collect.l2.beta
    }

  }

  # imputation
  for(ii in 1:m){
    if(!silent){
      cat("Creating imputed data set (",ii,"/",m,") ...\n")
      flush.console()
    }

    gy.imp <- as.list(unique(group))
    for(gg in unique(group)){

      gi <- group==gg
      gprior <- prior[[gg]]

      # last state (imp)
      last.imp <- if(isL2 | ncat>0) glast[[gg]]$finimp.latnorm else glast[[gg]]$finimp
      if(ncon>0 & ncat==0 & !isL2)
        last.imp <- last.imp[(nrow(y[gi,,drop=F])+1):nrow(last.imp), 1:ncon, drop=F]
      last.imp.L2 <- if(isL2) glast[[gg]]$l2.finimp.latnorm else NULL

      # function arguments (group specific)
      gclus <- clus[gi]
      gclus <- matrix( match(gclus, unique(gclus))-1, ncol=1 )
      it <- dim(glast[[gg]]$collectbeta)[3]
      func.args <- list( Y=if(ncon>0 & ncat==0 & !isL2) y[gi,,drop=F] else NULL,
                         Y.con=if(ncon>0 & (ncat>0 | isL2)) y[gi,,drop=F] else NULL,
                         Y.cat=if(ncat>0) ycat[gi,,drop=F] else NULL,
                         Y.numcat=if(ncat>0) ynumcat[gg,] else NULL,
                         Y2.con=if(ncon.L2>0) y.L2[gi,,drop=F] else NULL,
                         Y2.cat=if(ncat.L2>0) ycat.L2[gi,,drop=F] else NULL,
                         Y2.numcat=if(ncat.L2>0) ynumcat.L2[gg,] else NULL,
                         X=pred[gi,xcol,drop=F],
                         X2=if(isL2) pred.L2[gi,xcol.L2,drop=F] else NULL,
                         Z=pred[gi,zcol,drop=F],
                         clus=gclus,
                         beta.start=.extractMatrix(glast[[gg]]$collectbeta,it),
                         l2.beta.start=.extractMatrix(glast[[gg]]$collect.l2.beta,it),
                         u.start=.extractMatrix(glast[[gg]]$collectu,it),
                         l1cov.start=.extractMatrix(glast[[gg]]$collectomega,it),
                         l2cov.start=.extractMatrix(glast[[gg]]$collectcovu,it),
                         start.imp=last.imp,
                         l2.start.imp=last.imp.L2,
                         l1cov.prior=gprior$Binv,
                         l2cov.prior=gprior$Dinv,
                         a=gprior$a,
                         meth=if(random.L1!="none") "random" else NULL,
                         nburn=n.iter,
                         output=0
      )
      func.args <- func.args[!sapply(func.args,is.null)]

      cur <- do.call( func, func.args )
      glast[[gg]] <- cur

      # save imputations
      ri <- (nrow(gclus)+1):nrow(cur$finimp)
      ci <- which(colnames(cur$finimp) %in% c(yvrs,yvrs.L2))
      gy.imp[[gg]] <- cur$finimp[ri,ci,drop=F]

      # current parameter dimensions (group-specific)
      bdim <- dim(cur$collectbeta)[1:2]
      pdim <- dim(cur$collectcovu)[1:2]
      sdim <- dim(cur$collectomega)[1:2]

      # save parameter chains
      iind <- (n.iter*(ii-1)+1):(n.iter*ii)
      ipar[["beta"]][1:bdim[1],1:bdim[2],iind,gg] <- cur$collectbeta
      if(keep.chains=="diagonal"){
        ipar[["psi"]][1:pdim[1],1,iind,gg] <- .adiag(cur$collectcovu)
      }else{
        ipar[["psi"]][1:pdim[1],1:pdim[2],iind,gg] <- cur$collectcovu
      }
      # ... random covariance matrices at L1
      if(random.L1=="mean"){
        tmp <- cur$collectomega
        dim(tmp) <- c(nr2[gg],nc[gg],nr2[gg],n.iter)
        if(keep.chains=="diagonal"){
          ipar[["sigma"]][1:sdim[2],1,iind,gg] <- .adiag(apply(tmp,c(1,3,4),mean))
        }else{
          ipar[["sigma"]][1:sdim[2],1:sdim[2],iind,gg] <- apply(tmp,c(1,3,4),mean)
        }
      }else{
        if(keep.chains=="diagonal"){
          ipar[["sigma"]][1:sdim[1],1,iind,gg] <- .adiag(cur$collectomega,
            stacked=(random.L1=="full"))
        }else{
          ipar[["sigma"]][1:sdim[1],1:sdim[2],iind,gg] <- cur$collectomega
        }
      }
      # ... L2 model
      if(isL2){
        bdim2 <- dim(cur$collect.l2.beta)[1:2]
        ipar[["beta2"]][1:bdim2[1],1:bdim2[2],iind,gg] <- cur$collect.l2.beta
      }

    }
    y.imp <- data.matrix(do.call(rbind,gy.imp))
    rpm[,ii] <- y.imp[,c(yvrs,yvrs.L2)][is.na(data.ord[,c(yvrs,yvrs.L2),drop=F])]

  }

  if(!silent){
    cat("Done!\n")
  }

  # clean up
  srt <- data.ord[,ncol(data.ord)]
  data.ord <- data.ord[,-ncol(data.ord)]

  # restore original seed (if seed was provided)
  if(!is.null(seed)){
    if(is.null(original.seed)){
      rm(".Random.seed", envir = .GlobalEnv)
    }else{
      assign(".Random.seed", original.seed, envir=.GlobalEnv)
    }
  }

  # *** prepare output
  #

  # save pred
  if( save.pred & !missing(formula) ){
    ps1 <- colnames(pred) %in% psave
    ps2 <- (colnames(pred.L2) %in% psave) & !(colnames(pred.L2) %in% colnames(pred)[ps1])
    data.ord <- cbind(data.ord, pred[,ps1,drop=F])
    if(isL2) cbind(data.ord, pred.L2[,ps2,drop=F])
  }
  # ordering
  attr(data.ord,"sort") <- srt
  attr(data.ord,"group") <- group.original
  # categorical variables
  if(ncat>0 | ncat.L2>0){
    attr(data.ord,"cvrs") <- names(ycat.labels)
    attr(data.ord,"levels") <- cbind(ynumcat,if(isL2) ynumcat.L2)
    attr(data.ord,"labels") <- ycat.labels
  }
  # model summary
  model <- list(clus=clname, yvrs=yvrs, pvrs=pvrs, qvrs=qvrs,
                yvrs.L2=if(isL2) yvrs.L2 else NULL,
                pvrs.L2=if(isL2) pvrs.L2 else NULL)
  attr(model,"is.L2") <- isL2
  attr(model,"full.names") <- list(pvrs=pnames, qvrs=qnames,
                                   pvrs.L2=if(isL2) pnames.L2 else NULL)

  out <- list(
    data=data.ord,
    replacement.mat=rpm,
    index.mat=ind,
    call=match.call(),
    model=model,
    random.L1=random.L1,
    prior=prior,
    iter=list(burn=n.burn, iter=n.iter, m=m),
    keep.chains=keep.chains,
    par.burnin=bpar,
    par.imputation=ipar
  )
  class(out) <- c("mitml","jomo")
  out

}

