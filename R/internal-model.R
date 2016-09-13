# prepare model input by formula
.model.byFormula <- function(data, formula, group, group.original,
                      method=c("pan","jomo","jomo.matrix")){

  # L2: separate model equations
  formula  <- .check.modelL2(formula)
  isL2 <- attr(formula,"is.L2")
  if(isL2){
    formula.L2 <- formula[[2]]
    formula <- formula[[1]]
  }

  method <- match.arg(method)

  # *** evaluate L1 model
  #

  ft <- terms(formula)
  tl <- attr(ft,"term.labels")
  vrs <- attr(ft,"variables")[-1]

  # responses
  yvrs <- as.character(vrs)[attr(ft,"response")]
  yvrs <- gsub("[\r\n]","",yvrs)
  y.fml <- as.formula(paste0("~",yvrs))
  yvrs <- attr(terms(y.fml), "term.labels")
  # check for untransformed yvrs
  err <- !(yvrs %in% colnames(data))
  if(any(err)) stop("Could not find: ", paste0(yvrs[err],collapse=", "), "). Target variables must be contained in the data set 'as is', and transformations must be applied beforehand.")

  # cluster id
  clt <- tl[grep("\\|",tl)]
  if(length(clt)==0) stop("Cluster indicator not found in formula\n\n",formula,"\n\nPlease specify the cluster indicator and at least one random term using the '|' operator.")
  clt <- strsplit( clt, split="[[:blank:]]*\\|[[:blank:]]*" )[[1]]
  clname <- clt[2]

  # order data
  data <- data[ order(group,data[,clname]), ]
  group.original <- group.original[ order(group) ]
  group <- group[ order(group) ]

  # predictors: fixed
  pvrs <- c(if(attr(ft,"intercept")){"(Intercept)"}, tl[!grepl("\\|",tl)])
  fe.fml <- c(if(attr(ft,"intercept")){"1"}else{"0"}, tl[!grepl("\\|",tl)])
  fe.fml <- as.formula(paste0("~", paste0(fe.fml,collapse="+")))
  # predictors: random
  cl.fml <- as.formula(paste("~",clt[1]))
  cl.ft <- terms(cl.fml)
  qvrs <- c(if(attr(cl.ft,"intercept")){"(Intercept)"}, attr(cl.ft,"term.labels"))

  # model matrix for fe and cl
  attr(data,"na.action") <- identity
  mmp <- suppressWarnings( model.matrix(fe.fml, data=data) )
  mmq <- suppressWarnings( model.matrix(cl.fml, data=data) )
  pnames <- colnames(mmp)
  qnames <- colnames(mmq)
  psave <- setdiff( c(pnames,qnames), c("(Intercept)",colnames(data)) )

  switch( method ,
    pan={ # panImpute (matrix input)
      y <- data.matrix(data[yvrs])
      ycat <- NULL
    },
    jomo={ # jomoImpute, for higher-level functions (data input)
      y <- data[yvrs]
      ycat <- NULL
    },
    jomo.matrix={ # jomoImpute, for lower-level versions (matrix input)
      y <- data.matrix(data[yvrs])
      cvrs <- sapply(data[,yvrs,drop=F], is.factor)
      ycat <- y[,cvrs,drop=F]
      y <- y[,!cvrs,drop=F]
    }
  )

  clus <- data[,clname]
  pred <- cbind(mmp, mmq[,!(qnames%in%pnames),drop=F])
  xcol <- which(colnames(pred)%in%pnames)
  zcol <- which(colnames(pred)%in%qnames)

  # assign to parent.frame
  inp <- list(
    y=y, ycat=ycat, clus=clus, pred=pred, xcol=xcol, zcol=zcol, data=data,
    group=group, group.original=group.original, psave=psave, clname=clname,
    yvrs=yvrs, pvrs=pvrs, qvrs=qvrs, pnames=pnames, qnames=qnames
  )
  for(i in names(inp)) assign(i, inp[[i]], pos=parent.frame())

  # *** evaluate L2 model
  #

  if(isL2){

    ft <- terms(formula.L2)
    tl <- attr(ft,"term.labels")
    vrs <- attr(ft,"variables")[-1]

    # responses
    yvrs <- as.character(vrs)[attr(ft,"response")]
    yvrs <- gsub("[\r\n]","",yvrs)
    y.fml <- as.formula(paste0("~",yvrs))
    yvrs <- attr(terms(y.fml), "term.labels")
    # check for untransformed yvrs
    err <- !(yvrs %in% colnames(data))
    if(any(err)) stop("Could not find: ", paste0(yvrs[err],collapse=", "), "). Target variables must be contained in the data set 'as is', and transformations must be applied beforehand.")

    # predictors: fixed only at L2
    pvrs <- c(if(attr(ft,"intercept")){"(Intercept)"}, tl[!grepl("\\|",tl)])
    fe.fml <- c(if(attr(ft,"intercept")){"1"}else{"0"}, tl[!grepl("\\|",tl)])
    fe.fml <- as.formula(paste0("~", paste0(fe.fml,collapse="+")))

    # model matrix for FE only
    attr(data,"na.action") <- identity
    mmp <- suppressWarnings( model.matrix(fe.fml, data=data) )
    pnames <- colnames(mmp)
    psave <- c( psave, setdiff( c(pnames), c("(Intercept)",colnames(data)) ) )

    switch( method ,
      jomo={ # jomoImpute, for higher-level functions (data input)
        y <- data[yvrs]
        ycat <- NULL
      },
      jomo.matrix={ # jomoImpute, for lower-level versions (matrix input)
        y <- data.matrix(data[yvrs])
        cvrs <- sapply(data[,yvrs,drop=F], is.factor)
        ycat <- y[,cvrs,drop=F]
        y <- y[,!cvrs,drop=F]
      }
    )

    pred <- mmp
    xcol <- which(colnames(pred) %in% pnames)

    # assign to parent.frame
    inp <- list(
      y.L2=y, ycat.L2=ycat, pred.L2=pred, xcol.L2=xcol, yvrs.L2=yvrs,
      pvrs.L2=pvrs, pnames.L2=pnames, psave=psave
    )
    for(i in names(inp)) assign(i, inp[[i]], pos=parent.frame())

  }

  invisible(NULL)

}


# prepare model input by type
.model.byType <- function(data, type, group, group.original,
                   method=c("pan","jomo","jomo.matrix")){

  # L2: separate model equations
  type <- .check.modelL2(type)
  isL2 <- attr(type,"is.L2")
  if(isL2){
    type.L2 <- type[[2]]
    type <- type[[1]]
  }

  # *** evaluate L1 model
  #

  if(ncol(data)!=length(type)) stop("Length of 'type' must be equal to the number of colums in 'data'.")
  if(sum(type==-2)<1) stop("Cluster indicator not found.")
  if(sum(type==-2)>1) stop("Only one cluster indicator may be specified.")

  data <- data[ order(group,data[,type==-2]), ]
  group.original <- group.original[ order(group) ]
  group <- group[ order(group) ]

  clname <- colnames(data)[type==-2]
  clus <- data[,clname]
  yvrs <- colnames(data)[type==1]

  switch( method ,
    pan={ # panImpute (matrix input)
      y <- data.matrix(data[yvrs])
      ycat <- NULL
    },
    jomo={ # jomoImpute, newer versions (data input)
      y <- data[yvrs]
      ycat <- NULL
    },
    jomo.matrix={ # jomoImpute, older versions (matrix input)
      y <- data.matrix(data[yvrs])
      cvrs <- sapply(data[,yvrs,drop=F], is.factor)
      ycat <- y[,cvrs,drop=F]
      y <- y[,!cvrs,drop=F]
    }
  )

  pred <- cbind(1,as.matrix(data[type%in%c(2,3)]))
  pvrs <- c("(Intercept)",colnames(data)[type%in%c(2,3)])
  qvrs <- c("(Intercept)",colnames(data)[type==3])
  colnames(pred) <- pvrs

  xcol <- 1:length(pvrs)
  zcol <- xcol[pvrs%in%qvrs]

  # assign to parent frame
  inp <- list(
    y=y, ycat=ycat, clus=clus, pred=pred, xcol=xcol, zcol=zcol, data=data,
    group=group, group.original=group.original, clname=clname, yvrs=yvrs,
    pvrs=pvrs, qvrs=qvrs, pnames=pvrs, qnames=qvrs
  )
  for(i in names(inp)) assign(i, inp[[i]], pos=parent.frame())

  # *** evaluate L2 model
  #

  if(isL2){

    if(ncol(data)!=length(type.L2)) stop("Length of 'type' must be equal to the number of colums in 'data'.")

    yvrs <- colnames(data)[type.L2==1]

    switch( method ,
      pan={ # panImpute (matrix input)
        y <- data.matrix(data[yvrs])
        ycat <- NULL
      },
      jomo={ # jomoImpute, newer versions (data input)
        y <- data[yvrs]
        ycat <- NULL
      },
      jomo.matrix={ # jomoImpute, older versions (matrix input)
        y <- data.matrix(data[yvrs])
        cvrs <- sapply(data[,yvrs,drop=F], is.factor)
        ycat <- y[,cvrs,drop=F]
        y <- y[,!cvrs,drop=F]
      }
    )

    pred <- cbind(1,as.matrix(data[type.L2%in%c(2,3)]))
    pvrs <- c("(Intercept)",colnames(data)[type.L2%in%c(2,3)])
    colnames(pred) <- pvrs

    xcol <- 1:length(pvrs)

    # assign to parent frame
    inp <- list(
      y.L2=y, ycat.L2=ycat, pred.L2=pred, xcol.L2=xcol, yvrs.L2=yvrs,
      pvrs.L2=pvrs, pnames.L2=pvrs
    )
    for(i in names(inp)) assign(i, inp[[i]], pos=parent.frame())

  }

  invisible(NULL)

}


