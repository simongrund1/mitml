# prepare model input by formula
.model.byFormula <- function(data, formula, group, group.original,
                      method = c("pan", "jomo", "jomo.matrix")){

  # check model, separate equations
  formula  <- .check.model(formula)

  isML <- attr(formula, "is.ML")
  isL2 <- attr(formula, "is.L2")

  if(isL2){
    formula.L2 <- formula[[2]]
    formula <- formula[[1]]
  }

  method <- match.arg(method)

  # *** evaluate L1 model
  #

  ft <- terms(formula)
  tl <- attr(ft, "term.labels")
  vrs <- attr(ft, "variables")[-1]
  nms <- colnames(data)

  # responses
  yvrs <- as.character(vrs)[attr(ft, "response")]
  yvrs <- gsub("[\r\n]", "", yvrs)
  y.fml <- as.formula(paste0("~", yvrs))
  yvrs <- attr(terms(y.fml), "term.labels")
  # check for untransformed yvrs
  err <- !(yvrs %in% nms)
  if(any(err)) stop("Could not find: ", paste0(yvrs[err], collapse = ", "), "). Target variables must be contained in the data set 'as is', and transformations must be applied beforehand.")

  # cluster id
  clt <- tl[grep("\\|", tl)]

  if(method == "pan" & !isML) stop("Cluster indicator not found in formula\n\n", .formula2char(formula), "\n\nPlease specify the cluster indicator and at least one random term using the '|' operator. Single-level imputation is supported by jomoImpute().")

  # extract and reorder
  if(isML){

    clt <- strsplit( clt, split = "[[:blank:]]*\\|[[:blank:]]*" )[[1]]
    clname <- clt[2]

    # order data and grouping
    data <- data[ order(group, data[,clname]), ]
    group.original <- group.original[ order(group) ]
    group <- group[ order(group) ]

  }else{
    clname <- NULL
  }

  # predictors: fixed
  pvrs <- c(if(attr(ft, "intercept")){"(Intercept)"}, tl[!grepl("\\|", tl)])
  fe.fml <- c(if(attr(ft, "intercept")){"1"}else{"0"}, tl[!grepl("\\|", tl)])
  fe.fml <- as.formula(paste0("~", paste0(fe.fml, collapse = "+")))

  # predictors: random
  if(isML){
    cl.fml <- as.formula(paste("~", clt[1]))
    cl.ft <- terms(cl.fml)
    qvrs <- c(if(attr(cl.ft, "intercept")){"(Intercept)"}, attr(cl.ft, "term.labels"))
  }else{
    cl.fml <- ~0
    qvrs <- NULL
  }

  # model matrix for fe and cl
  attr(data, "na.action") <- identity
  mmp <- suppressWarnings( model.matrix(fe.fml, data = data) )
  mmq <- suppressWarnings( model.matrix(cl.fml, data = data) )
  pnames <- colnames(mmp)
  qnames <- colnames(mmq)
  psave <- setdiff( c(pnames, qnames), c("(Intercept)", nms) )

  switch( method ,
    # panImpute (matrix input)
    pan={
      y <- data.matrix(data[yvrs])
      ycat <- NULL
    },
    # jomoImpute, for higher-level functions (data frames, uses jomo for preprocessing)
    jomo={
      y <- data[yvrs]
      ycat <- NULL
    },
    # jomoImpute, for higher- and lower-level versions (preprocessed matrix input)
    jomo.matrix={
      y <- data.matrix(data[yvrs])
      cvrs <- sapply(data[, yvrs, drop = F], is.factor)
      ycat <- y[,cvrs, drop = F]
      y <- y[,!cvrs, drop = F]
    }
  )

  clus <- if(isML) data[,clname] else NULL
  pred <- cbind(mmp, mmq[,!(qnames%in%pnames), drop = F])
  xcol <- which(colnames(pred)%in%pnames)
  zcol <- which(colnames(pred)%in%qnames)

  # assign to parent.frame
  inp <- list(
    y = y, ycat = ycat, clus = clus, pred = pred, xcol = xcol, zcol = zcol, data = data,
    group = group, group.original = group.original, psave = psave, clname = clname,
    yvrs = yvrs, pvrs = pvrs, qvrs = qvrs, pnames = pnames, qnames = qnames
  )

  for(i in names(inp)) assign(i, inp[[i]], pos = parent.frame())

  # *** evaluate L2 model
  #

  if(isL2){

    ft <- terms(formula.L2)
    tl <- attr(ft, "term.labels")
    vrs <- attr(ft, "variables")[-1]

    # responses
    yvrs <- as.character(vrs)[attr(ft, "response")]
    yvrs <- gsub("[\r\n]", "", yvrs)
    y.fml <- as.formula(paste0("~", yvrs))
    yvrs <- attr(terms(y.fml), "term.labels")
    # check for untransformed yvrs
    err <- !(yvrs %in% nms)
    if(any(err)) stop("Could not find: ", paste0(yvrs[err], collapse = ", "), "). Target variables must be contained in the data set 'as is', and transformations must be applied beforehand.")

    # predictors: fixed only at L2
    pvrs <- c(if(attr(ft, "intercept")){"(Intercept)"}, tl[!grepl("\\|", tl)])
    fe.fml <- c(if(attr(ft, "intercept")){"1"}else{"0"}, tl[!grepl("\\|", tl)])
    fe.fml <- as.formula(paste0("~", paste0(fe.fml, collapse = "+")))

    # model matrix for FE only
    attr(data, "na.action") <- identity
    mmp <- suppressWarnings( model.matrix(fe.fml, data = data) )
    pnames <- colnames(mmp)
    psave <- c( psave, setdiff( c(pnames), c("(Intercept)", nms) ) )

    switch( method ,
      jomo={ # jomoImpute, for higher-level functions (data input)
        y <- data[yvrs]
        ycat <- NULL
      },
      jomo.matrix={ # jomoImpute, for lower-level versions (matrix input)
        y <- data.matrix(data[yvrs])
        cvrs <- sapply(data[,yvrs, drop = F], is.factor)
        ycat <- y[,cvrs, drop = F]
        y <- y[,!cvrs, drop = F]
      }
    )

    pred <- mmp
    xcol <- which(colnames(pred) %in% pnames)

    # assign to parent.frame
    inp <- list(
      y.L2 = y, ycat.L2 = ycat, pred.L2 = pred, xcol.L2 = xcol, yvrs.L2 = yvrs,
      pvrs.L2 = pvrs, pnames.L2 = pnames, psave = psave
    )

    for(i in names(inp)) assign(i, inp[[i]], pos = parent.frame())

  }

  invisible(NULL)

}

# convert formula to character
.formula2char <- function(x){

  chr <- as.character(x)
  paste(chr[c(3, 1, 2)])

}

.check.model <- function(x){
# check model type and number of levels

  xnew <- x

  # ensure proper list format
  if(is.list(x) & length(x) > 2) stop("Cannot determine the number of levels. The 'formula' or 'type' argument must indicate either a single-level model, a model for responses at level 1, or two models for responses at level 1 and 2.")

  if(!is.list(x)) x <- list(x)

  # check cluster specification and model type
  clt <- lapply(x, function(z){
    if(is.language(z)){
      tl <- attr(terms(z), "term.labels")
      tl[grep("\\|", tl)]
    }else{
      which(z == -2)
    }
  })
  isML <- length(clt[[1]]) > 0
  isL2 <- length(x) == 2

  if(isL2 & !isML) stop("No cluster variable found. Imputation models for responses at level 1 and 2 require the specification of a cluster variable in the level-1 equation.")

  attr(xnew, "is.ML") <- isML
  attr(xnew, "is.L2") <- isL2
  xnew

}

.check.variablesL2 <- function(x, clus){
# check for variables at L2 (constant at L1)

  apply(x, 2, function(a) all( abs(a-clusterMeans(a, clus)) < sqrt(.Machine$double.eps),
                               na.rm = T))

}

# convert type to formula
.type2formula <- function(data, type){

  # L2: separate model equations
  type <- .check.model(type)
  isML <- attr(type, "is.ML")
  isL2 <- attr(type, "is.L2")
  if(isL2){
    type.L2 <- type[[2]]
    type <- type[[1]]
  }

  nms <- colnames(data)

  # grouping
  grp <- if(any(type == -1)) nms[type == -1] else NULL
  if(isL2 & is.null(grp)){
    if(any(type.L2 == -1)) grp <- nms[type.L2 == -1]
  }

  # L1 model
  if(ncol(data) != length(type)) stop("Length of 'type' must be equal to the number of colums in 'data'.")
  if(sum(type == -2)>1) stop("Only one cluster indicator may be specified.")
   
  cls <- nms[type == -2]

  yvrs <- paste( nms[type == 1], collapse = "+" )
  pvrs <- paste( c(1, nms[type%in%c(2, 3)]), collapse = "+" )
  qvrs <- if(isML) paste( c(1, nms[type == 3]), collapse = "+" ) else NULL
  
  # build L1 formula
  cls.fml <- if(isML) paste("+ (", qvrs, "|", cls, ")") else NULL
  fml <- formula( paste(yvrs, "~", pvrs, cls.fml) )

  # L2 model
  if(isL2){

    if(ncol(data) != length(type.L2)) stop("Length of 'type' must be equal to the number of colums in 'data'.")

    yvrs <- paste( nms[type.L2 == 1], collapse = "+" )
    pvrs <- paste( c(1, nms[type.L2%in%c(2, 3)]), collapse = "+" )

    # build formula (make list)
    fml <- list( fml, formula( paste(yvrs, "~", pvrs) ) )

  }

  attr(fml, "group") <- grp
  attr(fml, "is.ML") <- isML
  attr(fml, "is.L2") <- isL2

  return(fml)

}

