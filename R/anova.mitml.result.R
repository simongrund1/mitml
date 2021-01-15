anova.mitml.result <- function(object, ..., method = "D3", data = NULL, ariv = c("default", "positive", "robust")){

  # create list of models
  mod.list <- c(list(object), list(...))

  # ***
  # check input
  #

  # check lists
  m <- length(object)

  if(length(mod.list) == 1) stop("Comparison requires at least two lists of fitted statistical models.")
  if(any(!sapply(mod.list, is.list))) stop("The 'object' and '...' arguments must be lists of fitted statistical models.")
  if(any(sapply(mod.list[-1], length) != m)) stop("The 'object' and '...' arguments must be lists with the same length.")

  # check method
  method.choices <- c("D3", "D4", "D2")
  method <- original.method <- match.arg(method, method.choices)

  # check model classes
  cls.list <- lapply(mod.list, function(x) class(x[[1]]))

  if(any(sapply(cls.list[-1], "[", 1) != cls.list[[1]][1])) warning("The 'object' and '...' arguments appear to include objects of different classes. Results may not be trustworthy.")

  .checkNamespace(unique(unlist(cls.list)))

  # check for REML and refit (if needed)
  reml.list <- lapply(mod.list, function(x) sapply(x, .checkREML))
  reml <- any(unlist(reml.list))

  if(reml){
    for(ii in seq_along(mod.list)){
      mod.list[[ii]][reml.list[[ii]]] <- lapply(mod.list[[ii]][reml.list[[ii]]], .updateML)
    }
  }

  # ***
  # check method and possible fallback methods
  #

  # find user-defined method and possible fallback options
  try.method <- method.choices[seq.int(which(method.choices == original.method), length(method.choices))]
  error.msg <- character()

  # try logLik evaluation methods until working method is found
  for(mm in seq_along(try.method)){

    if(try.method[mm] == "D3") try.fun <- .evaluateUserLogLik
    if(try.method[mm] == "D4") try.fun <- .evaluateStackedLogLik
    if(try.method[mm] == "D2") try.fun <- .evaluateLogLik

    # check if method can be applied to specified objects
    res <- lapply(mod.list, function(x, fun){
      tryCatch(expr = suppressMessages(suppressWarnings(fun(x[1]))),
               error = function(e) e
      )
    }, fun = try.fun)


    # if applicable, proceed; otherwise, save error message and try next method (if any)
    notApplicable <- sapply(res, inherits, what = "error")
    if(any(notApplicable)){

      # save error message
      err <- as.character(res[[which(notApplicable)[1]]])
      error.msg[try.method[mm]] <- sub("^Error in .*: ", "", err)

      # try next method (if any)
      if(mm < length(try.method)){
        next()
      }else{
        stop("The '", original.method, "' method is not supported for the specified models, and no valid alternative was found. Problems were due to:\n", paste(error.msg, collapse = ""))
      }

    }else{

      # set method, print warning if needed
      method <- try.method[mm]
      if(method != original.method) warning("The '", original.method, "' method is not supported for the specified models. Switching to '", method, "'.")
      break()

    }

  }

  # ***
  # find order of models
  #

  # try to determine (numerator) degrees of freedom for each model
  df.list <- lapply(lapply(mod.list, "[[", 1), .getDFs)

  # check if models can be ordered
  reorderModels <- FALSE
  if(all(!sapply(df.list, is.null))){

    df.method <- sapply(df.list, attr, which = "type")

    # check if extraction method was consistent across models
    if(all(df.method[-1] == df.method[1])){
      reorderModels <- TRUE
    }

  }

  # re-order models (if possible)
  if(reorderModels){
    mod.list <- mod.list[order(unlist(df.list), decreasing = TRUE)]
  }else{
    warning("Could not determine the order of models in 'object' and '...'. The order is therefore assumed to be as specified (with decreasing complexity). Please check whether this was intended, and see '?testModels' for specific comparisons between models.")
  }

  # ***
  # perform model comparisons
  #

  # model comparisons
  nmod <- length(mod.list)
  out.list <- vector("list", nmod-1)

  for(ii in seq_len(nmod-1)){

    # make call
    cll <- call("testModels", model = quote(mod.list[[ii]]), null.model = quote(mod.list[[ii+1]]))
    cll[["method"]] <- method

    if(method == "D2") cll[["use"]] <- "likelihood"
    if(method == "D4"){
      if(!is.null(data)) cll[["data"]] <- data
      cll[["ariv"]] <- ariv
    }

    # evaluate call
    out.list[[ii]] <- eval(cll)

  }

  # try to get model formulas
  fml <- character(nmod)
  for(ii in seq_len(nmod)){

    f <- .getFormula(mod.list[[ii]][[1]])
    fml[ii] <- f

  }

  out <- list(
    call = match.call(),
    test = out.list,
    m = m,
    method = method,
    use = "likelihood",
    ariv = ariv,
    data = !is.null(data),
    formula = fml,
    order.method = ifelse(reorderModels, df.method[1], NULL),
    reml = reml
  )

  class(out) <- "mitml.anova"
  out

}

