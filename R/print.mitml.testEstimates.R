print.mitml.testEstimates <- function(x, digits = 3, sci.limit = 5, ...){
# print method for MI estimates

  cll <- x$call
  est <- x$estimates
  ep <- x$extra.pars
  m <- x$m
  adj.df <- x$adj.df
  df.com <- x$df.com

  # print header
  cat("\nCall:\n", paste(deparse(cll)), sep = "\n")
  cat("\nFinal parameter estimates and inferences obtained from", m, "imputed data sets.\n")
  cat("\n")

  # print results
  if(!is.null(est)){

    # format numeric results
    pl <- attr(est, "par.labels")
    out <- .formatTable(est, digits = digits, sci.limit = sci.limit, labels = pl)
    for(i in seq_len(nrow(out))) cat(out[i,], "\n")

  }

  # print other results
  if(!is.null(ep)){

    if(!is.null(est)) cat("\n")

    # format numeric results
    pl <- attr(ep, "par.labels")
    out <- .formatTable(ep, digits = digits, sci.limit = sci.limit, labels = pl)
    for(i in seq_len(nrow(out))) cat(out[i,], "\n")

  }

  cat("\n")
  
  # print footer
  if(adj.df){
    cat(c("Hypothesis test adjusted for small samples with",
          paste0("df=[", paste(df.com, collapse = ","), "]\ncomplete-data degrees of freedom.")))
  }else{
    cat("Unadjusted hypothesis test as appropriate in larger samples.")
  }

  cat("\n\n")

  invisible()

}

summary.mitml.testEstimates <- function(object, ...){
# summary method for objects of class mitml.testEstimates

  print.mitml.testEstimates(object, ...)

}
