print.mitml.testModels <- function(x, digits = 3, sci.limit = 5, ...){
# print method for MI estimates

  cll <- x$call
  test <- x$test
  method <- x$method
  use <- x$use
  reml <- x$reml
  m <- x$m
  data <- x$data
  ariv <- x$ariv
  adj.df <- x$adj.df
  df.com <- x$df.com

  # print header
  cat("\nCall:\n", paste(deparse(cll)), sep = "\n")
  cat("\nModel comparison calculated from", m, "imputed data sets.")

  # print method
  cat("\nCombination method:", method)
  if(method == "D2") cat(" (", use, ")", sep = "")
  if(method == "D4" && ariv == "robust") cat(" (robust)", sep = "")
  cat("\n\n")

  # print test results
  test.digits <- c(digits, 0, rep(digits, ncol(test)-2))
  out <- .formatTable(test, digits = test.digits, sci.limit = sci.limit)
  for(i in seq_len(nrow(out))) cat("  ", out[i,], "\n")

  cat("\n")

  # print footer (if any)
  footer <- FALSE

  if(method == "D1"){
    footer <- TRUE
    if(adj.df){
      cat("Hypothesis test adjusted for small samples with ",
          paste0("df=[", paste(df.com, collapse = ","), "]\ncomplete-data degrees of freedom."),
          "\n", sep = "")
    }else{
      cat("Unadjusted hypothesis test as appropriate in larger samples.\n")
    }
  }

  if(method == "D4"){
    footer <- TRUE
    if(data){
      cat("Data for stacking were extracted from the `data` argument.\n")
    }else{
      cat("Data for stacking were automatically extracted from the fitted models.\n")
    }
  }

  if(reml){
    footer <- TRUE
    cat("Models originally fit with REML were automatically refit using ML.\n")
  }

  if(footer) cat("\n")

  invisible()

}

summary.mitml.testModels <- function(object, ...){
# summary method for objects of class mitml.testModels

  print.mitml.testModels(object, ...)

}
