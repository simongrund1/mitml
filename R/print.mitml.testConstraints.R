print.mitml.testConstraints <- function(x, digits = 3, sci.limit = 5, ...){
# print method for MI estimates

  cll <- x$call
  test <- x$test
  constraints <- x$constraints
  method <- x$method
  m <- x$m
  adj.df <- x$adj.df
  df.com <- x$df.com

  # print header
  cat("\nCall:\n", paste(deparse(cll)), sep = "\n")
  cat("\nHypothesis test calculated from", m, "imputed data sets. The following\nconstraints were specified:\n\n")

  # print constrained estimates
  est <- cbind(x$Qbar, sqrt(diag(x$T)))
  colnames(est) <- c("Estimate", "Std. Error")
  rownames(est) <- paste0(constraints, ":")

  out <- .formatTable(est, digits = digits, sci.limit = sci.limit)
  for(i in seq_len(nrow(out))) cat("  ", out[i,], "\n")

  # print method
  cat("\nCombination method:", method, "\n\n")

  # print test results
  test.digits <- c(digits, 0, rep(digits, ncol(test)-2))
  out <- .formatTable(test, digits = test.digits, sci.limit = sci.limit)
  for(i in seq_len(nrow(out))) cat("  ", out[i,], "\n")

  # print footer
  if(method == "D1"){
    cat("\n")
    if(adj.df){
      cat(c("Hypothesis test adjusted for small samples with",
            paste0("df=[", paste(df.com, collapse = ","), "]\ncomplete-data degrees of freedom.")))
    }else{
      cat("Unadjusted hypothesis test as appropriate in larger samples.")
    }
  cat("\n")
  }

  cat("\n")

  invisible()

}

summary.mitml.testConstraints <- function(object, ...){
# summary method for objects of class mitml.testConstraints

  print.mitml.testConstraints(object, ...)

}

