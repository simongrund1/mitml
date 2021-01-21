print.mitml.anova <- function(x, digits = 3, sci.limit = 5, ...){
# print method for anova method

  cll <- x$call
  test <- x$test
  fml <- x$formula
  method <- x$method
  data <- x$data
  ariv <- x$ariv
  order.method <- x$order.method
  use <- x$use
  reml <- x$reml
  m <- x$test[[1]]$m

  n.tests <- length(fml)

  # print header
  cat("\nCall:\n", paste(deparse(cll)), sep = "\n")
  cat("\nModel comparison calculated from", m, "imputed data sets.")

  # print method
  cat("\nCombination method:", method)
  if(method == "D2") cat(" (", use, ")", sep = "")
  if(method == "D4" && ariv == "robust") cat(" (robust)", sep = "")
  cat("\n")

  # print model formulas
  cat("\n")
  for(mm in seq.int(1, n.tests)) cat("Model ", mm, ": ", fml[mm], "\n", sep = "")
  cat("\n")

  # combine multiple tests in one table
  test.tab <- lapply(test, "[[", "test")
  test.tab <- do.call(rbind, test.tab)
  rn <- paste0(seq.int(1, n.tests - 1), " vs ", seq.int(2, n.tests), " ")
  rownames(test.tab) <- rn

  # format table
  test.digits <- c(digits, 0, rep(digits, ncol(test.tab)-2))
  out <- .formatTable(test.tab, digits = test.digits, sci.limit = sci.limit)
  for(i in seq_len(nrow(out))) cat("  ", out[i,], "\n")

  cat("\n")

  # print footer
  if(is.null(order.method)){
    cat("Models were ordered as provided by the user (by decreasing complexity).\n")
  }else{
    cat("Models were automatically ordered via '", order.method, "' (by decreasing complexity).\n", sep = "")
  }

  if(method == "D4"){
    if(data){
      cat("Data for stacking were extracted from the `data` argument.\n")
    }else{
      cat("Data for stacking were automatically extracted from the fitted models.\n")
    }
  }

  if(reml){
    cat("Models originally fit with REML were automatically refit using ML.\n")
  }

  cat("\n")

  invisible()

}

