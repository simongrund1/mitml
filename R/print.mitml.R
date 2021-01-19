print.mitml <- function(x, ...){
# print method for objects of class "mitml"

  cl <- x$call
  vrs <-x$model 
  itr <- x$iter
  ngr <- length(unique(attr(x$data, "group")))
  isML <- attr(x$model, "is.ML")
  isL2 <- attr(x$model, "is.L2")

  cat("\nCall:\n", paste(deparse(cl)), sep = "\n")
  cat("\n")

  if(isL2) cat("Level 1:\n", collapse = "\n")
  if(isML) cat(formatC("Cluster variable:", width=-25), vrs$clus, sep = " ", collapse = "\n")
  cat(formatC("Target variables:", width=-25), vrs$yvrs, collapse = "\n")
  cat(formatC("Fixed effect predictors:", width=-25), vrs$pvrs, collapse = "\n")
  if(isML) cat(formatC("Random effect predictors:", width=-25), vrs$qvrs, collapse = "\n")

  if(isL2){
    cat("\n")
    cat(formatC("Level 2:\n", width=-25), collapse = "\n")
    cat(formatC("Target variables:", width=-25), vrs$yvrs.L2, collapse = "\n")
    cat(formatC("Fixed effect predictors:", width=-25), vrs$pvrs.L2, collapse = "\n")
  }

  cat("\nPerformed", sprintf("%.0f", itr$burn), "burn-in iterations, and generated", sprintf("%.0f", itr$m),
      "imputed data sets,\neach", sprintf("%.0f", itr$iter), "iterations apart.",
      if(ngr>1){c("\nImputations were carried out seperately within", sprintf("%.0f", ngr), "groups.\n")}, "\n")

  invisible(NULL)

}
