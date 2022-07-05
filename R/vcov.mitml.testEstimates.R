vcov.mitml.testEstimates <- function(object, ...){
# extract pooled variance-covariance matrix of parameter estimates

  est <- object$estimates
  v <- attr(est, "T")

  # full variance-covariance matrix
  if(!is.null(v)){
    out <- v
    rownames(out) <- colnames(out) <- rownames(est)
  } else
  # diagonal (squared SEs)
  if(any(colnames(est) == "Std.Error")){
    warning("Could find only diagonal elements of the pooled variance-covariance matrix.")
    p <- nrow(est)
    out <- matrix(NA_real_, nrow = p, ncol = p)
    diag(out) <- est[, "Std.Error", drop = TRUE]^2
    rownames(out) <- colnames(out) <- rownames(est)
  }else{
    stop("Could not find the pooled variance-covariance matrix.")
  }
  
  return(out)

}

