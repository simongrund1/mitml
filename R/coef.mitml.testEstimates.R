coef.mitml.testEstimates <- function(object, ...){
# extract pooled parameter estimates

  est <- object$estimates
  out <- est[, 1, drop = TRUE]
  if(is.null(names(out))) names(out) <- rownames(est)
  
  return(out)

}

