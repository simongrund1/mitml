confint.mitml.testEstimates <- function(object, parm, level=0.95, ...){
# calculate confidence intervals from pooled estimates

  est <- object$estimates

  pnames <- rownames(est)
  if(missing(parm)) parm <- pnames
  if(is.numeric(parm)) parm <- pnames[parm]

  cf <- est[parm,1]
  se <- est[parm,2]
  df <- est[parm,4]

  a <- (1-level)/2
  fac <- qt(1-a, est[parm,"df"])
  pct <- paste(format(100*c(a,1-a), trim=TRUE, scientific=FALSE, digits=3), "%")


  ci <- matrix(NA_real_, length(parm), 2, dimnames=list(parm,pct))
  ci[,1] <- cf - se*fac
  ci[,2] <- cf + se*fac
  
  return(ci)

}

