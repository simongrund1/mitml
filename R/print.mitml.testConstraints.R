print.mitml.testConstraints <- function(x,...){
# print method for MI estimates

  cl <- x$call
  test <- x$test
  cons <- x$constraints
  mth <- x$method
  m <- x$m
  adj <- x$adj.df
  dfc <- x$df.com

  # header
  cat("\nCall:\n", paste(deparse(cl)), sep="\n")
  cat("\nHypothesis test calculated from",m,"imputed data sets. The following\nconstraints were specified:\n\n")

  # print constraint table
  est <- matrix(c(x$Qbar, sqrt(diag(x$T))), ncol=2)
  colnames(est) <- c("Estimate", "Std. Error")
  rownames(est) <- paste0(cons, ":")

  out <- .formatTable.helper(est)
  for(i in 1:nrow(out)) cat("  ", out[i,],"\n")

  cat("\nCombination method:",mth,"\n")

  # print test table
  fmt <- c("%.3f","%.0f","%.3f","%.3f","%.3f")
  fmt[test>=10^5] <- "%.3e" # large values
  out <- sprintf(fmt,test)

  # table
  cat("\n")
  w <- max(sapply(c(out,colnames(test)),nchar))
  cat("  ",format(colnames(test),justify="right",width=w),"\n")
  cat("  ",format(out,justify="right",width=w),"\n")

  if(mth=="D1"){
  cat(if(adj){c("\nHypothesis test adjusted for small samples with",
              paste("df=[",paste(dfc,collapse=","),"]\ncomplete-data degrees of freedom.",sep=""))
      }else{"\nUnadjusted hypothesis test as appropriate in larger samples."},"\n")
  }

  cat("\n")

  invisible()
}

summary.mitml.testConstraints <- function(object,...){
# summary method for objects of class mitml.testConstraints

  print.mitml.testConstraints(object,...)

}

.formatTable.helper <- function(x){

  f <- sprintf("%.3f",x)
  w <- max(sapply(c(colnames(x),f),nchar))
  out <- matrix("",nrow(x)+1,ncol(x)+1)
  out[,1] <- format(c("",rownames(x)))
  out[1,-1] <- format(colnames(x),justify="right",width=w)
  out[-1,-1] <- format(f,justify="right",width=w)

  return(out)
}

