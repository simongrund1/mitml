mids2mitml.list <- function(x){
# convert mids to mitml.list

  if(!requireNamespace("mice", quietly = TRUE)) stop("The 'mice' package must be installed to use this function.")
  m <- x$m

  out <- list()
  length(out) <- m
  for(ii in 1:m){
    out[[ii]] <- mice::complete(x, action = ii)
  }

  class(out) <- c("mitml.list", "list")
  return(out)

}

