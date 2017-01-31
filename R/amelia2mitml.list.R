amelia2mitml.list <- function(x){
# convert amelia to mitml.list

  out <- unname(x$imputations)
  class(out) <- c("mitml.list","list")
  out

}

