sort.mitml.list <- function(x, decreasing = FALSE, by, ...){
# sort list of multiply imputed data sets

  expr <- substitute(by)
  args0 <- list(decreasing = decreasing, ...)

  res <- lapply(x, function(i){
    args <- eval(expr, i, parent.frame())
    if(!is.list(args)) args <- list(args)
    ind <- do.call("order", c(args, args0))
    i[ind,]
  })

  as.mitml.list(res)
  
}
