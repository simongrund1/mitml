with.mitml.list <- function(data, expr, use.data = FALSE, ...){
# evaluates an expression for a list of data sets

  expr <- substitute(expr)
  pf <- parent.frame()

  out <- lapply(data, function(d, expr, pf, use.data){
    if(use.data){
      expr[["data"]] <- substitute(d)
      eval(expr, parent.frame())
    }else{
      eval(expr, d, pf)
    }
  }, expr = expr, pf = pf, use.data = use.data)

  class(out) <- c("mitml.result", "list")
  return(out)

}
