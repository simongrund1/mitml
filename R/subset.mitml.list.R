subset.mitml.list <- function(x, subset, select, ...){
# subset list of multiply imputed data sets
# NOTE: code adapted from subset.data.frame (by Peter Dalgaard and Brian Ripley)

  rind <- if (missing(subset)) {
    lapply(x, function(i) rep(TRUE, nrow(i)))
  } else {
    ss <- substitute(subset)
    rind <- lapply(x, function(i) eval(ss, i, parent.frame()))
    if (!is.logical(unlist(rind))) stop("'subset' must be logical")
    lapply(rind, function(i) i & !is.na(i))
  }

  cind <- if (missing(select)) {
    lapply(x, function(i) TRUE)
  } else {
    nl <- lapply(x, function(i){
      l <- as.list(seq_along(i))
      names(l) <- names(i)
      l
    })
    se <- substitute(select)
    lapply(nl, function(i) eval(se, i, parent.frame()))
  }

  res <- lapply(seq_along(x), function(i) x[[i]][rind[[i]], cind[[i]], drop = FALSE])
  as.mitml.list(res)
  
}

