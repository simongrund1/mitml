mitml.list2mids <- function(x, data, fill = FALSE, where = NULL){
# convert objects of class "mitml.list" to "mids"

  # check for 'mice'
  if(!requireNamespace("mice", quietly = TRUE)) stop("The 'mice' package must be installed to use this function.")

  # check variable names
  nms.inc <- names(data)
  nms.imp <- unique(do.call(c, lapply(x, names)))

  if(any(c(".imp", ".id") %in% nms.inc)) stop("Columns named '.imp' or '.id' are not allowed in 'data'.")
  if(any(c(".imp", ".id") %in% nms.imp)) stop("Columns named '.imp' or '.id' are not allowed in 'x'.")

  nms.new <- nms.imp[!nms.imp %in% nms.inc]
  if(length(nms.new) > 0L){
    if(!fill) stop("Some variables in the imputed data ('x') are not present in the original data ('data') Use 'fill = TRUE' to include them.")
    data[, nms.new] <- NA
  }

  # prepare data
  z <- c(list(data), x)
  for(i in seq_along(z)){
    z[[i]] <- cbind(.imp = i - 1, .id = seq.int(1, nrow(z[[i]])), z[[i]])
  }

  return(mice::as.mids(long = do.call(rbind, z), where = where))

}

