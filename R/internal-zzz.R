.localUpdate <- function(object, envir = parent.frame(), ...){
# update call in parent frame

  cll <- getCall(object)
  if (is.null(call)) stop("Need an object with a call component.")

  # update call components based on additional arguments (...)
  extras <- match.call(expand.dots = FALSE)$...
  for(i in names(extras)) cll[[i]] <- extras[[i]]

  # update in local environment
  eval(cll, envir = envir)

}

.checkDeprecated <- function(x, arg.list, name){
# match argument list (arg.list, usually ...) by name against deprecated (name)
# and return matching value if match is found, otherwise return original value
# (x)

  cll <- match.call()

  nms <- names(arg.list)
  m <- sapply(nms, function(n, o){
    m <- try(match.arg(n, o), silent = TRUE)
    return(if(inherits(m, "try-error")) NA else m)
  }, o = name)

  # if match is found, print message and assign value to new name
  if(any(!is.na(m))){
    ans <- arg.list[[nms[1]]]
    msg <- paste0("The '", name, "' argument is deprecated. Please use '", as.character(cll[[2]]), "' instead.")
    warning(msg)
  }else{
    ans <- x
  }

  return(ans)

}

.checkNamespace <- function(x){
# check required packages for supported object types

  # specify class-package pairs
  cls.pkg <- list(
    "lme4" = "^g?l?merMod$",
    "nlme" = "^n?lme$",
    "geepack" = "^geeglm$",
    "survival" = "^coxph$",
    "MASS" = "^polr$"
  )

  # match class to package names
  req.pkg <- lapply(cls.pkg, function(p, x) grep(pattern = p, x = x, value = TRUE), x = x)
  req.pkg <- req.pkg[sapply(req.pkg, length) > 0]

  for(i in seq_along(req.pkg)){
    pkg.name <- names(req.pkg)[i]
    pkg.cls <- paste(req.pkg[[i]], collapse = "|")
    if(!requireNamespace(pkg.name, quietly = TRUE)) stop("The '", pkg.name, "' package must be installed in order to use this function with objects of class '", pkg.cls, "'.")
  }

  invisible(NULL)

}

.formatTable <- function(x, prefix = "%.", postfix = "f", digits = 3, sci.limit = 5, width,
                         col.names, row.names, labels = NULL, labels.sep = 3){
# format table with common format and fixed width

  # row and column names
  if(missing(col.names)) col.names <- colnames(x)
  if(missing(row.names)) row.names <- rownames(x)

  # fotmat
  fmt <- paste0(prefix, digits, postfix)
  if(ncol(x) %% length(fmt)) stop("Format and table dimensions do not match.")
  fmt <- rep_len(fmt, length.out = ncol(x))

  # format for large values
  isLarge <- apply(x, 2, function(z, a) any(z >= 10^a), a = sci.limit)
  fmt[isLarge] <- sub(paste0(postfix, "$"), "e", fmt[isLarge])

  # make formatted matrix
  y <- matrix("", nrow(x), ncol(x))
  for(i in seq_len(ncol(x))) y[,i] <- sprintf(fmt[i], x[,i] + 0)

  # find width
  if(missing(width)) width <- max(sapply(c(colnames(x), y), nchar))

  # fill table
  out <- matrix("", nrow(x)+1, ncol(x)+1)
  out[,1] <- format(c("", row.names), justify = "left")
  out[1, -1] <- format(col.names, justify = "right", width = width)
  out[-1, -1] <- format(y, justify = "right", width = width)

  # add labels (if any)
  if(!is.null(labels)){
    labels[nchar(labels) > 0] <- paste0("(", labels[nchar(labels) > 0], ")")
    pl <- format(labels, justify = "left")
    nc <- max(nchar(pl))
    out[-1, 1] <- paste0(out[-1, 1], paste0(rep(" ", labels.sep), collapse = ""), pl)
    out[1, 1] <- paste0(out[1, 1], paste0(rep(" ", nc + labels.sep), collapse = ""))
  }

  return(out)

}


.extractMatrix <- function(x, ...){
# extract submatrix from array (indexed by ...)

  if(is.null(dim(x))) return(x)

  out <- `[`(x, , , ...)
  dim(out) <- dim(x)[1:2]
  dimnames(out) <- dimnames(x)[1:2]

  out

}

.adiag <- function(x, stacked = FALSE){
# extract diagonal elements of first two dimensions in three-dimensional array
# containing either square (default) or stacked matrices

  d <- dim(x)

  # indices for diagonal entries (square or stacked-square)
  if(stacked){
    i <- seq_len(d[2]) + d[1]*(seq_len(d[2])-1)
    i <- outer(i, (seq_len(d[1]/d[2])-1)*d[2], `+`)
    i <- outer(i, (seq_len(d[3])-1)*d[1]*d[2], `+`)
  }else{
    i <- seq_len(d[1]) + d[1]*(seq_len(d[1])-1)
    i <- outer(i, (seq_len(d[3])-1)*d[1]^2, `+`)
  }

  x[as.vector(i)]

}

