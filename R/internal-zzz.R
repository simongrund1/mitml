.extractMatrix <- function(x, ...){
# extract submatrix from array (indexed by ...)

  if(is.null(dim(x))) return(x)

  out <- `[`(x,,,...)
  dim(out) <- dim(x)[1:2]
  dimnames(out) <- dimnames(x)[1:2]

  out

}

.adiag <- function(x, stacked=FALSE){
# extract diagonal elements of first two dimensions in three-dimensional array
# containing either square (default) or stacked matrices

  d <- dim(x)

  # indices for diagonal entries (square or stacked-square)
  if(stacked){
    i <- seq_len(d[2]) + d[1]*(seq_len(d[2])-1)
    i <- outer(i,(seq_len(d[1]/d[2])-1)*d[2],`+`)
    i <- outer(i,(seq_len(d[3])-1)*d[1]*d[2],`+`)
  }else{
    i <- seq_len(d[1]) + d[1]*(seq_len(d[1])-1)
    i <- outer(i,(seq_len(d[3])-1)*d[1]^2,`+`)
  }

  x[as.vector(i)]

}

.check.modelL2 <- function(x){
# check number of model equations

  if(is.list(x) & length(x)>2)
    stop("Cannot determine the number of levels. The 'formula' or 'type' argument must indicate either a single model for responses at level 1, or two models for responses at level 1 and 2.")

  if(is.list(x) & length(x)==1) x <- x[[1]] # unlist
  isL2 <- is.list(x) & length(x)==2

  attr(x,"is.L2") <- isL2
  x

}

.check.variablesL2 <- function(x,clus){
# check for variables at L2 (constant at L1)

  apply(x, 2, function(a) all( abs(a-clusterMeans(a,clus)) < sqrt(.Machine$double.eps),
                               na.rm=T))

}

