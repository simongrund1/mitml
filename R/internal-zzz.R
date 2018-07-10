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

