summary.mitml <- function(object, n.Rhat=3, goodness.of.appr=FALSE,
                          autocorrelation=FALSE, ...){

# summary method for objects of class "mitml"

  inc <- object$data
  ngr <- length(unique(attr(object$data,"group")))
  prm <- object$par.imputation
  iter <- dim(prm[[1]])[3]
  k <- object$iter$iter
  isL2 <- attr(object$model,"is.L2")

  # parameter chains (for backwards compatibility)
  if(is.null(object$keep.chains)) object$keep.chains <- "full"

  # percent missing
  mdr <- sapply(inc, FUN=function(x){mean(is.na(x))})
  mdr[] <- sprintf(mdr*100,fmt="%.1f")
  mdr <- gsub("^0.0$","0",mdr)

  # convergence for imputation phase
  conv <- NULL
  Rhat <- ifelse(is.null(n.Rhat), FALSE, n.Rhat >= 2)

  SDprop <- goodness.of.appr
  ACF <- autocorrelation
  if(Rhat|SDprop|ACF){

    conv <- c(list(beta=NULL), if(isL2) list(beta2=NULL), list(psi=NULL, sigma=NULL))
    for(pp in names(conv)){

      ni <- dim(prm[[pp]])[1]
      nj <- dim(prm[[pp]])[2]
      nl <- dim(prm[[pp]])[4]
      cmat <- matrix(NA_real_, ni*nj*nl, 3+Rhat+SDprop+3*ACF)
      cmat[,1] <- rep(1:ni,nj*nl)
      cmat[,2] <- rep(1:nj,each=ni,times=nl)
      cmat[,3] <- rep(1:nl,each=ni*nj)
      colnames(cmat) <- c("i1", "i2", "grp", if(Rhat) "Rhat", if(SDprop) "SDprop",
                          if(ACF) c("lag-1","lag-k","lag-2k"))

      for(ll in 1:nl){ # by group

        for(jj in 1:nj){
          for(ii in 1:ni){

            # check for redundant entries
            if(pp=="psi"){
              if(jj > ii) next
            }
            if(pp=="sigma"){
              if(jj > ((ii-1)%%nj)+1) next
            }
            ind <- ( cmat[,1]==ii & cmat[,2]==jj & cmat[,3]==ll )
            chn <- matrix(prm[[pp]][ii,jj,,ll], 1, iter)
            # potential scale reduction (Rhat)
            if(Rhat) cmat[ind,"Rhat"] <- .GelmanRubin(chn,n.Rhat)
            # goodness of approximation
            if(SDprop) cmat[ind,"SDprop"] <- .SDprop(chn)
            # autocorrelation
            if(ACF){
              cmat[ind,"lag-1"] <- .reducedACF(chn, lag=1, smooth=0)
              cmat[ind,"lag-k"] <- .reducedACF(chn, lag=k, smooth=2, sd=.5)
              cmat[ind,"lag-2k"] <- .reducedACF(chn, lag=2*k, smooth=2, sd=.5)
            }
          }
        }
      }
      conv[[pp]] <- cmat[!apply(cmat,1,function(x) any(is.na(x))),,drop=F]
    }

  attr(conv,"stats") <- c("Rhat","SDprop","ACF")[c(Rhat,SDprop,ACF)]
  }

  smr <- list(
    call=object$call,
    model=object$model,
    prior=object$prior,
    iter=object$iter,
    keep.chains=object$keep.chains,
    ngr=ngr,
    missing.rates=mdr,
    conv=conv
  )

  class(smr) <- "mitml.summary"
  smr

}

.reducedACF <- function(x, lag, smooth=0, sd=.5){

  # check NA
  if(all(is.na(x))) return(NA)

  n <- length(x)
  lag0 <- lag
  lag <- lag + (-smooth:smooth)

  ac <- numeric(length(lag))
  y <- x - mean(x)
  ss.y <- sum(y^2)

  for(li in 1:length(lag)){
    ll <- lag[li]
    # leave at 0 for constant value
    ac[li] <- if(ss.y>0) sum( y[1:(n-ll)] * y[1:(n-ll)+ll] ) / ss.y else 0
  }

  if(smooth>0){
    # weights based on normal density
    w <- dnorm(-smooth:smooth, 0, sd)
    ac <- sum( ac * (w/sum(w)) )
  }

  ac

}

