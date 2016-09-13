print.mitml.summary <- function(x,...){
# print method for objects of class "summary.mitml"

  cl <- x$call
  vrs <- x$model
  itr <- x$iter
  ngr <- x$ngr
  mdr <- x$missing.rates
  conv <- x$conv
  isL2 <- attr(x$model,"is.L2")

  # print general information
  cat("\nCall:\n", paste(deparse(cl)), sep="\n")
  cat("\n")

  if(isL2) cat("Level 1:\n", collapse="\n")
  cat(formatC("Cluster variable:",width=-25), vrs$clus, sep=" ", collapse="\n")
  cat(formatC("Target variables:",width=-25), vrs$yvrs, collapse="\n")
  cat(formatC("Fixed effect predictors:",width=-25), vrs$pvrs, collapse="\n")
  cat(formatC("Random effect predictors:",width=-25), vrs$qvrs, collapse="\n")

  if(isL2){
    cat("\n")
    cat(formatC("Level 2:\n",width=-25), collapse="\n")
    cat(formatC("Target variables:",width=-25), vrs$yvrs.L2, collapse="\n")
    cat(formatC("Fixed effect predictors:",width=-25), vrs$pvrs.L2, collapse="\n")
  }

  cat("\nPerformed", sprintf("%.0f",itr$burn), "burn-in iterations, and generated", sprintf("%.0f",itr$m),
      "imputed data sets,\neach", sprintf("%.0f",itr$iter), "iterations apart.",
      if(ngr>1){c("\nImputations were carried out seperately within", sprintf("%.0f",ngr), "groups.")},"\n")

  # print convergence diagnostics
  if(!is.null(conv)){

    for(cc in attr(conv,"stats")){

      # summary for Rhat and SDprop
      if(cc=="Rhat"|cc=="SDprop"){

        cout <- matrix(c( sapply(conv, function(z) min(z[,cc])),
                  sapply(conv, function(z) quantile(z[,cc],.25)),
                  sapply(conv, function(z) mean(z[,cc])),
                  sapply(conv, function(z) median(z[,cc])),
                  sapply(conv, function(z) quantile(z[,cc],.75)),
                  sapply(conv, function(z) max(z[,cc])) ), ncol=6 )
        rownames(cout) <- c("Beta:",if(isL2) "Beta2:","Psi:","Sigma:")
        colnames(cout) <- c("Min","25%","Mean","Median","75%","Max")
        clab <- switch(cc, Rhat="\nPotential scale reduction (Rhat, imputation phase):\n",
                           SDprop="\nGoodness of approximation (imputation phase):\n")
        cat(clab,"\n")
        print.table(round(cout,3))

        clab <- switch(cc, Rhat="\nLargest potential scale reduction:\n",
                           SDprop="\nPoorest approximation:\n")
        cat(clab)
        maxval <- lapply(conv, function(a) a[which.max(a[,cc]),1:2])
        cat("Beta: [", paste(maxval$beta,collapse=",") ,"], ",
            if(isL2) paste0("Beta2: [", paste(maxval$beta2,collapse=",") ,"], "),
            "Psi: [", paste(maxval$psi,collapse=",") ,"], ",
            "Sigma: [", paste(maxval$sigma,collapse=",") ,"]\n", sep="")

      }

      # summary for ACF
      if(cc=="ACF"){

        cout <- c( sapply(conv, function(z) mean(z[,"lag-1"])),
                   sapply(conv, function(z) mean(z[,"lag-k"])),
                   sapply(conv, function(z) mean(z[,"lag-2k"])),
                   sapply(conv, function(z) max(z[,"lag-1"])),
                   sapply(conv, function(z) max(z[,"lag-k"])),
                   sapply(conv, function(z) max(z[,"lag-2k"])) )
        neg <- cout<0
        cout <- sprintf(cout,fmt="%.3f")
        cout[neg] <- gsub("^-0","-",cout[neg])
        cout[!neg] <- gsub("^0"," ",cout[!neg])
        cout <- matrix(cout, 3+isL2, 6)
        cout <- rbind(c(" Lag1"," Lagk","Lag2k"," Lag1"," Lagk","Lag2k"), cout)
        rownames(cout) <- c("","Beta:",if(isL2) "Beta2:","Psi:","Sigma:")
        colnames(cout) <- c(" Mean","",""," Max","","")
        cat("\nAutocorrelation (ACF, imputation phase):\n\n")
        print.table(cout)

        cat("\nLargest autocorrelation at lag k:\n")
        maxval <- lapply(conv, function(a) a[which.max(abs(a[,"lag-k"])),1:2])
        cat("Beta: [", paste(maxval$beta,collapse=",") ,"], ",
            if(isL2) paste0("Beta2: [", paste(maxval$beta2,collapse=",") ,"], "),
            "Psi: [", paste(maxval$psi,collapse=",") ,"], ",
            "Sigma: [", paste(maxval$sigma,collapse=",") ,"]\n", sep="")

      }
    }
  }

  # missing data rates
  mdrout <- t(as.matrix(mdr))
  rownames(mdrout) <- "MD%"
  cat("\nMissing data per variable:\n")
  print.table(mdrout)

  cat("\n")

  invisible()
}

