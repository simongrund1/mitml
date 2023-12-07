.pool.estimates <- function(Qhat, Uhat, m, diagonal = FALSE, df.com = NULL, nms = NULL) {

  # pool point estimates
  Qbar <- apply(Qhat, 1, mean)

  # pool variances and inferences
  if (!is.null(Uhat)) {

    Ubar <- apply(Uhat, 1:2, mean)
    B <- tcrossprod(Qhat - Qbar) / (m-1)
    T <- Ubar + (1 + m^(-1)) * B

    se <- sqrt(diag(T))
    t <- Qbar/se

    r <- (1 + m^(-1)) * diag(B) / diag(Ubar)

    # compute degrees of freedom
    v <- vm <- (m-1) * (1 + r^(-1))^2
    if (!is.null(df.com)) {
      lam <- r / (r+1)
      vobs <- (1-lam) * ((df.com+1) / (df.com+3)) * df.com
      v <- (vm^(-1) + vobs^(-1))^(-1)
    }

    fmi <- (r + 2 / (v+3)) / (r+1)
    pval <- 2 * (1 - pt(abs(t), df = v))

    # create output
    out <- matrix(c(Qbar, se, t, v, pval, r, fmi), ncol = 7)
    colnames(out) <- c("Estimate", "Std.Error", "t.value", "df", "P(>|t|)", "RIV", "FMI")
    if(!diagonal) attr(out, "T") <- T

  } else {

    # create output
    out <- matrix(Qbar, ncol = 1)
    colnames(out) <- "Estimate"

  }

  # parameter names
  rownames(out) <- nms
  attr(out, "par.labels") <- attr(nms, "par.labels")

  return(out)

}


.D1 <- function(Qhat, Uhat, df.com) {
# pooling for multidimensional estimands (D1, Li et al., 1991; Reiter, 2007)

  k <- dim(Qhat)[1]
  m <- dim(Qhat)[2]

  # D1
  Qbar <- apply(Qhat, 1, mean)
  Ubar <- apply(Uhat, c(1, 2), mean)

  B <- cov(t(Qhat))
  r <- (1+m^(-1)) * sum(diag(B %*% solve(Ubar))) / k
  Ttilde <- (1 + r) * Ubar

  val <- t(Qbar) %*% solve(Ttilde) %*% Qbar / k

  # compute degrees of freedom (df2)
  t <- k*(m-1)
  if(!is.null(df.com)) {

    # warn about poor behavior for t<=4
    if (t <= 4) {
      warning("Degrees of freedom (df.com) may not be trustworthy, because the number of imputations is too low (m \u2264 5). To obtain trustworthy results, re-run the procedure with a larger number of imputations.")
    }

    # small-sample degrees of freedom (Reiter, 2007; Eq. 1-2)
    a <- r * t / (t-2)
    vstar <- ( (df.com+1) / (df.com+3) ) * df.com

    c0 <- 1 / (t-4)
    c1 <- vstar - 2 * (1+a)
    c2 <- vstar - 4 * (1+a)

    z <- 1 / c2 +
         c0 * (a^2 * c1 / ((1+a)^2 * c2)) +
         c0 * (8*a^2 * c1 / ((1+a) * c2^2) + 4*a^2 / ((1+a) * c2)) +
         c0 * (4*a^2 / (c2 * c1) + 16*a^2 * c1 / c2^3) +
         c0 * (8*a^2 / c2^2)
 
    v <- 4 + 1/z

  } else {

    if (t > 4) {
      v <- 4 + (t-4) * (1 + (1 - 2*t^(-1)) * (r^(-1)))^2
    } else {
      v <- t * (1 + k^(-1)) * ((1 + r^(-1))^2) / 2
    }

  }

  return(list(F = val, k = k, v = v, r = r))

}

.D2 <- function(d, k) {
# pooling for multidimensional estimands (D2, Li, Meng et al., 1991)

  m <- length(d)

  # D2
  dbar <- mean(d)

  r <- (1 + m^(-1)) * var(sqrt(d))

  val <- (dbar/k - (m+1)/(m-1) * r) / (1+r)

  # compute degrees of freedom (df2)
  v <- k^(-3/m) * (m-1) * (1 + r^(-1))^2

  return(list(F = val, k = k, v = v, r = r))

}

