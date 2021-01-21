# ***
# misc. methods
#

# * check for REML fit

.checkREML <- function(object, ...) UseMethod(".checkREML", object)
.checkREML.default <- function(object) return(FALSE)
.checkREML.merMod <- function(object) return(lme4::isREML(object))
.checkREML.lme <- function(object) return(object$method == "REML")

# * update REML fit with ML

.updateML <- function(object, ...) UseMethod(".updateML", object)
.updateML.default <- function(object) return(object)
.updateML.merMod <- function(object) return(.localUpdate(object, REML = FALSE))
.updateML.lme <- function(object) return(.localUpdate(object, data = object$data, method = "ML"))

# * determine degrees of freedom

.getDFs <- function(object, ...) UseMethod(".getDFs", object)

.getDFs.default <- function(object){

  df <- NULL

  # try logLik
  df.try <- try(attr(logLik(object), "df"), silent = TRUE)
  if(!inherits(df.try, "try-error")){
    df <- df.try
    attr(df, "type") <- "logLik"
  }

  # try df.residual and sample size (nobs, residuals)
  # NOTE: does not account for scale parameters (e.g., residual variance)
  if(is.null(df)){
    rdf <- try(df.residual(object), silent = TRUE)
    n <- try(nobs(object), silent = TRUE)
    if(inherits(n, "try-error")) n <- try(length(predict(object)), silent = TRUE)
    if(inherits(n, "try-error")) n <- try(length(residuals(object)), silent = TRUE)
    if(!inherits(rdf, "try-error") && !inherits(n, "try-error")){
      df <- n - rdf
      attr(df, "type") <- "df.residual"
    }
  }

  return(df)

}

.getDFs.lavaan <- function(object){

  df <- attr(lavaan::logLik(object), "df")
  attr(df, "type") <- "logLik"
  return(df)

}


# * extract model formula

.getFormula <- function(object, ...) UseMethod(".getFormula", object)

.getFormula.default <- function(object){

    fml <- try(deparse(formula(object)))
    if(inherits(fml, "try-error")) fml <- NULL
    fml <- Reduce(paste, fml)

    return(fml)

}

.getFormula.lme <- function(object){

  fe.fml <- deparse(formula(object))
  re.fml <- lapply(formula(object$modelStruct$reStruct), deparse)
  for(ff in seq_along(re.fml)) re.fml[[ff]] <- paste0(re.fml[[ff]], "|", names(re.fml)[ff])
  fml <- paste(c(fe.fml, unlist(re.fml)), collapse = ", ")

  return(fml)

}

.getFormula.lavaan <- function(object){

  cll <- getCall(object)
  fml <- deparse(cll[c(1, match("model", names(cll)))])
  fml <- sub(")$", ", ...)", fml)

  return(fml)

}
