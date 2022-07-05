\name{confint.mitml.testEstimates}
\alias{coef.mitml.testEstimates}
\alias{vcov.mitml.testEstimates}
\alias{confint.mitml.testEstimates}
\title{Summarize and extract pooled parameter estimates}
\description{
Functions for extracting results and computing confidence intervals from the pooled parameter estimates computed with \code{\link{testEstimates}}.
}

\usage{

\method{coef}{mitml.testEstimates}(object, ...)
\method{vcov}{mitml.testEstimates}(object, ...)
\method{confint}{mitml.testEstimates}(object, parm, level = 0.95, ...)

}
\arguments{

  \item{object}{An object of class \code{mitml.testEstimates} as produced by \code{testEstimates}.}
  \item{parm}{(optional) A reference to the parameters for which to calculate confidence intervals. Can be a character or integer vector denoting names or position of parameters, respectively.
              If missing, all parameters are considered (the default).}
  \item{level}{The confidence level. Default is to \code{0.95} (i.e., 95\%).}
  \item{\dots}{Not being used.}

}
\details{

These functions can be used to extract information and compute additional results from pooled parameter estimates.
The \code{coef} and \code{vcov} methods extract the pooled parameter estimates and their pooled variance-covariance matrix (with the squared standard errors in the diagonal).
The \code{confint} method computes confidence intervals with the given confidence level for the pooled parameters on the basis of a \eqn{t}-distribution, with estimates, standard errors, and degrees of freedom as returned by \code{\link{testEstimates}}.

}
\value{


For \code{coef}: A vector containing the pooled parameter estimates
For \code{vcov}: A matrix containing the pooled variance-covariance matrix of the parameter estimates.
For \code{confint}: A matrix containing the lower and upper bounds of the confidence intervals.

}

\author{Simon Grund}

\seealso{\code{\link{testEstimates}}}
\examples{
data(studentratings)

fml <- ReadDis ~ ReadAchiev + (1|ID)
imp <- panImpute(studentratings, formula = fml, n.burn = 500, n.iter = 100, m = 5)

implist <- mitmlComplete(imp)

# fit regression model
fit <- with(implist, lm(ReadDis ~ 1 + ReadAchiev))
est <- testEstimates(fit)

# extract parameter estimates and variance-covariance matrix
coef(est)
vcov(est)

# compute confidence intervals
confint(est)

# ... with different confidence levels
confint(est, level = 0.90)
confint(est, level = 0.999)
}

