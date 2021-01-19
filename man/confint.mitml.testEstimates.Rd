\name{confint.mitml.testEstimates}
\alias{confint.mitml.testEstimates}
\title{Compute confidence intervals}
\description{
Computes confidence intervals on the basis of the final parameter estimates and inferences given by \code{\link{testEstimates}}.
}

\usage{

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

This function computes confidence intervals with the given confidence level for the pooled parameters on the basis of a \eqn{t}-distribution, with estimates, standard errors, and degrees of freedom as returned by \code{\link{testEstimates}}.

}
\value{

A matrix containing the lower and upper bounds of the confidence intervals.

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

# compute confidence intervals
confint(est)

# ... with different confidence levels
confint(est, level = 0.90)
confint(est, level = 0.999)
}

