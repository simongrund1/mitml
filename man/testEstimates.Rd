\name{testEstimates}
\alias{testEstimates}
\title{Compute final estimates and inferences}
\description{
Computes final parameter estimates and inferences from multiply imputed data sets.
}

\usage{

testEstimates(model, qhat, uhat, var.comp=FALSE, df.com=NULL)

}
\arguments{

  \item{model}{A list of fitted statistical models (see examples).}
  \item{qhat, uhat}{Two matrices or lists containing point and variances estimates, respectively, for each imputed data set (see examples).}
  \item{var.comp}{A logical flag indicating if estimates for variance components should be calculated. Default is to \code{FALSE}.}
  \item{df.com}{(optional) A numeric vector denoting the complete-data degrees of freedom for the hypothesis test.}

}
\details{

This function calculates final parameter estimates and inferences as suggested by Rubin (1987, "Rubin's rules") for each parameter of the fitted model.
In other words, \code{testEstimates} aggregates estimates and standard errors across multiply imputed data sets.
The parameters and standard errors can either be supplied as fitted statistical models (\code{model}), or as two matrices or lists (\code{qhat} and \code{uhat}, see examples).

Rubin's original method assumes that the complete-data degrees of freedom are infinite, which is reasonable in larger samples.
Alternatively, the degrees of freedom can be adjusted for smaller samples by specifying \code{df.com} (Barnard & Rubin, 1999).
The \code{df.com} argument can either be a single number if the degrees of freedom are equal for all tests, or a numeric vector with one element per test.

Using the \code{var.comp} argument, final estimates for variance components and related parameters can be requested.
These will be shown as a separate table within the console output.
Accessing variance components is highly dependent on the model being estimated and not implemented for all models. 
Users may prefer calculating these estimates manually using \code{\link{with.mitml.list}} (see Example 3).
No inferences are calculated for variance components.

Currently, the procedure supports statistical models that define \code{coef} and \code{vcov} methods (e.g., \code{lm}), multilevel models estimated with \code{lme4} or \code{nlme}, and GEEs estimated with \code{geepack}.
The arguments \code{qhat} and \code{uhat} allow for more general aggregation of parameter estimates regardless of model class.
Support for further models may be added in future releases.

}
\value{

Returns a list containing the final parameter and inferences, the relative increase in variance due to nonresponse, and the fraction of missing information (Rubin, 1987).
A \code{print} method is used for better readable console output. 

}

\references{
Barnard, J., & Rubin, D. B. (1999). Small-sample degrees of freedom with multiple imputation. \emph{Biometrika, 86}, 948-955.

Rubin, D. B. (1987). \emph{Multiple imputation for nonresponse in surveys}. Hoboken, NJ: Wiley.
}


\author{Simon Grund}

\seealso{\code{\link{with.mitml.list}}, \code{\link{confint.mitml.testEstimates}}}
\examples{
data(studentratings)

fml <- ReadDis + SES ~ ReadAchiev + (1|ID)
imp <- panImpute(studentratings, formula=fml, n.burn=1000, n.iter=100, m=5)

implist <- mitmlComplete(imp, print=1:5)

# fit multilevel model using lme4
require(lme4)
fit.lmer <- with(implist, lmer(SES ~ (1|ID)))

# * Example 1: combine estimates using model recognition
# final estimates and inferences sperately for each parameter (Rubin's rules)
testEstimates(fit.lmer)

# ... adjusted df for finite samples
testEstimates(fit.lmer, df.com=49)

# ... with additional table for variance components and ICCs
testEstimates(fit.lmer, var.comp=TRUE)

# * Example 2: combine estimates using matrices or lists
fit.lmer <- with(implist, lmer(SES ~ ReadAchiev + (1|ID)))

qhat <- sapply(fit.lmer, fixef)
uhat <- sapply(fit.lmer, function(x) diag(vcov(x)))
testEstimates(qhat=qhat, uhat=uhat)
}
