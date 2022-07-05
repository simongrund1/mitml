\name{testEstimates}
\alias{testEstimates}
\title{Compute final estimates and inferences}
\description{
Computes final parameter estimates and inferences from multiply imputed data sets.
}

\usage{

testEstimates(model, qhat, uhat, extra.pars = FALSE, df.com = NULL, ...)

}
\arguments{

  \item{model}{A list of fitted statistical models as produced by \code{\link{with.mitml.list}} or similar.}
  \item{qhat}{A matrix or list containing the point estimates of the parameters for each imputed data set (see 'Details').}
  \item{uhat}{(optional) An array, matrix, or list containing the variance estimates (i.e., squared standard errors) of the parameters for each imputed data set (see 'Details').}
  \item{extra.pars}{A logical flag indicating if estimates of additional parameters (e.g., variance components) should be calculated. Default is \code{FALSE}.}
  \item{df.com}{(optional) A numeric vector denoting the complete-data degrees of freedom for the hypothesis tests (see 'Details').}
  \item{\dots}{Not used.}

}
\details{

This function calculates pooled parameter estimates and inferences as suggested by Rubin (1987, "Rubin's rules") for each parameter of the fitted model.
The parameters can either be extracted automatically from the fitted statistical models (\code{model}) or provided manually as matrices, arrays, or lists (\code{qhat} and \code{uhat}, see 'Examples').

Rubin's original method assumes that the complete-data degrees of freedom are infinite, which is reasonable in larger samples.
Alternatively, the degrees of freedom can be adjusted for smaller samples by specifying \code{df.com} (Barnard & Rubin, 1999).
The \code{df.com} argument can either be a single number if the degrees of freedom are equal for all parameters being tested, or a numeric vector with one element per parameter.

Using the \code{extra.pars} argument, pooled estimates for additional parameters can be requested (e.g., variance components).
This option is available for a number of models but may not provide estimates for all parameters in all model types.
In such a case, users may extract the estimates of additional parameters by hand and pool them with the \code{qhat} argument (see 'Examples').
No inferences are calculated for pooled additional parameters.

Currently, the procedure supports automatic extraction of model parameters from models that define \code{coef} and \code{vcov} methods (e.g., \code{lm}, \code{glm}, \code{lavaan} and others) as well as multilevel models estimated with \code{lme4} or \code{nlme} and GEEs estimated with \code{geepack}.
The arguments \code{qhat} and \code{uhat} provide a general method for pooling parameter estimates regardless of model type (see 'Examples').
Support for further models may be added in future releases.

}
\value{

A list containing the pooled parameter and inferences.
A \code{print} method is used for more readable output. 

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
imp <- panImpute(studentratings, formula = fml, n.burn = 1000, n.iter = 100, m = 5)

implist <- mitmlComplete(imp)

# fit multilevel model using lme4
require(lme4)
fit.lmer <- with(implist, lmer(SES ~ (1|ID)))

# * Example 1: pool estimates of fitted models (automatic)
# pooled estimates and inferences separately for each parameter (Rubin's rules)
testEstimates(fit.lmer)

# ... adjusted df for finite samples
testEstimates(fit.lmer, df.com = 49)

# ... with additional table for variance components and ICCs
testEstimates(fit.lmer, extra.pars = TRUE)

# * Example 2: pool estimates using matrices or lists (qhat, uhat)
fit.lmer <- with(implist, lmer(SES ~ ReadAchiev + (1|ID)))

qhat <- sapply(fit.lmer, fixef)
uhat <- sapply(fit.lmer, function(x) diag(vcov(x)))

testEstimates(qhat = qhat)
testEstimates(qhat = qhat, uhat = uhat)
}
