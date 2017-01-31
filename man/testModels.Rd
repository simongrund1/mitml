\name{testModels}
\alias{testModels}
\title{Test multiple parameters and compare nested models}
\description{
Performs multi-parameter hypothesis tests for a vector of statistical parameters and compares nested statistical models obtained from multiply imputed data sets.
}

\usage{

testModels(model, null.model, method=c("D1","D2","D3"), 
  use=c("wald","likelihood"), df.com=NULL)

}
\arguments{

  \item{model}{A list of fitted statistical models (see examples).}
  \item{null.model}{A list of fitted (more restrictive) statistical models.}
  \item{method}{A character string denoting the method by which the test is performed. Can be either \code{"D1"}, \code{"D2"} or \code{"D3"} (see details). Default is \code{"D1"}.}
  \item{use}{A character string denoting Wald- or likelihood-based based tests. Can be either \code{"wald"} or \code{"likelihood"}. Only used if \code{method="D2"}.}
  \item{df.com}{(optional) A single number or a numeric vector denoting the complete-data degrees of freedom for the hypothesis test. Only used if \code{method="D1"}.}

}
\details{

This function compares two nested statistical models which differ by one or more parameters.
In other words, the function performs Wald-like and likelihood-ratio hypothesis tests for the statistical parameters by which the two models differ.

The general approach to Wald-like inference for multi-dimensional estimands was introduced Rubin (1987) and further developed by Li, Raghunathan and Rubin (1991).
This procedure is commonly referred to as \eqn{D_1} and can be used by setting \code{method="D1"}.
\eqn{D_1} is the multi-parameter equivalent of \code{\link{testEstimates}}, that is, it tests multiple parameters simultaneously.
For \eqn{D_1}, the complete-data degrees of freedom are assumed to be infinite, but they can be adjusted for smaller samples by supplying \code{df.com} (Reiter, 2007).

An alternative method for Wald-like hypothesis tests was suggested by Li, Meng, Raghunathan and Rubin (1991).
The procedure is often called \eqn{D_2} and can be used by setting \code{method="D2"}.
\eqn{D_2} calculates the Wald-test directly for each data set and then aggregates the resulting \eqn{\chi^2} values.
The source of these values is specified by the \code{use} argument.
If \code{use="wald"} (the default), then a Wald-like hypothesis test similar to \eqn{D_1} is performed.
If \code{use="likelihood"}, then the two models are compared through their likelihood.

A third method relying on likelihood-based comparisons was suggested by Meng and Rubin (1992).
This procedure is referred to as \eqn{D_3} and can be used by setting \code{method="D3"}.
\eqn{D_3} compares the two models by aggregating the likelihood-ratio test across multiply imputed data sets.

In general, Wald-like hypothesis tests (\eqn{D_1} and \eqn{D_2}) are appropriate if the parameters can be assumed to follow a multivariate normal distribution (e.g., regression coefficients, fixed effects in multilevel models).
Likelihood-based comparisons (\eqn{D_2} and \eqn{D_3}) are also appropriate in such cases and may also be used for variance components.

The function supports different classes of statistical models depending on which \code{method} is chosen.
\eqn{D_1} supports quite general models as long as they define \code{coef} and \code{vcov} methods (or similar) for extracting the parameter estimates and their estimated covariance matrix. 
\eqn{D_2} can be used for the same models (if \code{use="wald"}, or alternatively, for models that define a \code{logLik} method (if \code{use="likelihood"}).
Finally, \eqn{D_3} supports linear models and linear mixed-effects models with a single cluster variable as estimated by \code{lme4} or \code{nlme} (see Laird, Lange, & Stram, 1987).
Support for other statistical models may be added in future releases.

}
\value{

Returns a list containing the results of the model comparison, and the relative increase in variance due to nonresponse (Rubin, 1987).
A \code{print} method is used for better readable console output. 

}

\note{

The \eqn{D_3} method and the likelihood-based \eqn{D_2} assume that models were fit using maximum likelihood (ML).
Models fit using REML are automatically refit using ML.

}

\references{
Meng, X.-L., & Rubin, D. B. (1992). Performing likelihood ratio tests with multiply-imputed data sets. \emph{Biometrika, 79}, 103-111. 

Laird, N., Lange, N., & Stram, D. (1987). Maximum likelihood computations with repeated measures: Application of the em algorithm. \emph{Journal of the American Statistical Association, 82}, 97-105.

Li, K.-H., Meng, X.-L., Raghunathan, T. E., & Rubin, D. B. (1991). Significance levels from repeated p-values with multiply-imputed data. \emph{Statistica Sinica, 1}, 65-92.

Li, K. H., Raghunathan, T. E., & Rubin, D. B. (1991). Large-sample significance levels from multiply imputed data using moment-based statistics and an F reference distribution. \emph{Journal of the American Statistical Association, 86}, 1065-1073.

Reiter, J. P. (2007). Small-sample degrees of freedom for multi-component significance tests with multiple imputation for missing data. \emph{Biometrika, 94}, 502-508.

Rubin, D. B. (1987). \emph{Multiple imputation for nonresponse in surveys}. Hoboken, NJ: Wiley.
}

\author{Simon Grund}

\seealso{\code{\link{anova.mitml.result}}, \code{\link{testEstimates}}, \code{\link{testConstraints}}, \code{\link{with.mitml.list}}}
\examples{
data(studentratings)

fml <- ReadDis + SES ~ ReadAchiev + (1|ID)
imp <- panImpute(studentratings, formula=fml, n.burn=1000, n.iter=100, m=5)

implist <- mitmlComplete(imp, print=1:5)

# * Example 1: multiparameter hypothesis test for 'ReadDis' and 'SES'
# This tests the hypothesis that both effects are zero.

require(lme4)
fit0 <- with(implist, lmer(ReadAchiev ~ (1|ID), REML=FALSE))
fit1 <- with(implist, lmer(ReadAchiev ~ ReadDis + (1|ID), REML=FALSE))

# apply Rubin's rules
testEstimates(fit1)

# multiparameter hypothesis test using D1 (default)
testModels(fit1, fit0)

# ... adjusting for finite samples
testModels(fit1, fit0, df.com=47)

# ... using D2 ("wald", using estimates and covariance-matrix)
testModels(fit1, fit0, method="D2")

# ... using D2 ("likelihood", using likelihood-ratio test)
testModels(fit1, fit0, method="D2", use="likelihood")

# ... using D3 (likelihood-ratio test, requires ML fit)
testModels(fit1, fit0, method="D3")

\dontrun{
# * Example 2: multiparameter test using D3 with nlme

# for D3 to be calculable, the 'data' argument for 'lme' must be
# can be constructed manually

require(nlme)
fit0 <- with(implist, lme(ReadAchiev~1, random=~1|ID,
  data=data.frame(ReadAchiev,ID), method="ML"))
fit1 <- with(implist, lme(ReadAchiev ~ 1 + ReadDis, random=~ 1|ID,
  data=data.frame(ReadAchiev,ReadDis,ID), method="ML"))

# multiparameter hypothesis test using D3
testModels(fit1, fit0, method="D3")
}
}
