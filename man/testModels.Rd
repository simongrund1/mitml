\name{testModels}
\alias{testModels}
\title{Test multiple parameters and compare nested models}
\description{
Performs multi-parameter hypothesis tests for a vector of statistical parameters and compares nested statistical models obtained from multiply imputed data sets.
}

\usage{

testModels(model, null.model, method = c("D1", "D2", "D3", "D4"), 
  use = c("wald", "likelihood"), ariv = c("default", "positive", "robust"),
  df.com = NULL, data = NULL)

}
\arguments{

  \item{model}{A list of fitted statistical models (``full'' model) as produced by \code{\link{with.mitml.list}} or similar.}
  \item{null.model}{A list of fitted statistical models (``restricted'' model) as produced by \code{\link{with.mitml.list}} or similar.}
  \item{method}{A character string denoting the method by which the test is performed. Can be \code{"D1"}, \code{"D2"}, \code{"D3"}, or \code{"D4"} (see 'Details'). Default is \code{"D1"}.}
  \item{use}{A character string denoting Wald- or likelihood-based based tests. Can be either \code{"wald"} or \code{"likelihood"}. Only used if \code{method = "D2"}.}
  \item{ariv}{A character string denoting how the ARIV is calculated. Can be \code{"default"}, \code{"positive"}, or \code{"robust"} (see 'Details').}
  \item{df.com}{(optional) A number denoting the complete-data degrees of freedom for the hypothesis test. Only used if \code{method = "D1"}.}
  \item{data}{(optional) A list of imputed data sets (see 'Details'). Only used if \code{method = "D4"}}

}
\details{

This function compares two nested statistical models fitted to multiply imputed data sets by pooling Wald-like or likelihood-ratio tests.

Pooling methods for Wald-like tests of multiple parameters were introduced by Rubin (1987) and further developed by Li, Raghunathan and Rubin (1991).
The pooled Wald test is referred to as \eqn{D_1} and can be used by setting \code{method = "D1"}.
\eqn{D_1} is the multi-parameter equivalent of \code{\link{testEstimates}}, that is, it tests multiple parameters simultaneously.
For \eqn{D_1}, the complete-data degrees of freedom are assumed to be infinite, but they can be adjusted for smaller samples by supplying \code{df.com} (Reiter, 2007).

An alternative method for Wald-like hypothesis tests was suggested by Li, Meng, Raghunathan and Rubin (1991).
The procedure is called \eqn{D_2} and can be used by setting \code{method = "D2"}.
\eqn{D_2} calculates the Wald-test directly for each data set and then pools the resulting \eqn{\chi^2} values.
The source of these values is specified by the \code{use} argument.
If \code{use = "wald"} (the default), then a Wald test similar to \eqn{D_1} is performed.
If \code{use = "likelihood"}, then the two models are compared with a likelihood-ratio test instead.

Pooling methods for likelihood-ration tests were suggested by Meng and Rubin (1992).
This procedure is referred to as \eqn{D_3} and can be used by setting \code{method = "D3"}.
\eqn{D_3} compares the two models by pooling the likelihood-ratio test across multiply imputed data sets.

Finally, an improved method for pooling likelihood-ratio tests was recommended by Chan & Meng (2019).
This method is referred to as \eqn{D_4} and can be used by setting \code{method = "D4"}.
\eqn{D_4} also compares models by pooling the likelihood-ratio test but does so in a more general and efficient manner.

The function supports different classes of statistical models depending on which \code{method} is chosen.
\eqn{D_1} supports models that define \code{coef} and \code{vcov} methods (or similar) for extracting the parameter estimates and their estimated covariance matrix. 
\eqn{D_2} can be used for the same models (if \code{use = "wald"} and models that define a \code{logLik} method (if \code{use = "likelihood"}).
\eqn{D_3} supports linear models, linear mixed-effects models (see Laird, Lange, & Stram, 1987) with an arbitrary cluster structed if estimated with \code{lme4} or a single cluster if estimated by \code{nlme}, and structural equation models estimated with \code{lavaan} (requires ML estimator, see 'Note').
Finally, \eqn{D_4} supports models that define a \code{logLik} method but can fail if the data to which the model was fitted cannot be found.
In such a case, users can provide the list of imputed data sets directly by specifying the \code{data} argument or refit with the \code{include.data} argument in \code{\link{with.mitml.list}}.
Support for other statistical models may be added in future releases.

The \eqn{D_4}, \eqn{D_3}, and \eqn{D_2} methods support different estimators of the relative increase in variance (ARIV), which can be specified with the \code{ariv} argument.
If \code{ariv = "default"}, the default estimators are used.
If \code{ariv = "positive"}, the default estimators are used but constrained to take on strictly positive values.
This is useful if the estimated ARIV is negative.
If \code{ariv = "robust"}, which is available only for \eqn{D_4}, the "robust" estimator proposed by Chan & Meng (2019) is used.
This method should be used with caution, because it requires much stronger assumptions and may result in liberal inferences if these assumptions are violated.

}
\value{

A list containing the results of the model comparison.
A \code{print} method is used for more readable output. 

}

\note{

The methods \eqn{D_4}, \eqn{D_3}, and the likelihood-based \eqn{D_2} assume that models were fit using maximum likelihood (ML).
Models fit using REML are automatically refit using ML.
Models fit in \code{'lavaan'} using the MLR estimator or similar techniques that require scaled \eqn{chi^2} difference tests are currently not supported.

}

\references{
Chan, K. W., & Meng, X.-L. (2019). Multiple improvements of multiple imputation likelihood ratio tests. ArXiv:1711.08822 [Math, Stat]. \url{http://arxiv.org/abs/1711.08822}

Laird, N., Lange, N., & Stram, D. (1987). Maximum likelihood computations with repeated measures: Application of the em algorithm. \emph{Journal of the American Statistical Association, 82}, 97-105.

Li, K.-H., Meng, X.-L., Raghunathan, T. E., & Rubin, D. B. (1991). Significance levels from repeated p-values with multiply-imputed data. \emph{Statistica Sinica, 1}, 65-92.

Li, K. H., Raghunathan, T. E., & Rubin, D. B. (1991). Large-sample significance levels from multiply imputed data using moment-based statistics and an F reference distribution. \emph{Journal of the American Statistical Association, 86}, 1065-1073.

Meng, X.-L., & Rubin, D. B. (1992). Performing likelihood ratio tests with multiply-imputed data sets. \emph{Biometrika, 79}, 103-111. 

Reiter, J. P. (2007). Small-sample degrees of freedom for multi-component significance tests with multiple imputation for missing data. \emph{Biometrika, 94}, 502-508.

Rubin, D. B. (1987). \emph{Multiple imputation for nonresponse in surveys}. Hoboken, NJ: Wiley.
}

\author{Simon Grund}
\seealso{\code{\link{testEstimates}}, \code{\link{testConstraints}}, \code{\link{with.mitml.list}}, \code{\link{anova.mitml.result}}}
\examples{
data(studentratings)

fml <- ReadDis + SES ~ ReadAchiev + (1|ID)
imp <- panImpute(studentratings, formula = fml, n.burn = 1000, n.iter = 100, m = 5)

implist <- mitmlComplete(imp)

# * Example 1: multiparameter hypothesis test for 'ReadDis' and 'SES'
# This tests the hypothesis that both effects are zero.

require(lme4)
fit0 <- with(implist, lmer(ReadAchiev ~ (1|ID), REML = FALSE))
fit1 <- with(implist, lmer(ReadAchiev ~ ReadDis + (1|ID), REML = FALSE))

# apply Rubin's rules
testEstimates(fit1)

# multiparameter hypothesis test using D1 (default)
testModels(fit1, fit0)

# ... adjusting for finite samples
testModels(fit1, fit0, df.com = 47)

# ... using D2 ("wald", using estimates and covariance-matrix)
testModels(fit1, fit0, method = "D2")

# ... using D2 ("likelihood", using likelihood-ratio test)
testModels(fit1, fit0, method = "D2", use = "likelihood")

# ... using D3 (likelihood-ratio test, requires ML fit)
testModels(fit1, fit0, method = "D3")

# ... using D4 (likelihood-ratio test, requires ML fit)
testModels(fit1, fit0, method = "D4")
}
