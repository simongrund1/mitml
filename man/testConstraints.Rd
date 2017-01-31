\name{testConstraints}
\alias{testConstraints}
\title{Test functions and constraints of model parameters}
\description{
Performs hypothesis tests for arbitrary functions of a parameter vector using the Delta method.
}

\usage{

testConstraints(model, qhat, uhat, constraints, method=c("D1","D2"), df.com=NULL)

}
\arguments{

  \item{model}{A list of fitted statistical models (see examples).}
  \item{qhat, uhat}{Two matrices/arrays or lists containing estimates of the parameter vector and its covariance matrix, respectively, for each imputed data set (see examples).}
  \item{constraints}{A character vector specifying constraints or functions of the vector of model parameters to be tested.}
  \item{method}{A character string denoting the method by which the test is performed. Can be either \code{"D1"} or \code{"D2"} (see details). Default is to \code{"D1"}.}
  \item{df.com}{(optional) A single number or a numeric vector denoting the complete-data degrees of freedom for the hypothesis test. Only used if \code{method="D1"}.}

}
\details{

This function is similar in functionality to \code{\link{testModels}} but extended to arbitrary functions (or constraints) of the model parameters.
The function is based on the Delta method (e.g., Casella & Berger, 2002) according to which any function of the parameters can be tested using Wald-like methods, assuming that their sampling distribution is approximately normal.
It is assumed that the parameters can be extracted using \code{coef} and \code{vcov} methods from the fitted models (or similar; e.g., regression coefficients, fixed effects in multilevel models)
In cases where this is not possible, hypothesis tests can be carried out using user-supplied matrices/arrays or lists (\code{qhat} and \code{uhat}, see examples).

Constraints and functions of the model parameters can be specified in the \code{constraints} argument.
The constraints must be supplied as a character vector, where each string denotes a function or a constraint to be tested (see examples).

The Wald-like tests that are carried out by \code{testConstraints} can be aggregated across data sets by methods \eqn{D_1} (Li, Raghunathan & Rubin, 1991) and \eqn{D_2} (Li, Meng, Raghunathan & Rubin, 1991), where \eqn{D_1} operates on the constrained estimates and standard errors, and \eqn{D_2} operates on the Wald-statistics (for an explanation, see \code{testModels}).

For \eqn{D_1}, the complete-data degrees of freedom can be adjusted for smaller samples by specifying \code{df.com}.

Currently, the procedure supports statistical models that define \code{coef} and \code{vcov} methods (e.g., \code{lm}), multilevel models estimated with \code{lme4} or \code{nlme}, and GEEs estimated with \code{geepack}.
The arguments \code{qhat} and \code{uhat} allow for more general hypothesis tests regardless of model class.
Support for further models may be added in future releases.

}
\value{

Returns a list containing the results of the model comparison, the constrained estimates and standard errors, and the relative increase in variance due to nonresponse (Rubin, 1987).
A \code{print} method is used for better readable console output. 

}

\references{
Casella, G., & Berger, R. L. (2002). \emph{Statistical inference (2nd. Ed.)}. Pacific Grove, CA: Duxbury.

Li, K.-H., Meng, X.-L., Raghunathan, T. E., & Rubin, D. B. (1991). Significance levels from repeated p-values with multiply-imputed data. \emph{Statistica Sinica, 1}, 65-92.

Li, K. H., Raghunathan, T. E., & Rubin, D. B. (1991). Large-sample significance levels from multiply imputed data using moment-based statistics and an F reference distribution. \emph{Journal of the American Statistical Association, 86}, 1065-1073.
}

\author{Simon Grund}

\seealso{\code{\link{testModels}}, \code{\link{with.mitml.list}}}
\examples{
data(studentratings)

fml <- MathDis + ReadDis + SchClimate ~ (1|ID)
imp <- panImpute(studentratings, formula=fml, n.burn=1000, n.iter=100, m=5)

implist <- mitmlComplete(imp, print=1:5)

# fit simple regression model
fit.lm <- with(implist, lm(SchClimate ~ ReadDis + MathDis))

# apply Rubin's rules
testEstimates(fit.lm)

# * Example 1: test 'identity' function of two parameters
# multi-parameter hypothesis test, equivalent to model comparison

cons <- c("ReadDis","MathDis")
testConstraints(fit.lm, constraints=cons)

# ... adjusting for finite samples
testConstraints(fit.lm, constraints=cons, df.com=749)

# ... using D2
testConstraints(fit.lm, constraints=cons, method="D2")

# * Example 2: test for equality of two effects
# tests the hypothesis that the effects of 'ReadDis' and 'MathDis'
# are equal (ReadDis=MathDis)

cons <- c("ReadDis-MathDis")
testConstraints(fit.lm, constraints=cons)

# * Example 3: test against a fixed value
# tests the hypothesis that the effect of "ReadDis" is one (ReadDis=1)

cons <- c("ReadDis-1")
testConstraints(fit.lm, constraints=cons)

# * Example 4: test 'identity' using arrays and list

fit.lm <- with(implist, lm(SchClimate ~ ReadDis + MathDis))

cons <- c("ReadDis","MathDis")
qhat <- sapply(fit.lm, coef)
uhat <- sapply(fit.lm, function(x) vcov(x), simplify="array")
testConstraints(qhat=qhat, uhat=uhat, constraints=cons)
}
