\name{testConstraints}
\alias{testConstraints}
\title{Test functions and constraints of model parameters}
\description{
Performs hypothesis tests for arbitrary functions of the model parameters using the Delta method.
}

\usage{

testConstraints(model, qhat, uhat, constraints, method = c("D1", "D2"),
  ariv = c("default", "positive"), df.com = NULL)

}
\arguments{

  \item{model}{A list of fitted statistical models as produced by \code{\link{with.mitml.list}} or similar.}
  \item{qhat}{A matrix or list containing the point estimates of the parameters for each imputed data set (see 'Details').}
  \item{uhat}{An array or list containing the variance-covariance matrix of the parameters for each imputed data set (see 'Details').}
  \item{constraints}{A character vector specifying constraints or functions of the vector of model parameters to be tested (see 'Details').}
  \item{method}{A character string denoting the method by which the test is performed. Can be \code{"D1"} or \code{"D2"} (see 'Details'). Default is \code{"D1"}.}
  \item{ariv}{A character string denoting how the ARIV is calculated. Can be \code{"default"} or \code{"positive"} (see 'Details').}
  \item{df.com}{(optional) A single number or a numeric vector denoting the complete-data degrees of freedom for the hypothesis test (see 'Details'). Only used if \code{method = "D1"}.}

}
\details{

This function performs tests of arbitrary functions (or constraints) of the model parameters using similar methods as \code{\link{testModels}}.
The function relies on the Delta method (e.g., Casella & Berger, 2002) for testing functions of the parameters and assumes that their sampling distribution is approximately normal.
The parameters can either be extracted automatically from the fitted statistical models (\code{model}) or provided manually as matrices, arrays, or lists (\code{qhat} and \code{uhat}, see 'Examples').

Constraints and other functions of the model parameters are specified in the \code{constraints} argument.
The constraints must be supplied as a character vector, where each string denotes a function or a constraint to be tested (see 'Examples').

The Wald-like tests that are carried out by \code{testConstraints} are pooled across the imputed data sets with the \eqn{D_1} (Li, Raghunathan & Rubin, 1991) or \eqn{D_2} (Li, Meng, Raghunathan & Rubin, 1991) method, where \eqn{D_1} operates on the constrained point and variance estimates, and \eqn{D_2} operates on the Wald-statistics (for additional details, see \code{testModels}).
The pooled estimates and standard errors reported in the output are always based on \eqn{D_1}.

For \eqn{D_1}, the complete-data degrees of freedom can be adjusted for smaller samples by specifying \code{df.com} (see \code{testModels}).

This function supports general statistical models that define \code{coef} and \code{vcov} methods (e.g., \code{lm}, \code{glm}, \code{lavaan} and others) as well as multilevel models estimated with \code{lme4} or \code{nlme} and GEEs estimated with \code{geepack}.
The arguments \code{qhat} and \code{uhat} provide a general method for pooling parameter estimates regardless of model type (see 'Examples').
Support for further models may be added in future releases.

The \code{ariv} argument determines how the average relative increase in variance (ARIV) is calculated (see \code{testModels}).
If \code{ariv = "default"}, the default estimators are used.
If \code{ariv = "positive"}, the default estimators are used but constrained to take on strictly positive values.

}
\value{

A list containing the results of the model comparison.
A \code{print} method is used for more readable output. 

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
imp <- panImpute(studentratings, formula = fml, n.burn = 1000, n.iter = 100, m = 5)

implist <- mitmlComplete(imp)

# fit simple regression model
fit.lm <- with(implist, lm(SchClimate ~ ReadDis + MathDis))

# apply Rubin's rules
testEstimates(fit.lm)

# * Example 1: test 'identity' function of two parameters (automatic)
# test equivalent to model comparison with a restricted model (without 'ReadDis'
# and 'MathDis')

cons <- c("ReadDis", "MathDis")
testConstraints(fit.lm, constraints = cons)

# ... adjusting for finite samples
testConstraints(fit.lm, constraints = cons, df.com = 749)

# ... using D2
testConstraints(fit.lm, constraints = cons, method = "D2")

# * Example 2: test for equality of two parameters
# tests the hypothesis that the coefficients pertaining to 'ReadDis' and 'MathDis'
# are equal (ReadDis = MathDis)

cons <- c("ReadDis-MathDis")
testConstraints(fit.lm, constraints = cons)

# * Example 3: test against a fixed value
# tests the hypothesis that the coefficient of "ReadDis" is equal to one
# (i.e., 'ReadDis' - 1 == 0)

cons <- c("ReadDis-1")
testConstraints(fit.lm, constraints = cons)

# * Example 4: test 'identity' function of two parameters (qhat, uhat)

fit.lm <- with(implist, lm(SchClimate ~ ReadDis + MathDis))

qhat <- sapply(fit.lm, coef)
uhat <- sapply(fit.lm, function(x) vcov(x), simplify = "array")

cons <- c("ReadDis", "MathDis")
testConstraints(qhat = qhat, uhat = uhat, constraints = cons)
}
