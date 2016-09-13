\name{anova.mitml.result}
\alias{anova.mitml.result}
\title{Compare several nested models}
\description{
Performs model comparisons for a series of nested statistical models fitted using \code{with.mitml.list}.
}

\usage{

\method{anova}{mitml.result}(object, ...)

}

\arguments{

  \item{object}{An object of class \code{mitml.result} as produced by \code{with.mitml.list}.}
  \item{\dots}{Additional objects of class \code{mitml.result} to be included in the comparison.}

}

\details{

  This function performs several model comparisons between models fitted using \code{with.mitml.list}.
  If possible, the models are compared using the \eqn{D_3} statistic (Meng & Rubin, 1992).
  If this method is unavailable, the \eqn{D_2} statistic is used instead (Li, Meng, Raghunathan, & Rubin, 1991).
  The \eqn{D_3} method currently supports linear models and linear mixed-effects models with a single cluster variable as estimated by \code{lme4} or \code{nlme} (see Laird, Lange, & Stram, 1987).

  This function is essentially a wrapper for \code{\link{testModels}} with the advantage that several models can be compared simultaneously.
  All model comparisons are likelihood-based.
  For further options for model comparisons (e.g., Wald-based procedures) and finer control, see \code{testModels}.

}

\value{

  Returns a list containing the results of each model comparison.
  A \code{print} method is used for better readable console output. 

}

\author{Simon Grund}
\seealso{\code{\link{with.mitml.list}}, \code{\link{testModels}}}

\examples{
require(lme4)
data(studentratings)

fml <- ReadDis + SES ~ ReadAchiev + (1|ID)
imp <- panImpute(studentratings, formula=fml, n.burn=1000, n.iter=100, m=5)

implist <- mitmlComplete(imp, print=1:5)

# simple comparison (same as testModels)
fit0 <- with(implist, lmer(ReadAchiev ~ (1|ID), REML=FALSE))
fit1 <- with(implist, lmer(ReadAchiev ~ ReadDis + (1|ID), REML=FALSE))
anova(fit1,fit0)

\dontrun{
# multiple comparisons
fit2 <- with(implist, lmer(ReadAchiev ~ ReadDis + (1+ReadDis|ID), REML=FALSE))
anova(fit2,fit1,fit0)
}
}

\references{
Meng, X.-L., & Rubin, D. B. (1992). Performing likelihood ratio tests with multiply-imputed data sets. \emph{Biometrika, 79}, 103-111. 

Laird, N., Lange, N., & Stram, D. (1987). Maximum likelihood computations with repeated measures: Application of the em algorithm. \emph{Journal of the American Statistical Association, 82}, 97-105.

Li, K. H., Raghunathan, T. E., & Rubin, D. B. (1991). Large-sample significance levels from multiply imputed data using moment-based statistics and an F reference distribution. \emph{Journal of the American Statistical Association, 86}, 1065-1073.
}
\keyword{methods}
