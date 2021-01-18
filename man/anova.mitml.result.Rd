\name{anova.mitml.result}
\alias{anova.mitml.result}
\title{Compare several nested models}
\description{
Performs model comparisons for a series of nested statistical models fitted using \code{with.mitml.list}.
}

\usage{

\method{anova}{mitml.result}(object, ..., method=c("D3","D4","D2"), ariv=c("default","positive","robust"), data=NULL)

}

\arguments{

  \item{object}{An object of class \code{mitml.result} as produced by \code{with.mitml.list}.}
  \item{\dots}{Additional objects of class \code{mitml.result} to be included in the comparison.}
  \item{method}{A character string denoting the method used for the model comparisons. Can be either \code{"D3"}, \code{"D4"} or \code{"D2"} (see details). Default is \code{"D3"}.}
  \item{ariv}{A character string denoting how the ARIV is calculated. Can be either \code{"default"}, \code{"positive"}, or \code{"robust"} (see details).}
  \item{data}{(optional) A list of imputed data sets (see details).}

}

\details{

  This function performs likelihood-based comparisons between multiple statistical models fitted with \code{with.mitml.list}.
  If possible, the models are compared using the \eqn{D_3} statistic (Meng & Rubin, 1992).
  If this method is unavailable, the \eqn{D_4} or \eqn{D_2} statistic is used instead (Chan & Meng, 2019; Li, Meng, Raghunathan, & Rubin, 1991).

  This function is essentially a wrapper for \code{\link{testModels}} with the advantage that several models can be compared simultaneously.
  For a list of supported models and further options for more specific model comparisons, see \code{\link{testModels}}.

  The \code{ariv} argument affects how the average relative increase in variance is calculated (see also \code{\link{testModels}}).
  Note that the \eqn{D_4} method can fail if the data to which the model was fitted cannot be found.
  In such a case, the \code{data} argument can be used to specify the list of imputed data sets directly (see also \code{\link{testModels}}).

}

\value{

  Returns a list containing the results of each model comparison.
  A \code{print} method is used for better readable output.

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
