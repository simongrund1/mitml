\name{multilevelR2}
\alias{multilevelR2}
\title{Calculate R-squared measures for multilevel models}
\description{
Calculates several measures for the proportion of explained variance in a fitted linear mixed-effects (i.e., ) multilevel model (or a list of fitted models).
}

\usage{

multilevelR2(model, print = c("RB1", "RB2", "SB", "MVP"))

}

\arguments{

  \item{model}{Either a fitted linear mixed-effects model as produced by \code{lme4} or \code{nlme}, or a list of fitted models as produced by \code{with.mitml.list}.}
  \item{print}{A character vector denoting which measures should be calculated (see details). Default is to printing all measures.}

}

\details{

  This function calculates several measures of explained variance (\eqn{R^2}) for linear-mixed effects models.
  It can be used with a single model, as produced by the packages \code{lme4} or \code{nlme}, or a list of fitted models produced by \code{with.mitml.list}.
  In the latter case, the \eqn{R^2} measures are calculated separately for each imputed data set and then averaged across data sets.

  Different \eqn{R^2} measures can be requested using the \code{print} argument.
  Specifying \code{RB1} and \code{RB2} will return the explained variance at level 1 and level 2, respectively, according to Raudenbush and Bryk (2002, pp. 74 and 79).
  Specifying \code{SB} will return the total variance explained according to Snijders and Bosker (2012, p. 112).
  Specifying \code{MVP} will return the total variance explained based on ``multilevel variance partitioning'' as proposed by LaHuis, Hartman, Hakoyama, and Clark (2014).

}

\value{

Returns a numeric vector containing the \eqn{R^2} measures requested in \code{print}.

}

\note{

  Calculating \eqn{R^2} measures is currently only supported for two-level models with a single cluster variable.

}

\author{Simon Grund}
\references{
LaHuis, D. M., Hartman, M. J., Hakoyama, S., & Clark, P. C. (2014). Explained variance measures for multilevel models. \emph{Organizational Research Methods}, 17, 433-451.

Raudenbush, S. W., & Bryk, A. S. (2002). Hierarchical linear models: Applications and data analysis methods (2nd ed.). Thousand Oaks, CA: Sage.

Snijders, T. A. B., & Bosker, R. J. (2012). Multilevel analysis: An introduction to basic and advanced multilevel modeling. Thousand Oaks, CA: Sage.
}
\examples{
require(lme4)
data(studentratings)

fml <- MathAchiev + ReadAchiev + CognAbility ~ 1 + (1|ID)
imp <- panImpute(studentratings, formula = fml, n.burn = 1000, n.iter = 100, m = 5)

implist <- mitmlComplete(imp)

fit <- with(implist, lmer(MathAchiev ~ 1 + CognAbility + (1|ID)))
multilevelR2(fit)
}
