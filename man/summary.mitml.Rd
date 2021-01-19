\name{summary.mitml}
\alias{summary.mitml}
\title{Summary measures for imputation models}
\description{
Provides summary statistics and additional information on imputations in objects of class \code{mitml}.
}

\usage{

\method{summary}{mitml}(object, n.Rhat = 3, goodness.of.appr = FALSE, autocorrelation = FALSE, ...)

}

\arguments{

  \item{object}{An object of class \code{mitml} as produced by \code{panImpute} or \code{jomoImpute}.}
  \item{n.Rhat}{(optional) An integer denoting the number of segments used for calculating the potential scale reduction factor. Default is \code{3}.}
  \item{goodness.of.appr}{(optional) A logical flag indicating if the goodness of approximation should be printed. Default is \code{FALSE} (see 'Details').}
  \item{autocorrelation}{(optional) A logical flag indicating if the autocorrelation should be printed. Default is \code{FALSE} (see 'Details').}
  \item{\dots}{Not used.}

}

\details{

The \code{summary} method calculates summary statistics for objects of class \code{mitml} as produced by \code{\link{panImpute}} or \code{\link{jomoImpute}}.
The output includes the potential scale reduction factor (PSRF, or \eqn{\hat{R}}) and (optionally) the goodness of approximation and autocorrelation. 

The PSRF is calculated for each parameter of the imputation model and can be used as a convergence diagnostic (Gelman and Rubin, 1992).
Calculation of the PSRFs can be suppressed by setting \code{n.Rhat = NULL}.
The PSRFs are not computed from different chains but by dividing each chain from the imputation phase into a number of segments as denoted by \code{n.Rhat}.
This is slightly different from the original method proposed by Gelman and Rubin.

The goodness of approximation measure indicates what proportion of the posterior standard deviation is due to simulation error.
This is useful for assessing the accuracy of the posterior summaries (e.g., the EAP).
The autocorrelation includes estimates of the autocorrelation in the chains at lag 1 (i.e., for consecutive draws) and for lags \eqn{k} and \eqn{2k}, where \eqn{k} is the number of iterations between imputations.
For lag \eqn{k} and \eqn{2k}, the autocorrelation is slightly smoothed to reduce the influence of noise on the estimates (see \code{\link{plot.mitml}}).

}

\value{

An object of class \code{summary.mitml}.
A print method is used for more readable output.

}

\references{
Gelman, A., and Rubin, D. B. (1992). Inference from iterative simulation using multiple sequences. \emph{Statistical Science, 7}, 457-472.

Hoff, P. D. (2009). \emph{A first course in Bayesian statistical methods}. New York, NY: Springer.
}

\author{Simon Grund}
\seealso{\code{\link{panImpute}}, \code{\link{jomoImpute}}, \code{\link{plot.mitml}}}
\examples{
data(studentratings)

fml <- ReadDis + SES ~ ReadAchiev + (1|ID)
imp <- panImpute(studentratings, formula = fml, n.burn = 1000, n.iter = 100, m = 5)

# print summary
summary(imp)
}
\keyword{methods}
