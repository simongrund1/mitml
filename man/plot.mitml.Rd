\name{plot.mitml}
\alias{plot.mitml}
\title{Print diagnostic plots}
\description{
Generates diagnostic plots for assessing the convergence and autocorrelation behavior of \code{pan}'s and \code{jomo}'s MCMC algorithms.
}

\usage{

\method{plot}{mitml}(x, print = c("beta", "beta2", "psi", "sigma"), pos = NULL, group = "all",
  trace = c("imputation", "burnin", "all"), thin = 1, smooth = 3, n.Rhat = 3,
  export = c("none", "png", "pdf"), dev.args = list(), ...)

}

\arguments{

  \item{x}{An object of class \code{mitml} as produced by \code{panImpute} or \code{jomoImpute}.}
  \item{print}{A character vector containing one or several of \code{"beta"}, \code{"beta2"}, \code{"psi"} or \code{"sigma"} denoting which parameters should be plotted. Default is to plot all parameters.}
  \item{pos}{Either \code{NULL} or an integer vector denoting a specific entry in \code{"beta"}, \code{"beta2"}, \code{"psi"} or \code{"sigma"}. Default is \code{NULL}, which plots all entries.}
  \item{group}{Either \code{"all"} or an integer denoting for which group the plots should be generated. Used only when groupwise imputation was used. Default is \code{"all"}.}
  \item{trace}{One of \code{"imputation"}, \code{"burnin"} or \code{"all"} denoting which part of the chain should be used for the trace plot. Default is \code{"imputation"}, which plots only the iterations after burn-in.}
  \item{thin}{An integer denoting the thinning factor that is applied before plotting. Default is \code{1}, which corresponds to no thinning.}
  \item{smooth}{A numeric value denoting the smoothing factor for the trend line in trace plots. Higher values correspond to less smoothing. Default is \code{3}. If set to \code{0} or \code{NULL}, the trend line is suppressed.}
  \item{n.Rhat}{An integer denoting the number of segments of each chain used for calculating the potential scale reduction factor. Default is \code{3}.}
  \item{export}{(optional) A character string specifying if plots should be exported to a file. If \code{"png"} or \code{"pdf"}, then plots are printed into a folder named "mitmlPlots" in the current directory using either the \code{png} or the \code{pdf} device. Default is \code{"none"}, which does not export files.}
  \item{dev.args}{(optional) A named list containing additional arguments that are passed to the graphics device.}
  \item{\dots}{Parameters passed to the plotting functions.}

}

\details{

The \code{plot} method generates a series of plots for the parameters of the imputation model which can be used for diagnostic purposes.
In addition, a short summary of the parameter chain is displayed.

Setting \code{print} to \code{"beta"}, \code{"beta2"}, \code{"psi"} and \code{"sigma"} will plot the fixed effects, the variances and covariances of the random effects, and the variances and covariances of the residuals, respectively.
In this context, \code{"beta2"} refers to the fixed effects for target variables at level 2 and is only used when a two-part model was specified in (\code{\link{jomoImpute}}).
Each plotting window contains a trace plot (upper left), an autocorrelation plot (lower left), a kernel density approximation of the posterior distribution (upper right), and a posterior summary (lower right).
The summary includes the following quantities:
\describe{
  \item{\code{EAP}:}{Expected value a posteriori (i.e., the mean of the parameter chain)}
  \item{\code{MAP}:}{Mode a posteriori (i.e., the mode of the parameter chain)}
  \item{\code{SD}:}{Standard deviation of the parameter chain}
  \item{\code{2.5\%}:}{The 2.5\% quantile of parameter values}
  \item{\code{97.5\%}:}{The 97.5\% quantile of parameter values}
  \item{\code{Rhat}:}{Estimated potential scale reduction factor (\eqn{\hat{R}})}
  \item{\code{ACF-k}:}{Smoothed autocorrelation at lag \eqn{k}, where \eqn{k} is the number of iterations between imputations (see \code{\link{summary.mitml}})}
}
The \code{trace} and \code{smooth} arguments can be used to influence how the trace plot is drawn and what part of the chain should be used for it.
The \code{thin} argument can be used for thinning the chain before plotting, in which case the number of data points is reduced in the trace plot, and the autocorrelation is calculated up to lag \eqn{k/thin} (see above).
The \code{n.Rhat} argument controls the number of segments that are used for calculating the potential scale reduction factor (\eqn{\hat{R}}) in each plot (see \code{summary.mitml}).
Further aguments to the graphics device are supplied using the \code{dev.args} argument.

The \code{plot} function computes and displays diagnostic information primarily for the imputation phase (i.e., for iterations after burn-in).
This is the default in the \code{plot} function and the recommended method for most users.
If \code{trace = "all"}, the full chain is displayed with emphasis on the imputation phase, and the posterior summary is calculated based on only the iterations after burn-in (as recommended).
If \code{trace = "burnin"}, the posterior summary and the trace plots are calculated based on only the burn-on interations, which is generally not sufficient to establish convergence and should be used with caution.

}
\note{

The plots are presented on-screen one at a time.
To proceed with the next plot, the user may left-click in the plotting window or press the "enter" key while in the R console, depending on the operating system.
No plots are displayed when exporting to file.

}

\value{

None (invisible \code{NULL}).

}

\author{Simon Grund}
\seealso{\code{\link{panImpute}}, \code{\link{jomoImpute}}, \code{\link{summary.mitml}}}
\examples{
\dontrun{
data(studentratings)

# * Example 1: simple imputation

fml <- ReadDis + SES ~ ReadAchiev + (1|ID)
imp <- panImpute(studentratings, formula = fml, n.burn = 1000, n.iter = 100, m = 5)

# plot all parameters (default)
plot(imp)

# plot fixed effects only
plot(imp, print = "beta")

# export plots to file (using pdf device)
plot(imp, export = "pdf", dev.args = list(width = 9, height = 4, pointsize = 12))

# * Example 2: groupwise imputation

fml <- ReadDis + SES ~ ReadAchiev + (1|ID)
imp <- panImpute(studentratings, formula = fml, group = FedState, n.burn = 1000,
       n.iter = 100, m = 5)

# plot fixed effects for all groups (default for 'group')
plot(imp, print = "beta", group = "all")

# plot fixed effects for first group only
plot(imp, print = "beta", group = 1)
}
}
\keyword{methods}
