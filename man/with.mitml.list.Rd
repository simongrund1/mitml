\name{with.mitml.list}
\alias{with.mitml.list}
\alias{within.mitml.list}
\title{Evaluate an expression in a list of imputed data sets}
\description{
The functions \code{with} and \code{within} evaluate R expressions in a list of multiply imputed data sets.
}

\usage{

\method{with}{mitml.list}(data, expr, ...)
\method{within}{mitml.list}(data, expr, ignore=NULL, ...)

}

\arguments{

  \item{data}{A list of imputed data sets with class \code{mitml.list} as produced by \code{mitmlComplete}.}
  \item{expr}{An R expression to be evaluated for each data set.}
  \item{ignore}{A character vector denoting objects not to be saved.}
  \item{\dots}{Not being used.}

}

\details{

The two functions are defined as \code{with} and \code{within} methods for objects of class \code{mitml.list}. Both \code{with} and \code{within} evaluate an R expression in each of the imputed data sets.
However, the two functions return different values:
\code{with} returns the evaluated expression, whereas \code{within} returns the resuling data sets.
The \code{ignore} argument may be used to declare objects that are not to be saved within \code{within}.

}

\value{

% \itemize{
%   \item{\code{with}}{: Returns the evaluated expression from each data set as a list (class \code{mitml.result}. This is useful for fitting statistical models to multiply imputed data. The list of fitted models can be analyzed using \code{testEstimates}, \code{testModels}, \code{testConstraints}, or \code{anova}.}
%   \item{\code{within}}{: Evaluates the R expression for each data set and returns the altered data sets as a list \code{mitml.list}. This is useful for manipulating the data prior to data analysis (e.g., centering, calculating cluster means, etc.).}
% }

\code{with}: Returns the evaluated expression from each data set as a list (class \code{mitml.result}). This is useful for fitting statistical models to multiply imputed data. The list of fitted models can be analyzed using \code{\link{testEstimates}}, \code{\link{testModels}}, \code{\link{testConstraints}}, or \code{\link{anova}}.

\code{within}: Evaluates the R expression for each data set and returns the altered data sets as a list \code{mitml.list}. This is useful for manipulating the data prior to data analysis (e.g., centering, calculating cluster means, etc.).

}

\author{Simon Grund}
\seealso{\code{\link{mitmlComplete}}, \code{\link{anova.mitml.result}}, \code{\link{testEstimates}}, \code{\link{testModels}}, \code{\link{testConstraints}}}
\examples{
data(studentratings)

fml <- ReadDis + SES ~ ReadAchiev + (1|ID)
imp <- panImpute(studentratings, formula=fml, n.burn=1000, n.iter=100, m=5)

implist <- mitmlComplete(imp, print=1:5)

# * Example 1: data manipulation

# calculate and save cluster means
new1.implist <- within(implist, Means.ReadAchiev <- clusterMeans(ReadAchiev, ID))

# center variables, calculate interaction terms, ignore byproducts
new2.implist <- within(implist,{
  M.SES <- mean(SES)
  M.CognAbility <- mean(CognAbility)
  C.SES <- SES - M.SES
  C.CognAbility <- CognAbility - M.CognAbility
  SES.CognAbility <- C.SES * C.CognAbility
}, ignore=c("M.SES", "M.CognAbility"))

# * Example 2: fitting statistical models

# fit regression model
fit.lm <- with(implist, lm(ReadAchiev ~ ReadDis))

# fit multilevel model using lme4
require(lme4)
fit.lmer <- with(implist, lmer(ReadAchiev ~ ReadDis + (1|ID)))

# * Example 3: manual extraction of variance estimates
require(lme4)
fit.lmer <- with(implist, lmer(SES ~ (1|ID)))

# extract level-1 and level-2 variances
var.l1 <- sapply(fit.lmer, function(z) attr(VarCorr(z),"sc")^2)
var.l2 <- sapply(fit.lmer, function(z) VarCorr(z)$ID[1,1])

# calculate final estimate of the intraclass correlation
ICC <- mean( var.l2 / (var.l2+var.l1) )
}
\keyword{methods}
