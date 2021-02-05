\name{with.mitml.list}
\alias{with.mitml.list}
\alias{within.mitml.list}
\title{Evaluate an expression in a list of imputed data sets}
\description{
The \code{with} and \code{within} methods evaluate R expressions in a list of multiply imputed data sets.
}

\usage{

\method{with}{mitml.list}(data, expr, include.data = FALSE, ...)
\method{within}{mitml.list}(data, expr, ignore = NULL, ...)

}

\arguments{

  \item{data}{A list of imputed data sets with class \code{mitml.list} as produced by \code{mitmlComplete} or \code{as.mitml.list}.}
  \item{expr}{An R expression to be evaluated for each data set.}
  \item{include.data}{Either a logical flag or a character string denoting how the data are included when \code{expr} is evaluated (see 'Details'). If \code{FALSE}, an environment is created from \code{data}, and \code{expr} is evaluated therein. If \code{TRUE}, a call is constructed from \code{expr} and evaluated with the imputed data in the \code{"data"} slot. If character, a call is constructed from \code{expr} and evaluated with the imputed data in the slot named by \code{include.data}. Default is \code{FALSE}.}
  \item{ignore}{A character vector naming objects that are created but should not be saved (see 'Details').}
  \item{\dots}{Not used.}

}

\details{

The two functions provide \code{with} and \code{within} methods for objects of class \code{mitml.list}.
They evaluate an R expression repeatedly for each of the imputed data sets but return different values: \code{with} returns the result of the evaluated expression; \code{within} returns the resulting data sets.

The \code{within} function is useful for transforming and computing variables in the imputed data (e.g., centering, calculating cluster means, etc.).
The \code{with} function is useful, for example, for fitting statistical models.
The list of fitted models can be analyzed using \code{\link{testEstimates}}, \code{\link{testModels}}, \code{\link{testConstraints}}, or \code{\link[=anova.mitml.result]{anova}}.

The \code{include.data} argument can be used to include the imputed data sets in the call to fit statistical models (\code{expr}) using \code{with}.
This is useful for fitting models that require that the fitting function be called with a proper \code{data} argument (e.g., \code{lavaan} or \code{nlme}; see 'Examples'). 
Setting \code{include.data = TRUE} will fit the model with the imputed data sets used as the \code{data} argument.
Setting \code{include.data = "df"} (or similar) will fit the model with the imputed data sets as the \code{df} argument (useful if the function refers to the data by a nonstandard name, such as \code{"df"}).

The \code{ignore} argument can be used to declare objects that are not to be saved in the data sets created by \code{within}.

}

\value{

\code{with}: A list of class \code{mitml.results} containing the evaluated expression for each data set. 

\code{within}: A list of class \code{mitml.list} containing the imputed data modified by the evaluated expression.

}

\author{Simon Grund}
\seealso{\code{\link{mitmlComplete}}, \code{\link{anova.mitml.result}}, \code{\link{testEstimates}}, \code{\link{testModels}}, \code{\link{testConstraints}}}
\examples{
data(studentratings)

fml <- ReadDis + SES ~ ReadAchiev + (1|ID)
imp <- panImpute(studentratings, formula = fml, n.burn = 1000, n.iter = 100, m = 5)

implist <- mitmlComplete(imp)

# * Example 1: data transformation

# calculate and save cluster means
new1.implist <- within(implist, Means.ReadAchiev <- clusterMeans(ReadAchiev, ID))

# center variables, calculate interaction terms, ignore byproducts
new2.implist <- within(implist, {
  M.SES <- mean(SES)
  M.CognAbility <- mean(CognAbility)
  C.SES <- SES - M.SES
  C.CognAbility <- CognAbility - M.CognAbility
  SES.CognAbility <- C.SES * C.CognAbility
}, ignore = c("M.SES", "M.CognAbility"))

# * Example 2: fitting statistical models

# fit regression model
fit.lm <- with(implist, lm(ReadAchiev ~ ReadDis))

# fit multilevel model with lme4
require(lme4)
fit.lmer <- with(implist, lmer(ReadAchiev ~ ReadDis + (1|ID)))

\dontrun{
# fit structural equation model with lavaan (with include.data = TRUE)
require(lavaan)
mod <- "ReadAchiev ~ ReadDis"
fit.sem <- with(implist,
                sem(model = mod, cluster = "ID", estimator = "MLR"),
                include.data = TRUE)
}
}
\keyword{methods}
