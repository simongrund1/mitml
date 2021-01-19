\name{sort.mitml.list}
\alias{sort.mitml.list}
\title{Sort a list of imputed data sets}
\description{
Sorts a list of multiply imputed data sets according to an R expression.
}

\usage{

\method{sort}{mitml.list}(x, decreasing = FALSE, by, ...)

}

\arguments{

  \item{x}{A list of imputed data sets with class \code{mitml.list} as produced by \code{mitmlComplete} (or similar).}
  \item{decreasing}{Logical flag indicating if data sets should be sorted in decreasing (i.e., reversed) order. Default is `\code{FALSE}`.}
  \item{by}{An R expression or a list of multiple expressions by which to sort the imputed data sets (see 'Examples').}
  \item{\dots}{Further arguments to `\code{order}' (see 'Details').}

}

\details{

This function sorts a list of imputed data sets according to the R expression given in the \code{by} argument.
The function is based on \code{order} and works in a similar manner.
Note that sorting is performed individually for each data set.
For this reason, the order of cases may differ across data sets if the variables used for sorting contain different values.

}

\value{

A list of imputed data sets with class \code{mitml.list}.

}

\author{Simon Grund}
\examples{
data(studentratings)

fml <- ReadDis + SES ~ ReadAchiev + (1|ID)
imp <- panImpute(studentratings, formula = fml, n.burn = 1000, n.iter = 100, m = 5)

implist <- mitmlComplete(imp)

# * Example 1: sort by ID
sort(implist, by = ID)

# * Example 2: sort by combination of variables
sort(implist, by = list(FedState, ID, -SES))
}
\keyword{methods}
