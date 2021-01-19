\name{mids2mitml.list}
\alias{mids2mitml.list}
\title{Convert objects of class \code{mids} to \code{mitml.list}}
\description{
This function converts a \code{mids} class object (as produced by the \code{mice} package) to \code{mitml.list}. The resulting object may be used in further analyses.
}

\usage{

mids2mitml.list(x)

}

\arguments{

  \item{x}{An object of class \code{mids} as produced by \code{mice} (see the \code{mice} package).}

}

\value{

A list of imputed data sets with class \code{mitml.list}.

}

\author{Simon Grund}
\seealso{\code{\link{mitmlComplete}}}

\examples{
data(studentratings)

# imputation using mice
require(mice)
imp <- mice(studentratings)

implist <- mids2mitml.list(imp)
}
