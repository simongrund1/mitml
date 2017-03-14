\name{amelia2mitml.list}
\alias{amelia2mitml.list}
\title{Convert objects of class \code{amelia} to \code{mitml.list}}
\description{This function converts a \code{amelia} class object (as produced by the \code{Amelia} package) to \code{mitml.list}. The resulting object may be used in further analyses.}
\usage{

amelia2mitml.list(x)

}
\arguments{

  \item{x}{An object of class \code{amelia} as produced by \code{amelia} (see the \code{Amelia} package).}

}
\value{

A list of imputed data sets with an additional class attribute \code{mitml.list}.

}
\author{Simon Grund}
\seealso{\code{\link{mitmlComplete}}}

\examples{
data(studentratings)

require(Amelia)
imp <- amelia(x=studentratings[,c("ID","MathAchiev","ReadAchiev")], cs="ID")

implist <- amelia2mitml.list(imp)
}

