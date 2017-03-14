\name{long2mitml.list}
\alias{long2mitml.list}
\alias{jomo2mitml.list}
\title{Convert imputations from long format to \code{mitml.list}}
\description{These functions convert data sets containing multiple imputations in long format to objects of class \code{mitml.list}. The resulting object can be used in further analyses.}
\usage{

long2mitml.list(x, split, exclude=NULL)

jomo2mitml.list(x)

}
\arguments{

  \item{x}{A data frame in long format containing multiple imputations (see details).}
  \item{split}{A character string denoting the column in \code{x} that identifies different imputations (see details).}
  \item{exclude}{A vector denoting values of \code{split} which should be excluded from the list.}

}
\details{

The function \code{long2mitml.list} is intedended for converting data frames from the long format to \code{mitml.list} (i.e., a list of imputed data sets).
In this format, imputations are enclosed in a single data frame, and a \code{split} variable is used to identify different imputations.

Similarly, \code{jomo2mitml.list} is a special case of \code{long2mitml.list} intended for converting imputations that have been generated with \code{jomo} directly.
}
\value{

A list of imputed data sets with an additional class attribute \code{mitml.list}.

}
\author{Simon Grund}
\seealso{\code{\link{mitmlComplete}}}

\examples{
data(studentratings)
require(jomo)

# impute data using jomo (native functions)
clus <- studentratings[,"ID"]
Y <- studentratings[,c("ReadAchiev","ReadDis")]
imp <- jomo(Y=Y, clus=clus, nburn=1000, nbetween=100, nimp=5)

# split imputations
impList <- long2mitml.list(imp, split="Imputation", exclude=0)
impList <- jomo2mitml.list(imp)
}
