\name{long2mitml.list}
\alias{long2mitml.list}
\alias{jomo2mitml.list}
\title{Convert imputations from long format to \code{mitml.list}}
\description{
These functions convert data sets containing multiple imputations in long format to objects of class \code{mitml.list}. The resulting object can be used in further analyses.
}

\usage{

long2mitml.list(x, split, exclude = NULL)

jomo2mitml.list(x)

}

\arguments{

  \item{x}{A data frame in long format containing multiple imputations (see 'Details').}
  \item{split}{A character string denoting the column in \code{x} that identifies different imputations (see 'Details').}
  \item{exclude}{A vector denoting the values of \code{split} that should be excluded.}

}

\details{

The function \code{long2mitml.list} converts data frames from the long format to \code{mitml.list} (i.e., a list of imputed data sets).
In long format, all imputations are contained in a single data frame, where different imputations are denoted by \code{split}.
This function splits the data frame into a list of imputed data sets according to \code{split}, excluding the values specified by \code{exclude} (see the 'Examples').

The \code{jomo2mitml.list} function is a special case of \code{long2mitml.list} which converts imputations that have been generated with \code{jomo} (see the \code{jomo} package)).
}

\value{

A list of imputed data sets with class \code{mitml.list}.

}

\author{Simon Grund}
\seealso{\code{\link{mitmlComplete}}}

\examples{
data(studentratings)
require(jomo)

# impute data using jomo (native functions)
clus <- studentratings[, "ID"]
Y <- studentratings[, c("ReadAchiev", "ReadDis")]
imp <- jomo(Y = Y, clus = clus, nburn = 1000, nbetween = 100, nimp = 5)

# split imputations
impList <- long2mitml.list(imp, split = "Imputation", exclude = 0)
impList <- jomo2mitml.list(imp)
}
