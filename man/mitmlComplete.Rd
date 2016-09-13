\name{mitmlComplete}
\alias{mitmlComplete}
\title{Extract imputed data sets}
\description{This function extracts imputed data sets from \code{mitml} class objects as produced by \code{panImpute} and \code{jomoImpute}.}
\usage{

mitmlComplete(x, print="all", force.list=FALSE)

}
\arguments{

  \item{x}{An object of class \code{mitml} as produced by \code{panImpute} and \code{jomoImpute}.}
  \item{print}{Either an integer vector, \code{"list"}, or \code{"all"} denoting which data sets to extract. If set to \code{"list"} or \code{"all"}, then all imputed data sets will be returned as a list. Negative values and zero will return the original (incomplete) data set. Default is to \code{"all"}.}
  \item{force.list}{(optional) Logical flag indicating if single data sets should be enclosed in a list. Default is to \code{FALSE}.}

}
\value{

Single data sets are returned as a data frame unless \code{force.list=TRUE}.
If several data sets are extracted, the result is always a list of data sets with an additional class attribute \code{mitml.list}.

}
\author{Simon Grund}

\seealso{\code{\link{panImpute}}, \code{\link{jomoImpute}}}
\examples{
data(studentratings)

fml <- ReadDis + SES ~ ReadAchiev + (1|ID)
imp <- panImpute(studentratings, formula=fml, n.burn=1000, n.iter=100, m=5)

# extract original (incomplete) data set
mitmlComplete(imp, print=0)

# extract first imputed data set (returned as mitml.list)
mitmlComplete(imp, print=1, force.list=TRUE)

# extract all imputed data sets at once
implist <- mitmlComplete(imp, print="all")

\dontrun{
# ... alternatives with same results
implist <- mitmlComplete(imp, print=1:5)
implist <- mitmlComplete(imp, print="list")
}
}
