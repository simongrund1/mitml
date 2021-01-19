\name{mitmlComplete}
\alias{mitmlComplete}
\title{Extract imputed data sets}
\description{
This function extracts imputed data sets from \code{mitml} class objects as produced by \code{panImpute} and \code{jomoImpute}.
}

\usage{

mitmlComplete(x, print = "all", force.list = FALSE)

}

\arguments{

  \item{x}{An object of class \code{mitml} as produced by \code{panImpute} or \code{jomoImpute}.}
  \item{print}{Either an integer vector, \code{"list"}, or \code{"all"} denoting which data sets to extract. If set to \code{"list"} or \code{"all"}, all imputed data sets will be returned as a list. Negative values and zero return the original (incomplete) data set. Default is \code{"all"}.}
  \item{force.list}{(optional) Logical flag indicating if single data sets should be enclosed in a list. Default is \code{FALSE}.}

}

\value{

Usually a list of imputed data with class \code{mitml.list}
If only one data set is extracted: a data frame unless \code{force.list = TRUE}.

}

\author{Simon Grund}
\seealso{\code{\link{panImpute}}, \code{\link{jomoImpute}}}

\examples{
data(studentratings)

fml <- ReadDis + SES ~ ReadAchiev + (1|ID)
imp <- panImpute(studentratings, formula = fml, n.burn = 1000, n.iter = 100, m = 5)

# extract original (incomplete) data set
mitmlComplete(imp, print = 0)

# extract first imputed data set (returned as mitml.list)
mitmlComplete(imp, print = 1, force.list = TRUE)

# extract all imputed data sets at once
implist <- mitmlComplete(imp, print = "all")

\dontrun{
# ... alternatives with same results
implist <- mitmlComplete(imp, print = 1:5)
implist <- mitmlComplete(imp, print = "list")
}
}
