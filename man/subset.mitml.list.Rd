\name{subset.mitml.list}
\alias{subset.mitml.list}
\title{Subset a list of imputed data sets}
\description{
The functions can be used for creating subsets for a list of multiply imputed data sets.
}

\usage{

\method{subset}{mitml.list}(x, subset, select, ...)

}

\arguments{

  \item{x}{A list of imputed data sets with class \code{mitml.list} as produced by \code{mitmlComplete} (or similar).}
  \item{subset}{An R expression by which to subset each data set.}
  \item{select}{An R expression by which to select columns.}
  \item{\dots}{Not being used.}

}

\details{

This function can be used to create subsets and select variables for a list of multiply imputed data sets according to the R expressions given in the \code{subset} and \code{select} arguments.
The function is similar to and adapted from the \code{subset} function for regular data sets.
Note that subsetting is performed individually for each data set.
Thus, the cases included may differ across data sets if the variables used for subsetting contain different values.

}

\value{

A list of imputed data sets with an additional class attribute \code{mitml.list}.

}

\author{Simon Grund}
\examples{
data(studentratings)

fml <- ReadDis + SES ~ ReadAchiev + (1|ID)
imp <- panImpute(studentratings, formula=fml, n.burn=1000, n.iter=100, m=5)

implist <- mitmlComplete(imp,"all")

# * Example 1: subset by SES, select variables by name
subset(implist, SES < 25, select = c(ID, FedState, Sex, SES, ReadAchiev, ReadDis))

# * Example 2: subset by FedState, select variables by column number
subset(implist, FedState == "SH", select = -c(6:7,9:10))

\dontrun{
# * Example 3: subset by ID and Sex
subset(implist, ID %in% 1001:1005 & Sex == "Boy")

# * Example 4: select variables by name range
subset(implist, select = ID:Sex)
}
}
\keyword{methods}
