\name{mitml.list2mids}
\alias{mitml.list2mids}
\title{Convert objects of class \code{mitml.list} to \code{mids}}
\description{
This function converts a \code{mitml.list} class object to \code{mids} (as used in the \code{mice} package).
}

\usage{
mitml.list2mids(x, data, fill = FALSE, where = NULL)
}

\arguments{


  \item{x}{A list of imputed data sets with class \code{mitml.list} (as produced by \code{\link{mitmlComplete}}, \code{\link{mids2mitml.list}}, or similar).}
  \item{data}{A data frame containing to original (incomplete) data (see 'Details').}
  \item{fill}{A logical flag indicating whether variables in the imputed data that are not in the original data should be added and filled with \code{NA} (default is \code{FALSE}).}
  \item{where}{(optional) A data frame or matrix of logicals indicating the location of missing values (see 'Details').}

}

\details{

This function converts objects of class \code{mitml.list} into \code{mids} objects (as used in the \code{mice} package).
The conversion requires a list of imputed data sets and the original (incomplete) data set.
If the imputed data sets have been appended with new variables (e.g., by \code{\link{within.mitml.list}}), the new variables can be added to the original data set by setting \code{fill = TRUE}.

This function is essentially a wrapper around \code{\link[mice:as.mids]{as.mids}} that sets the case and imputation identifiers automatically and and passes the \code{where} argument as is (see also the documentation of \code{\link[mice:as.mids]{as.mids}}).

}

\value{

An object of class \code{mids}.

}

\author{Simon Grund}
\seealso{\code{\link{mitmlComplete}}, \code{\link{mids2mitml.list}}, \code{\link{within.mitml.list}}}

\examples{
data(studentratings)

fml <- ReadDis + SES ~ ReadAchiev + (1|ID)
imp <- panImpute(studentratings, formula = fml, n.burn = 1000, n.iter = 100, m = 5)

implist <- mitmlComplete(imp)

# * Example 1: simple conversion

# convert to 'mids'
impmids <- mitml.list2mids(implist, data = studentratings)

# * Example 2: conversion with additional variables

# compute new variables
implist <- within(implist, {
  M.ReadAchiev <- clusterMeans(ReadAchiev, ID)
  C.ReadAchiev <- ReadAchiev - M.ReadAchiev
})

# convert to 'mids'
impmids <- mitml.list2mids(implist, data = studentratings, fill = TRUE)
}
