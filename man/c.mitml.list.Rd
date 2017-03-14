\name{c.mitml.list}
\alias{c.mitml.list}
\alias{rbind.mitml.list}
\alias{cbind.mitml.list}
\title{Concatenate lists of imputed data sets}
\description{
These functions allow concatenating lists of imputed data sets by data set, row, or column.
}

\usage{

\method{c}{mitml.list}(...)
\method{rbind}{mitml.list}(...)
\method{cbind}{mitml.list}(...)

}

\arguments{

  \item{\dots}{One or several lists of imputed data sets with class \code{mitml.list} as produced by \code{mitmlComplete} (or similar).}

}

\details{

These function allow concatenating multiple lists of imputed data sets.
The function \code{c} concatenates by data set (i.e., by appending additional data sets to the list), \code{rbind} concatenates by row (i.e., appending additional rows to each data set), and \code{cbind} concatenates by column (i.e., by appending additional columns to each data set).

These functions are intended for experienced users and should be used with caution.
Appending rows or columns from multiple imputation procedures is usually unsafe unless in special applications (see Examples).
}

\value{

A list of imputed data sets with an additional class attribute \code{mitml.list}.

}

\author{Simon Grund}
\examples{
# Example 1: manual imputation by grouping variable

data(studentratings)
fml <- ReadDis + SES ~ ReadAchiev + (1|ID)

imp <- panImpute(subset(studentratings, FedState=="SH"), formula=fml, n.burn=1000, n.iter=100, m=5)
implist <- mitmlComplete(imp, print="all")

imp2 <- panImpute(subset(studentratings, FedState=="B"), formula=fml, n.burn=1000, n.iter=100, m=5)
implist2 <- mitmlComplete(imp2, print="all")

rbind(implist, implist2)

# Example 2: predicted values from linear model

pred <- with(implist, predict(lm(ReadDis~ReadAchiev)))
cbind(implist, pred.ReadDis=pred)
}
\keyword{methods}
