\name{write.mitmlSAV}
\alias{write.mitmlSAV}
\title{Write \code{mitml} objects to native SPSS format}
\description{
Saves objects of class \code{mitml} in the \code{.sav} format used by the statistical software SPSS (IBM Corp., 2013).
The function serves as a front-end for \code{write_sav} from the \code{haven} package.
}

\usage{

write.mitmlSAV(x, filename)

}

\arguments{

  \item{x}{An object of class \code{mitml} or \code{mitml.list} as produced by \code{panImpute}, \code{jomoImpute}, \code{mitmlComplete}, or similar).}
  \item{filename}{Name of the destination file, to be specified with or without file extension. The file extension (\code{.sav}) is appended if necessary.}

}

\details{

This function exports multiply imputed data sets to a single \code{.sav} file, in which an \code{Imputation_} variable separates the original data and the various imputed data sets.
Thus, \code{write.mitmlSAV} exports directly to the native SPSS format.

Alternatively, \code{\link{write.mitmlSPSS}} may be used for creating separate text and SPSS syntax files; an option that offers more control over the data format.

}

\value{

None (invisible \code{NULL}).

}

\references{
IBM Corp. (2013). \emph{IBM SPSS Statistics for Windows, Version 22.0}. Armonk, NY: IBM Corp
}

\author{Simon Grund}
\seealso{\code{\link{panImpute}}, \code{\link{jomoImpute}}, \code{\link{mitmlComplete}}, \code{\link{write.mitmlSPSS}}}
\examples{
data(studentratings)

fml <- ReadDis + SES ~ ReadAchiev + (1|ID)
imp <- panImpute(studentratings, formula = fml, n.burn = 1000, n.iter = 100, m = 5)

# write data file and SPSS syntax
write.mitmlSAV(imp, filename = "imputation")
}
