\name{write.mitmlMplus}
\alias{write.mitmlMplus}
\title{Write \code{mitml} objects to Mplus format}
\description{
Saves objects of class \code{mitml} as a series of text files which can be processed by the statistical software M\emph{plus} (Muthen & Muthen, 2012).
}

\usage{

write.mitmlMplus(x, filename, suffix = "list", sep = "\t", dec = ".",
  na.value = -999)

}

\arguments{

  \item{x}{An object of class \code{mitml} or \code{mitml.list} as produced by \code{panImpute}, \code{jomoImpute}, \code{mitmlComplete}, or similar).}
  \item{filename}{File base name for the text files containing the imputed data sets, specified without file extension.}
  \item{suffix}{File name suffix for the index file.}
  \item{sep}{The field separator.}
  \item{dec}{The decimal separator.}
  \item{na.value}{A numeric value coding the missing data in the resulting data files.}

}

\details{

The M\emph{plus} format for multiply imputed data sets comprises a set of text files, each containing one imputed data set, and an index file containing the names of all data files.
During export, factors and character variables are converted to numeric.
To make this more transparent, \code{write.mitmlMplus} produces a log file which contains information about the data set and the factors that have been converted.

In addition, a basic M\emph{plus} input file is generated that can be used for setting up subsequent analysis models.

}

\value{

None (invisible \code{NULL}).

}

\references{
Muthen, L. K., & Muthen, B. O. (2012). \emph{Mplus User's Guide. Seventh Edition.} Los Angeles, CA: Muthen & Muthen.
}

\author{Simon Grund}
\seealso{\code{\link{panImpute}}, \code{\link{jomoImpute}}, \code{\link{mitmlComplete}}}
\examples{
\dontrun{
data(studentratings)

fml <- ReadDis + SES ~ ReadAchiev + (1|ID)
imp <- panImpute(studentratings, formula = fml, n.burn = 1000, n.iter = 100, m = 5)

# write imputation files, index file, and log file
write.mitmlMplus(imp, filename = "imputation", suffix = "list", na.value = -999)
}
}
