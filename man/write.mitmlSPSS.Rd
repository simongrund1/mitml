\name{write.mitmlSPSS}
\alias{write.mitmlSPSS}
\title{Write \code{mitml} objects to SPSS compatible format}
\description{
Saves objects of class \code{mitml} as a text and a syntax file which can be processed by the statistical software SPSS (IBM Corp., 2013).
}

\usage{

write.mitmlSPSS(x, filename, sep = "\t", dec = ".", na.value = -999, syntax = TRUE,
  locale = NULL)

}

\arguments{

  \item{x}{An object of class \code{mitml} or \code{mitml.list} as produced by \code{panImpute}, \code{jomoImpute}, \code{mitmlComplete}, or similar).}
  \item{filename}{File base name of the data and syntax files, specified without file extension.}
  \item{sep}{The field separator.}
  \item{dec}{The decimal separator.}
  \item{na.value}{A numeric value coding the missing data in the resulting data file.}
  \item{syntax}{A logical flag indicating if an SPSS syntax file should be generated. This file contains instructions for SPSS for reading in the data file. Default is \code{TRUE}.}
  \item{locale}{(optional) A character string specifying the localization to be used in SPSS (e.g., \code{"en_US"}, \code{"de_DE"}; see 'Details').}

}

\details{

In SPSS, multiply imputed data are contained in a single file, in which an \code{Imputation_} variable separates the original data and the various imputed data sets.
During export, factors are converted to numeric, whereas character variables are left ``as is''.

By default, \code{write.mitmlSPSS} generates a raw text file containing the data, along with a syntax file containing instructions for SPSS.
This syntax file mimics SPSS's functionality to read text files with sensible settings.
In order to read in the data, the syntax file must be opened and executed using SPSS, or open using the GUI.
Manual changes to the syntax file can be required, for example, if the file path of the data file is not correctly represented in the syntax.
The \code{locale} argument can be used to ensure that SPSS reads the data in the correct locale.

Alternatively, \code{\link{write.mitmlSAV}} may be used for exporting directly to the SPSS native \code{.sav} format.
}

\value{

None (invisible \code{NULL}).

}

\references{
IBM Corp. \emph{IBM SPSS Statistics for Windows}. Armonk, NY: IBM Corp
}

\author{Simon Grund}
\seealso{\code{\link{panImpute}}, \code{\link{jomoImpute}}, \code{\link{mitmlComplete}}, \code{\link{write.mitmlSAV}}}
\examples{
data(studentratings)

fml <- ReadDis + SES ~ ReadAchiev + (1|ID)
imp <- panImpute(studentratings, formula = fml, n.burn = 1000, n.iter = 100, m = 5)

# write data file and SPSS syntax
write.mitmlSPSS(imp, filename = "imputation", sep = "\t", dec = ".",
                na.value = -999, locale = "en_US")
}
