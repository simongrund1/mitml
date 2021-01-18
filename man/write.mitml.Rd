\name{write.mitml}
\alias{write.mitml}
\title{Write \code{mitml} objects to file}
\description{
This function saves objects of class \code{mitml} in R binary formats (similar to \code{?save}).
}

\usage{

write.mitml(x, filename, drop = FALSE)

}

\arguments{

  \item{x}{An object of class \code{mitml} as produced by \code{panImpute} and \code{jomoImpute}.}
  \item{filename}{Name of the destination file, to be specified with file extension (e.g., \code{.R}, \code{.Rdata}).}
  \item{drop}{Logical flag indicating if the parameters of the imputation model should be dropped in favor for lower file size. Default is to \code{FALSE}.}

}

\value{

None (invisible \code{NULL}).

}

\author{Simon Grund}
\seealso{\code{\link{panImpute}}, \code{\link{jomoImpute}}, \code{\link{read.mitml}}}
\examples{
data(studentratings)

fml <- ReadDis + SES ~ ReadAchiev + (1|ID)
imp <- panImpute(studentratings, formula = fml, n.burn = 1000, n.iter = 100, m = 5)

# write full 'mitml' object (default)
write.mitml(imp, filename = "imputation.R")

# drop parameters of the imputation model
write.mitml(imp, filename = "imputation.R", drop = TRUE)
}
