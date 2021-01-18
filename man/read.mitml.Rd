\name{read.mitml}
\alias{read.mitml}
\title{Read \code{mitml} objects from file}
\description{
This function loads \code{mitml} class objects from R binary formats (similar to \code{?load}), usually produced by \code{write.mitml}.
}

\usage{

read.mitml(filename)

}

\arguments{

  \item{filename}{Name of the file to read, to be specified with file extension (e.g., \code{.R}, \code{.Rdata}).}

}

\value{

Returns the saved \code{mitml} class object.

}

\author{Simon Grund}
\seealso{\code{\link{panImpute}}, \code{\link{jomoImpute}}, \code{\link{write.mitml}}}
\examples{
data(studentratings)

fml <- ReadDis + SES ~ ReadAchiev + (1|ID)
imp <- panImpute(studentratings, formula = fml, n.burn = 1000, n.iter = 100, m = 5)

# write 'mitml' object
write.mitml(imp, filename = "imputation.R")

# read previously saved 'mitml' object
previous.imp <- read.mitml("imputation.R")

class(previous.imp)
previous.imp
}
