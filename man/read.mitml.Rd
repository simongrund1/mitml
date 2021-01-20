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

An object of class \code{mitml}.

}

\author{Simon Grund}
\seealso{\code{\link{panImpute}}, \code{\link{jomoImpute}}, \code{\link{write.mitml}}}
\examples{
\dontrun{
data(studentratings)

fml <- ReadDis + SES ~ ReadAchiev + (1|ID)
imp <- panImpute(studentratings, formula = fml, n.burn = 1000, n.iter = 100, m = 5)

# write 'mitml' object
write.mitml(imp, filename = "imputation.R")

# read previously saved 'mitml' object
old.imp <- read.mitml("imputation.R")

class(old.imp)
old.imp
}
}
