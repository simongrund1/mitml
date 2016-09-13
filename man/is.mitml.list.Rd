\name{is.mitml.list}
\alias{is.mitml.list}
\title{Check if an object is of class \code{mitml.list}}
\description{This function checks if its argument is a list of class \code{mitml.list}.}
\usage{

is.mitml.list(x)

}
\arguments{

  \item{x}{An R object.}

}
\value{

Either \code{TRUE} or \code{FALSE}.
A warning message is displayed if the contents of \code{x} do not appear to be data frames.

}
\author{Simon Grund}
\seealso{\code{\link{as.mitml.list}}}

\examples{
l <- list(data.frame(x=rnorm(20)))
l <- as.mitml.list(l)
is.mitml.list(l)
# TRUE

l <- as.list(1:10)
is.mitml.list(l)
# FALSE

class(l) <- "mitml.list"
is.mitml.list(l)
# TRUE, with a warning
}
