\name{clusterMeans}
\alias{clusterMeans}
\title{Calculate cluster means}
\description{
Calculates the mean of a given variable within each cluster, possibly conditioning on an additional grouping variable.
}

\usage{

clusterMeans(x, cluster, adj=FALSE, group=NULL)

}

\arguments{

  \item{x}{A numeric vector for which cluster means should be calculated. Can also be supplied as a character string denoting a variable in the current environment (see details).}
  \item{cluster}{A numeric vector or a factor denoting the cluster membership of each unit in \code{x}. Can also be supplied as a character string (see details).}
  \item{adj}{Logical flag indicating if person-adjusted group means should be calculated. The cluster mean is then calculated for each unit by excluding that unit from calculating the cluster mean. Default is to \code{FALSE}.}
  \item{group}{(optional) An grouping factor or a variable that can be interpreted as such. If specified, then cluster means are calculated conditionally on the grouping variable, that is, separately within sub-groups. Can also be supplied as a character string (see details).}

}

\details{

This function calculates the mean of a variable within each level of a cluster variable.
Any \code{NA} are omitted during calculation.

The three main arguments of the function can also be supplied as (single) character strings, denoting the name of the respective variables in the current environment.
This is especially useful for calculating several cluster means simultaneously, for example using \code{\link{within.mitml.list}} (see Example 2 below).

}

\value{

Returns a numeric vector with the same length as \code{x} containing the cluster mean for all units.

}

\author{Simon Grund, Alexander Robitzsch}
\seealso{\code{\link{within.mitml.list}}}
\examples{
data(studentratings)

fml <- ReadDis + SES ~ ReadAchiev + (1|ID)
imp <- panImpute(studentratings, formula=fml, n.burn=1000, n.iter=100, m=5)

implist <- mitmlComplete(imp, print=1:5)

# * Example 1: single cluster means

# calculate cluster means (for each data set)
with(implist, clusterMeans(ReadAchiev, ID))

# ... person-adjusted cluster means
with(implist, clusterMeans(ReadAchiev, ID, adj=TRUE))

# ... groupwise cluster means
with(implist, clusterMeans(ReadAchiev, ID, group=Sex))

# * Example 2: automated cluster means using 'for' and 'assign'

# calculate multiple cluster means within multiply imputed data sets
within(implist,{
  vars <- c("ReadAchiev","MathAchiev","CognAbility")
  for(i in vars) assign(paste(i,"Mean",sep="."), clusterMeans(i,ID))
  rm(i,vars)
})
}
