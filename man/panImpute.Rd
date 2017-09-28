\name{panImpute}
\alias{panImpute}
\title{Impute multilevel missing data using \code{pan}}
\description{This function provides an interface to the \code{pan} package for multiple imputation of multilevel data (Schafer & Yucel, 2002).
Imputations can be generated using \code{type} or \code{formula}, which offer different options for model specification.}

\usage{

panImpute(data, type, formula, n.burn=5000, n.iter=100, m=10, group=NULL, 
  prior=NULL, seed=NULL, save.pred=FALSE, keep.chains=c("full","diagonal"),
  silent=FALSE)

}

\arguments{

  \item{data}{A data frame containing incomplete and auxiliary variables, the cluster indicator variable, and any other variables that should be present in the imputed datasets.}
  \item{type}{An integer vector specifying the role of each variable in the imputation model (see details).}
  \item{formula}{A formula specifying the role of each variable in the imputation model. The basic model is constructed by \code{model.matrix}, thus allowing to include derived variables in the imputation model using \code{I()} (see details and examples).}
  \item{n.burn}{The number of burn-in iterations before any imputations are drawn. Default is to 5,000.}
  \item{n.iter}{The number of iterations between imputations. Default is to 100.}
  \item{m}{The number of imputed data sets to generate.}
  \item{group}{(optional) A character string denoting the name of an additional grouping variable to be used with the \code{formula} argument. When specified, the imputation model is run separately within each of these groups.}
  \item{prior}{(optional) A list with components \code{a}, \code{Binv}, \code{c}, and \code{Dinv} for specifying prior distributions for the covariance matrix of random effects and the covariance matrix of residuals (see details). Default is to using least-informative priors.}
  \item{seed}{(optional) An integer value initializing \code{pan}'s random number generator for reproducible results. Default is to using random seeds.}
  \item{save.pred}{(optional) Logical flag indicating if variables derived using \code{formula} should be included in the imputed data sets. Default is to \code{FALSE}.}
  \item{keep.chains}{(optional) A character string denoting which parameter chains to save. Default is to save all chains (see details).}
  \item{silent}{(optional) Logical flag indicating if console output should be suppressed. Default is to \code{FALSE}.}

}

\details{

This function serves as an interface to the \code{pan} algorithm.
The imputation model can be specified using either the \code{type} or the \code{formula} argument.

The \code{type} interface is designed to provide quick-and-easy imputations using \code{pan}.
The \code{type} argument must be an integer vector denoting the role of each variable in the imputation model:
\itemize{
  \item{\code{1}: target variables containing missing data}
  \item{\code{2}: predictors with fixed effect on all targets (completely observed)}
  \item{\code{3}: predictors with random effect on all targets (completely observed)}
  \item{\code{-1}: grouping variable within which the imputation is run separately}
  \item{\code{-2}: cluster indicator variable}
  \item{\code{0}: variables not featured in the model}
}
At least one target variable and the cluster indicator must be specified.
The intercept is automatically included both as a fixed and random effect.
If a variable of type \code{-1} is found, then separate imputations are performed within each level of that variable.

The \code{formula} argument is intended as more flexible and feature-rich interface to \code{pan}. Specifying the \code{formula} argument is similar to specifying other formulae in R.
Given below is a list of operators that \code{panImpute} currently understands:
\itemize{
  \item{\code{~}: separates the target (left-hand) and predictor (right-hand) side of the model}
  \item{\code{+}: adds target or predictor variables to the model}
  \item{\code{*}: adds an interaction term of two or more predictors}
  \item{\code{|}: denotes cluster-specific random effects and specifies the cluster indicator (e.g., \code{1|ID})}
  \item{\code{I()}: defines functions to be interpreted by \code{model.matrix}}
}
Predictors are allowed to have fixed effects, random effects, or both on all target variables.
The intercept is automatically included both as a fixed and a random effect, but it can be constrained if necessary (see examples).
Note that, when specifying random effects other than the intercept, these will \emph{not} be automatically added as fixed effects and must be included explicitly.
Any predictors defined by \code{I()} will be used for imputation but not included in the data set unless \code{save.pred=TRUE}.

In order to run separate imputations for each level of an additional grouping variable, the \code{group} argument may be used.
The name of the grouping variable must be given in quotes.

As a default prior, \code{panImpute} uses "least informative" inverse-Wishart priors for the covariance matrix of random effects and the covariance matrix of residuals, that is, with minimum degrees of freedom (largest dispersion) and identity matrices for scale.
For better control, the \code{prior} argument may be used for specifying alternative prior distributions.
These must be supplied as a list containing the following components:
\itemize{
  \item{\code{a}: degrees of freedom for the covariance matrix of residuals}
  \item{\code{Binv}: scale matrix for the covariance matrix of residuals}
  \item{\code{c}: degrees of freedom for the covariance matrix of random effects}
  \item{\code{Dinv}: scale matrix for the covariance matrix of random effects}
}
A sensible choice for a diffuse non-default prior is to set the degrees of freedom to the lowest value possible, and the scale matrices according to a prior guess of the corresponding covariance matrices (see Schafer & Yucel, 2002).

In imputation models with many parameters, the number of parameter chains being saved can be reduced with the \code{keep.chains} argument.
If set to \code{full} (the default), all chains are saved.
If set to \code{diagonal}, only chains pertaining to fixed effects and the diagonal entries of the covariance matrices are saved.
This setting influences the storage mode of parameters (e.g., dimensions and indices of arrays) and should be used with caution.

}

\value{

Returns an object of class \code{mitml}, containing the following components:

  \item{data}{The original (incomplete) data set, sorted according to the cluster variable and (if given) the grouping variable, with several attributes describing the original row order (\code{"sort"}) and grouping (\code{"group"}.}
  \item{replacement.mat}{A matrix containing the multiple replacements (i.e., imputations) for each missing value. The replacement matrix contains one row for each missing value and one one column for each imputed data set.}
  \item{index.mat}{A matrix containing the row and column index for each missing value. The index matrix is used to \emph{link} the missing values in the data set with their corresponding rows in the replacement matrix.}
  \item{call}{The matched function call.}
  \item{model}{A list containing the names of the cluster variable, the target variables, and the predictor variables with fixed and random effects, respectively.}
  \item{random.L1}{A character string denoting the handling of random residual covariance matrices (not used here; see \code{jomoImpute}).}
  \item{prior}{The prior parameters used in the imputation model.}
  \item{iter}{A list containing the number of burn-in iterations, the number of iterations between imputations, and the number of imputed data sets.}
  \item{par.burnin}{A multi-dimensional array containing the parameters of the imputation model from the burn-in phase.}
  \item{par.imputation}{A multi-dimensional array containing the parameters of the imputation model from the imputation phase.}

}

\note{

For objects of class \code{mitml}, methods for the generic functions \code{print}, \code{summary}, and \code{plot} have been defined.
\code{mitmlComplete} is used for extracting the imputed data sets.

}

\references{
Schafer, J. L., and Yucel, R. M. (2002). Computational strategies for multivariate linear mixed-effects models with missing values. \emph{Journal of Computational and Graphical Statistics, 11}, 437-457.
}

\author{Simon Grund, Alexander Robitzsch, Oliver Luedtke}
\seealso{\code{\link{jomoImpute}}, \code{\link{mitmlComplete}}, \code{\link{summary.mitml}}, \code{\link{plot.mitml}}}

\examples{
# NOTE: The number of iterations in these examples is much lower than it
# should be! This is done in order to comply with CRAN policies, and more
# iterations are recommended for applications in practice!

data(studentratings)

# *** ................................
# the 'type' interface
# 

# * Example 1.1: 'ReadDis' and 'SES', predicted by 'ReadAchiev' and 
# 'CognAbility', with random slope for 'ReadAchiev'

type <- c(-2,0,0,0,0,0,3,1,2,0)
names(type) <- colnames(studentratings)
type

imp <- panImpute(studentratings, type=type, n.burn=1000, n.iter=100, m=5)

# * Example 1.2: 'ReadDis' and 'SES' groupwise for 'FedState',
# and predicted by 'ReadAchiev'

type <- c(-2,-1,0,0,0,0,2,1,0,0)
names(type) <- colnames(studentratings)
type

imp <- panImpute(studentratings, type=type, n.burn=1000, n.iter=100, m=5)

# *** ................................
# the 'formula' interface
# 

# * Example 2.1: imputation of 'ReadDis', predicted by 'ReadAchiev'
# (random intercept)

fml <- ReadDis ~ ReadAchiev + (1|ID)
imp <- panImpute(studentratings, formula=fml, n.burn=1000, n.iter=100, m=5)

# ... the intercept can be suppressed using '0' or '-1' (here for fixed intercept)
fml <- ReadDis ~ 0 + ReadAchiev + (1|ID)
imp <- panImpute(studentratings, formula=fml, n.burn=1000, n.iter=100, m=5)

# * Example 2.2: imputation of 'ReadDis', predicted by 'ReadAchiev'
# (random slope)

fml <- ReadDis ~ ReadAchiev + (1+ReadAchiev|ID)
imp <- panImpute(studentratings, formula=fml, n.burn=1000, n.iter=100, m=5)

# * Example 2.3: imputation of 'ReadDis', predicted by 'ReadAchiev',
# groupwise for 'FedState'

fml <- ReadDis ~ ReadAchiev + (1|ID)
imp <- panImpute(studentratings, formula=fml, group="FedState", n.burn=1000,
n.iter=100, m=5)

# * Example 2.4: imputation of 'ReadDis', predicted by 'ReadAchiev'
# including the cluster mean of 'ReadAchiev' as an additional predictor

fml <- ReadDis ~ ReadAchiev + I(clusterMeans(ReadAchiev,ID)) + (1|ID)
imp <- panImpute(studentratings, formula=fml, n.burn=1000, n.iter=100, m=5)

# ... using 'save.pred' to save the calculated cluster means in the data set
fml <- ReadDis ~ ReadAchiev + I(clusterMeans(ReadAchiev,ID)) + (1|ID)
imp <- panImpute(studentratings, formula=fml, n.burn=1000, n.iter=100, m=5,
save.pred=TRUE)

head(mitmlComplete(imp,1))
}
\keyword{models}
