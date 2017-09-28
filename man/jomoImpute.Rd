\name{jomoImpute}
\alias{jomoImpute}
\title{Impute multilevel missing data using \code{jomo}}
\description{This function provides an interface to the \code{jomo} package, which uses the MCMC algorithms presented in Carpenter & Kenward (2013).
Through this wrapper function, imputations can be generated for (mixed) categorical and continuous variables (Goldstein et al., 2009) at level 1 and level 2 as well as imputation using random (residual) covariance matrices (Yucel, 2011).
Imputations can be generated using \code{type} or \code{formula}, which offer different options for model specification.}

\usage{

jomoImpute(data, type, formula, random.L1=c("none","mean","full"), n.burn=5000,
  n.iter=100, m=10, group=NULL, prior=NULL, seed=NULL, save.pred=FALSE,
  keep.chains=c("full","diagonal"), silent=FALSE)

}

\arguments{

  \item{data}{A data frame containing incomplete and auxiliary variables, the cluster indicator variable, and any other variables that should be present in the imputed datasets.}
  \item{type}{An integer vector specifying the role of each variable in the imputation model or a list of two vectors specifying a two-level model (see details).}
  \item{formula}{A formula specifying the role of each variable in the imputation model or a list of two formulas specifying a two-level model. The basic model is constructed by \code{model.matrix}, thus allowing to include derived variables in the imputation model using \code{I()} (see details and examples).}
  \item{random.L1}{A character string denoting if the covariance matrix of residuals should vary across groups and how the values of these matrices are stored (see details). Default is to \code{"none"}, assuming a common covariance matrix across clusters.}
  \item{n.burn}{The number of burn-in iterations before any imputations are drawn. Default is to 5,000.}
  \item{n.iter}{The number of iterations between imputations. Default is to 100.}
  \item{m}{The number of imputed data sets to generate. Default is to 10.}
  \item{group}{(optional) A character string denoting the name of an additional grouping variable to be used with the \code{formula} argument. When specified, the imputation model is run separately within each of these groups.}
  \item{prior}{(optional) A list with components \code{Binv}, \code{Dinv}, and \code{a} for specifying prior distributions for the covariance matrix of random effects and the covariance matrix of residuals (see details). Default is to using least-informative priors.}
  \item{seed}{(optional) An integer value initializing R's random number generator for reproducible results. Default is to using the global seed.}
  \item{save.pred}{(optional) Logical flag indicating if variables derived using \code{formula} should be included in the imputed data sets. Default is to \code{FALSE}.}
  \item{keep.chains}{(optional) A character string denoting which parameter chains to save. Default is to save all chains (see details).}
  \item{silent}{(optional) Logical flag indicating if console output should be suppressed. Default is to \code{FALSE}.}

}

\details{

This function serves as an interface to the \code{jomo} package and supports imputation of multilevel continuous and categorical data.
In order for categorical variables to be detected correctly, these must be formatted as a \code{factor} variables (see examples).
The imputation model can be specified using either the \code{type} or the \code{formula} argument.

The \code{type} interface is designed to provide quick-and-easy imputations using \code{jomo}. The \code{type} argument must be an integer vector denoting the role of each variable in the imputation model:
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

The \code{formula} argument is intended as more flexible and feature-rich interface to \code{jomo}. Specifying the \code{formula} argument is similar to specifying other formulae in R.
Given below is a list of operators that \code{jomoImpute} currently understands:
\itemize{
  \item{\code{~}: separates the target (left-hand) and predictor (right-hand) side of the model}
  \item{\code{+}: adds target or predictor variables to the model}
  \item{\code{*}: adds an interaction term of two or more predictors}
  \item{\code{|}: denotes cluster-specific random effects and specifies the cluster indicator (e.g., \code{1|ID})}
  \item{\code{I()}: defines functions to be interpreted by \code{model.matrix}}
}
Predictors are allowed to have fixed effects, random effects, or both on all target variables.
The intercept is automatically included both as a fixed and a random effect, but it can be constrained if necessary (see \code{\link{panImpute}}).
Note that, when specifying random effects other than the intercept, these will \emph{not} be automatically added as fixed effects and must be included explicitly.
Any predictors defined by \code{I()} will be used for imputation but not included in the data set unless \code{save.pred=TRUE}.

If missing data occur at both levels of the sample (level 1 and level 2), then a list of two \code{formula}s or \code{type}s may be provided.
The first element of this list denotes the imputation model for variables at level 1.
The second element denotes the imputation model for variables at level 2.
In such a case, missing values are imputed jointly at both levels (see examples, see also Carpenter and Kenward, 2013; Goldstein et al., 2009).

It is possible to model the covariance matrix of residuals at level 1 as random across clusters (Yucel, 2011; Carpenter & Kenward, 2013).
The \code{random.L1} argument determines this behavior and how the values of these matrices are stored.
If set to \code{"none"}, a common covariance matrix is assumed across groups (similar to \code{panImpute}).
If set to \code{"mean"}, the covariance matrices are random, but only the average covariance matrix is stored at each iteration.
If set to \code{"full"}, the covariance matrices are random, and all variances and covariances from all clusters are stored.

In order to run separate imputations for each level of an additional grouping variable, the \code{group} argument may be used.
The name of the grouping variable must be given in quotes.

As a default prior, \code{jomoImpute} uses "least informative" inverse-Wishart priors for the covariance matrix of random effects (and residuals at level 2) and the covariance matrix of residuals at level 1, that is, with minimum degrees of freedom (largest dispersion) and identity matrices for scale.
For better control, the \code{prior} argument may be used for specifying alternative prior distributions.
These must be supplied as a list containing the following components:
\itemize{
  \item{\code{Binv}: scale matrix for the covariance matrix of residuals at level 1}
  \item{\code{Dinv}: scale matrix for the covariance matrix of random effects and residuals at level 2}
  \item{\code{a}: starting value for the degrees of freedom of random covariance matrices of residuals (only used with \code{random.L1="mean"} or \code{random.L1="full"})}
}
Note that \code{jomo} does not allow for the degrees of freedom for the inverse-Wishart prior to be specified by the user.
These are always set to the lowest value possible (largest dispersion) or determined iteratively if the residuals at level 1 are modeled as random (see above).

In imputation models with many parameters, the number of parameter chains being saved can be reduced with the \code{keep.chains} argument (see \code{\link{panImpute}}).
This setting influences the storage mode of parameters (e.g., dimensions and indices of arrays) and should be used with caution.

}

\value{

Returns an object of class \code{mitml}, containing the following components:

  \item{data}{The original (incomplete) data set, sorted according to the cluster variable and (if given) the grouping variable, with several attributes describing the original order (\code{"sort"}), grouping (\code{"group"}) and factor levels of categorical variables.}
  \item{replacement.mat}{A matrix containing the multiple replacements (i.e., imputations) for each missing value. The replacement matrix contains one row for each missing value and one one column for each imputed data set.}
  \item{index.mat}{A matrix containing the row and column index for each missing value. The index matrix is used to \emph{link} the missing values in the data set with their corresponding rows in the replacement matrix.}
  \item{call}{The matched function call.}
  \item{model}{A list containing the names of the cluster variable, the target variables, and the predictor variables with fixed and random effects, at level 1 and level 2, respectively.}
  \item{random.L1}{A character string denoting the handling of the (random) covariance matrix of residuals at level 1 (see details).}
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
Carpenter, J. R., & Kenward, M. G. (2013). \emph{Multiple imputation and its application}. Hoboken, NJ: Wiley.

Goldstein, H., Carpenter, J., Kenward, M. G., & Levin, K. A. (2009). Multilevel models with multivariate mixed response types. \emph{Statistical Modelling}, 9, 173-197.

Yucel, R. M. (2011). Random covariances and mixed-effects models for imputing multivariate multilevel continuous data. \emph{Statistical Modelling}, 11, 351-370.
}

\author{Simon Grund, Alexander Robitzsch, Oliver Luedtke}
\seealso{\code{\link{panImpute}}, \code{\link{mitmlComplete}}, \code{\link{summary.mitml}}, \code{\link{plot.mitml}}}

\examples{
# NOTE: The number of iterations in these examples is much lower than it
# should be! This is done in order to comply with CRAN policies, and more
# iterations are recommended for applications in practice!

data(studentratings)
data(leadership)

# ***
# for further examples, see "panImpute"
#

?panImpute

# *** ................................
# the 'type' interface
# 

# * Example 1.1 (studentratings): 'ReadDis' and 'SES', predicted by 'ReadAchiev'
# (random slope)

type <- c(-2,0,0,0,0,1,3,1,0,0)
names(type) <- colnames(studentratings)
type

imp <- jomoImpute(studentratings, type=type, n.burn=100, n.iter=10, m=5)

# * Example 1.2 (leadership): all variables (mixed continuous and categorical
# data with missing values at level 1 and level 2)

type.L1 <- c(-2,1,0,1,1)   # imputation model at level 1
type.L2 <- c(-2,0,1,0,0)   # imputation model at level 2
names(type.L1) <- names(type.L2) <- colnames(leadership)

type <- list(type.L1, type.L2)
type

imp <- jomoImpute(leadership, type=type, n.burn=100, n.iter=10, m=5)


# *** ................................
# the 'formula' interface
# 

# * Example 2.1 (studentratings): 'ReadDis' and 'SES' predicted by 'ReadAchiev'
# (random slope)

fml <- ReadDis + SES ~ ReadAchiev + (1|ID)
imp <- jomoImpute(studentratings, formula=fml, n.burn=100, n.iter=10, m=5)

# * Example 2.2 (studentratings): 'ReadDis' predicted by 'ReadAchiev' and the
# the cluster mean of 'ReadAchiev'

fml <- ReadDis ~ ReadAchiev + I(clusterMeans(ReadAchiev,ID)) + (1|ID)
imp <- jomoImpute(studentratings, formula=fml, n.burn=100, n.iter=10, m=5)

# * Example 2.3 (studentratings): 'ReadDis' predicted by 'ReadAchiev', groupwise
# for 'FedState'

fml <- ReadDis ~ ReadAchiev + (1|ID)
imp <- jomoImpute(studentratings, formula=fml, group="FedState", n.burn=100, n.iter=10, m=5)

# * Example 2.4 (leadership): all variables (mixed continuous and categorical
# data with missing values at level 1 and level 2)

fml <- list( JOBSAT + NEGLEAD + WLOAD ~ 1 + (1|GRPID) , COHES ~ 1 )
imp <- jomoImpute(leadership, formula=fml, n.burn=100, n.iter=10, m=5)
}
\keyword{models}
