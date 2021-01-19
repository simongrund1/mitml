\name{mitml-package}
\alias{mitml-package}
\docType{package}
\title{mitml: Tools for multiple imputation in multilevel modeling}
\description{
Provides tools for multiple imputation of missing data in multilevel modeling.
This package includes a user-friendly interface to the algorithms implemented in the R packages \code{pan} and \code{jomo} as well as several functions for visualizing, managing, and analyzing multiply imputed data sets.

The main interface to \code{pan} is the function \code{\link{panImpute}}, which allows specifying imputation models for continuous variables with missing data at level 1.
In addition, the function \code{\link{jomoImpute}} provides an interface to \code{jomo}, which allows specifying imputation models for both continuous and categorical variables with missing data at levels 1 and 2 as well as single-level data.
The imputations and parameter values are stored in objects of class \code{mitml}.
To obtain the completed (i.e., imputed) data sets, \code{\link{mitmlComplete}} is used, producing a list of imputed data sets of class \code{mitml.list} that can be used in further analyses.

Several additional functions allow for convenient analysis of multiply imputed data sets including (bot not limited to) multilevel analyses with \code{lme4} or \code{nlme} and structural equation modeling with \code{lavaan}.
The functions \code{\link[=with.mitml.list]{within}}, \code{\link[=sort.mitml.list]{sort}}, and \code{\link[=subset.mitml.list]{subset}} can be used to manage and manipulate multiply imputed data sets.
Statistical models are fitted using \code{\link[=with.mitml.list]{with}}.
Pooled parameter estimates can be extracted with \code{\link{testEstimates}}, and model comparisons as well as single- and multi-parameter hypotheses tests can be performed using the functions \code{\link{testModels}} and \code{\link{testConstraints}}.
In addition, the \code{\link[=anova.mitml.result]{anova}} method provides a simple interface to model comparisons.

Data sets can be imported and exported from or to different statistical software packages.
Currently, \code{\link{mids2mitml.list}}, \code{\link{amelia2mitml.list}}, \code{\link{jomo2mitml.list}}, and \code{\link{long2mitml.list}} can be used for importing imputations for other packages in R.
In addition, \code{\link{write.mitmlMplus}}, \code{\link{write.mitmlSAV}}, and \code{\link{write.mitmlSPSS}} export data sets to M\emph{plus} and SPSS, respectively.

Finally, the package provides tools for summarizing and visualizing imputation models, which is useful for the assessment of convergence and the reporting of results.

The data sets contained in this package are published under the same license as the package itself.
They contain simulated data and may be used by anyone free of charge as long as reference to this package is given.
}
\author{
Authors: Simon Grund, Alexander Robitzsch, Oliver Luedtke

Maintainer: Simon Grund <grund@ipn.uni-kiel.de>
}
\keyword{package}
