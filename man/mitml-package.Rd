\name{mitml-package}
\alias{mitml-package}
\docType{package}
\title{mitml: Tools for multiple imputation in multilevel modeling}
\description{Provides tools for multiple imputation of missing data in multilevel modeling.
This package includes a user-friendly interface to the algorithms implemented in the R packages \code{pan} and \code{jomo}, as well as several functions for visualizing, managing and analyzing multiply imputed data sets.

The main interface to \code{pan} is the function \code{\link{panImpute}}, which allows specification of imputation models for continuous variables with missing data at level 1.
In addition, the function \code{\link{jomoImpute}} provides an interface to \code{jomo}, which extends the functionality of \code{pan} to continuous and categorical variables with missing data at level 1 and level 2.
Imputations and parameter chains are stored in objects of class \code{mitml}.
To obtain the completed (i.e., imputed) data sets, \code{\link{mitmlComplete}} is used, producing a list of imputed data sets of class \code{mitml.list} that can be used in further analyses.

Several additional functions allow for convenient analysis of multiply imputed data sets, especially when using the R packages \code{lme4} and \code{nlme}.
The functions \code{\link[=with.mitml.list]{with}} and \code{\link[=with.mitml.list]{within}} can be used for manipulating the data sets and for model fitting.
Final parameter estimates can be extracted using \code{\link{testEstimates}}.
Single- and multi-parameter hypotheses tests can be performed using the functions \code{\link{testConstraints}} and \code{\link{testModels}}.
In addition, the \code{\link{anova}} method provides a simple interface to model comparisons with automatic refitting of statistical models.

Data sets can be imported and exported from or to different statistical software packages.
Currently, \code{\link{mids2mitml.list}}, \code{\link{jomo2mitml.list}}, and \code{\link{long2mitml.list}} can be used for importing imputations for other packages in R (e.g., \code{mice} and \code{jomo}).
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
