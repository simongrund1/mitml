# mitml
#### Tools for multiple imputation in multilevel modeling
---

This [R](https://www.r-project.org/) package provides tools for multiple imputation of missing data in multilevel modeling. It includes a user-friendly interface to the packages `pan` and `jomo`, and several functions for visualization, data management, and the analysis of multiply imputed data sets.

The purpose of `mitml` is to provide users with a set of effective and user-friendly tools for multiple imputation of multilevel data without requiring advanced knowledge of its statistical underpinnings. Examples and additional information can be found in the official [documentation](https://cran.r-project.org/web/packages/mitml/mitml.pdf) of the package. If you use `mitml` and have suggestions for improvement, please email me (see [here](https://cran.r-project.org/web/packages/mitml/)) or file an [issue](https://github.com/simongrund1/mitml/issues) at the GitHub repository.

#### CRAN version

The official version of `mitml` is hosted on CRAN and may be found [here](https://cran.r-project.org/web/packages/mitml/). The CRAN version can be install from within R using:

```R
install.packages("mitml")
```

#### GitHub version
The version hosted here is essentially a snapshot of the CRAN version, allowing better tracking of [issues](https://github.com/simongrund1/mitml/issues) and requests for future release. In general, however, it may be different from the CRAN version and might contain intended changes in advance. The GitHub version can be installed using `devtools` as:

```R
install.packages("devtools")
devtools::install_github("simongrund1/mitml")
```
