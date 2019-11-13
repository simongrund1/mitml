
n_imputations.mitml <- function(imps) {
  imps$iter$m
}

n_imputations.mitml.list <- function(imps) {
  length(imps)
}

n_imputations.mids <- function(imps) {
  imps$m
}
n_imputations.amelia <- function(imps) {
  length(imps$imputations)
}
n_imputations = function(imps) { UseMethod("n_imputations") }

complete.mitml.list <- function(imps, action = "long") {
  if (action == "long") {
    dplyr::bind_rows(imps, .id = "m")
  } else if (is.numeric(action)) {
    imps[[action]]
  }
}

complete.amelia <- function(imps, action = "long") {
  if (action == "long") {
    dplyr::bind_rows(imps$imputations, .id = "m")
  } else if (is.numeric(action)) {
    imps$imputations[[action]]
  }
}
complete.mids <- mice::complete
complete <- function(...) { UseMethod("complete") }
