
mutate.mitml.list <- function(.data, ...) {
  r <- lapply(.data, FUN = dplyr::mutate, ...)
  class(r) = class(.data)
  r
}

mutate_at.mitml.list <- function(.data, ...) {
  r <- lapply(.data, FUN = dplyr::mutate_at, ...)
  class(r) = class(.data)
  r
}

mutate_if.mitml.list <- function(.data, ...) {
  r <- lapply(.data, FUN = dplyr::mutate_if, ...)
  class(r) = class(.data)
  r
}

tbl_df.mitml.list <- function(.data, ...) {
  r <- lapply(.data, FUN = dplyr::tbl_df, ...)
  class(r) = class(.data)
  r
}

select.mitml.list <- function(.data, ...) {
  r <- lapply(.data, FUN = dplyr::select, ...)
  class(r) = class(.data)
  r
}

arrange.mitml.list <- function(.data, ...) {
  r <- lapply(.data, FUN = dplyr::arrange, ...)
  class(r) = class(.data)
  r
}

group_by.mitml.list <- function(.data, ...) {
  r <- lapply(.data, FUN = dplyr::group_by, ...)
  class(r) = class(.data)
  r
}


ungroup.mitml.list <- function(.data, ...) {
  r <- lapply(.data, FUN = dplyr::ungroup, ...)
  class(r) = class(.data)
  r
}

summarise.mitml.list <- function(.data, ...) {
  r <- lapply(.data, FUN = dplyr::summarise, ...)
  class(r) = class(.data)
  r
}

