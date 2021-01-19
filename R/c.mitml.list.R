c.mitml.list <- function(...){
# merges two objects of class "mitml.list" by appending list entries

  as.mitml.list(unlist(list(...), recursive = FALSE))

}
