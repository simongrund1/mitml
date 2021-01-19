cbind.mitml.list <- function(...){
# merges two objects of class "mitml.list" by appending columns of list entries

  Map(cbind, ...)

}
