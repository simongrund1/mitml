rbind.mitml.list <- function(...){
# merges two objects of class "mitml.list" by appending rows of list entries

  Map(rbind, ...)

}
