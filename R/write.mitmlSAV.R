write.mitmlSAV <- function(x, filename){
# write to native SPSS format

  if(!("mitml" %in% class(x)) & !("mitml.list" %in% class(x))) stop("'x' must be of class 'mitml' or 'mitml.list'.")
  if(!grepl(".sav$", tolower(filename))) filename <- paste(filename, ".sav", sep = "")

  # convert mitml to mitml.list
  if("mitml" %in% class(x)){
    x <- mitmlComplete(x, "all", force.list = TRUE)
  }

  # add imputation indicator
  for(ii in 1:length(x)){
    x[[ii]] <- cbind(ii-1, x[[ii]])
    colnames(x[[ii]])[1] <- "Imputation_"
  }

  # write to file
  out <- do.call(rbind, x)
  haven::write_sav(out, filename)

  invisible()

}

