write.mitmlSAV <- function(x, filename){
# write to native SPSS format

  if(!inherits(x, "mitml") && !inherits(x, "mitml.list")) stop("'x' must be of class 'mitml' or 'mitml.list'.")
  if(!grepl(".sav$", tolower(filename))) filename <- paste(filename, ".sav", sep = "")

  # convert mitml to mitml.list
  if(inherits(x, "mitml")){
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

