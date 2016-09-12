crossValidate <- function(bugsModelList, dataValuesList, initialValuesList, dimsList){}




crossValidateOne <- function(model, leaveOutIndex){
  ## fill in each element of data along leaveOutIndex with NAs.  
  ## then, data na values will be filled in as each mcmc runs.  These estimated data values can be compared to 
  ## known data values, and the average loss (0/1) can be taken over all MCMC runs.  then take the average of these for all data points? woo!

  
  
}