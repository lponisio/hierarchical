library(nimble)
library(coda)
library(parallel)

crossValCalculate <- function(row, MCMCOut, dataDimensions, saveData){
  simDataArray <- array(MCMCOut[row,], dim = c(dataDimensions))
  discrep <- sum((simDataArray - saveData)^2)
  return(discrep)
}

crossValidateOne <- function(model,
                             dataNames,
                             MCMCIter,
                             burnInProp, thin,
                             leaveOutIndex,
                             MCMCdefs=NULL,
                             returnChains=TRUE){
  simpSetMCMCDefs <- function(Rmodel, MCMCdefs, MCMCname) {
    eval(MCMCdefs[[MCMCname]])
  }
  if(!is.null(MCMCdefs)){
    ## eval(MCMCdefs.opt2[['nimbleOpt2']])
    ## customSpec <- simpSetMCMCDefs(occ.R.model, MCMCdefs.opt2, 'nimbleOpt2')
  }
  ## fill in each element of data along leaveOutIndex with NAs.  then,
  ## data na values will be filled in as each mcmc runs.  These
  ## estimated data values can be compared to known data values, and
  ## the average loss (0/1) can be taken over all MCMC runs.  then
  ## take the average of these for all data points? woo!
  if(!inherits(model, "RmodelBaseClass")){
    stop("model is not an Rmodel")
  }
    if(thin = 0 | thin > MCMCIter){
    stop("thin needs to be an integer > 0 and < MCMCIter")
  }
  if(burnInProp >= 1 | burnInProp < 0){
    stop("burnInProp needs to be between 0 and 1")
  }
  testDataNames <- try(model[[dataNames]], silent=TRUE)
  if(inherits(testDataNames, "try-error")){
    stop(paste("dataNames", dataNames,
                "is not the name of the data in model"))
  } else{
    test2DataNames <- all(model$expandNodeNames(dataNames) %in%
                          model$getNodeNames(dataOnly=TRUE))
    if(test2DataNames == FALSE){
      stop(paste("dataNames", dataNames,
                  "is not the name of the data in model"))
    }
  }
  prepModel <- model$newModel()
  numBlocks <- prepModel$getVarInfo(dataNames)[['maxs']][leaveOutIndex]
  dataDimensions <- prepModel$getVarInfo(dataNames)[['maxs']]
  if(leaveOutIndex > length(dataDimensions)){
    stop("leaveOutIndex is not an index in the data")
  }
  saveData <- array(nimble:::values(prepModel, dataNames),
                    dim = c(dataDimensions))
  calcCrossVal <- function(i,
                           saveData,
                           prepModel,
                           leaveOutIndex,
                           dataDimensions,
                           dataNames){
    newModel <- prepModel$newModel()
    compileNimble(newModel)
    tempData <- saveData
    newModel$resetData()
    evalCode1 <- paste0("tempData[",paste(rep(",", leaveOutIndex - 1),
    collapse=""), i, paste(rep(",", length(dataDimensions) -
                                   leaveOutIndex), collapse=""),"] <- NA")
    eval(parse(text = evalCode1))
    ## create a list containing the data with the appropriate dataname
    evalCode2 <- paste0("modelDataList <- list(", dataNames, "= tempData)")
    eval(parse(text=evalCode2))
    newModel$setData(modelDataList)
    modelMCMCConf <- configureMCMC(newModel,
                                   monitors = dataNames, thin = thin)
    modelMCMC <- buildMCMC(modelMCMCConf)
    
    C.modelMCMC <- compileNimble(modelMCMC,
                                 project = newModel,
                                 resetFunctions = (i != 1))    
    C.modelMCMC$run(MCMCIter)
    MCMCout <- as.matrix(C.modelMCMC$mvSamples)
    sampNum <- dim(MCMCout)[1]
    startIndex <- max(c(ceiling(burnInProp*MCMCIter/thin), 1))
    message(paste("dropping data block", i))
    crossValValue <- sapply(startIndex:sampNum,
                                     crossValCalculate, MCMCout,
                                     dataDimensions, saveData)
    crossValAverage <- log(mean(crossValValue))
    nimble:::clearCompiled(C.modelMCMC)
    return(list(crossValAverage =crossValAverage,
                samples= MCMCout))
  }
  crossValOut <- mclapply(1:numBlocks, calcCrossVal,
                           saveData,
                           prepModel,
                           leaveOutIndex,
                           dataDimensions,
                           dataNames)
  crossVal <- sum(sapply(crossValOut, function(x) x$crossValAverage),
                  na.rm=TRUE)
  samples <- lapply(crossValOut, function(x) x$samples)
  chain.diag <-  do.call(cbind,
                         lapply(lapply(samples, geweke.diag),
                                function(x) x$z))
  out <- list(crossVal=crossVal,
              chain.diag=chain.diag,
              samples=samples)
    if(!returnChains){
    out$samples <- NULL
  }
  return(out)
}
