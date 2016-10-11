library(nimble)
library(parallel)

crossValCalculate <- function(row, MCMCOut, dataDimensions, saveData){
  simDataArray <- array(MCMCOut[row,], dim = c(dataDimensions))
  discrep <- sum((simDataArray - saveData)^2)
  return(discrep)
}

crossValidateOne <- function(model,
                             dataNames,
                             MCMCIter,
                             burnIn, thin,
                             leaveOutIndex,
                             MCMCdefs=NULL){
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
  newModel <- model$newModel()
  compileNimble(newModel)
  numBlocks <- newModel$getVarInfo(dataNames)[['maxs']][leaveOutIndex]
  dataDimensions <- newModel$getVarInfo(dataNames)[['maxs']]
  saveData <- array(nimble:::values(newModel, dataNames),
                    dim = c(dataDimensions))

  calcCrossVal <- function(i,
                           saveData,
                           newModel,
                           leaveOutIndex,
                           dataDimensions,
                           dataNames){
    tempData <- saveData
    compileNimble(newModel)
    newModel$resetData()
    evalCode1 <- paste0("tempData[",paste(rep(",", leaveOutIndex - 1),
    collapse=""), i, paste(rep(",", length(dataDimensions) -
                                   leaveOutIndex), collapse=""),"] <- NA")
    eval(parse(text = evalCode1))
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
    crossValValue <- unlist(mclapply(ceiling(burnIn/thin):sampNum,
                                     crossValCalculate, MCMCout,
                                     dataDimensions, saveData))
    crossValAverage <- log(mean(crossValValue))
    nimble:::clearCompiled(C.modelMCMC)
    return(crossValAverage)
  }
  crossVal <- sum(sapply(1:numBlocks, calcCrossVal,
                           saveData,
                           newModel,
                           leaveOutIndex,
                           dataDimensions,
                           dataNames))

  return(crossVal)
}
