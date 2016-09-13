library(nimble)

dyesCode <- nimbleCode({
  for (i in 1:BATCHES) {
    for (j in 1:SAMPLES) {
      y[i,j] ~ dnorm(mu[i], sd = sigma.within);
    }
    mu[i] ~ dnorm(theta, sd = sigma.between);
  }
  
  theta ~ dnorm(0.0, 1.0E-10);
  sigma.within ~ dunif(0, 100)
  sigma.between ~ dunif(0, 100)
})

dyesModel <- nimbleModel(dyesCode,
                         constants = list(BATCHES = 6, SAMPLES = 5))

data <- matrix(c(1545, 1540, 1595, 1445, 1595, 1520, 1440, 1555, 1550,
                 1440, 1630, 1455, 1440, 1490, 1605, 1595, 1515, 1450,
                 1520, 1560, 1510, 1465, 1635, 1480, 1580, 1495, 1560,
                 1545, 1625, 1445), nrow = 6)

dyesModel$setData(list(y = data))


## crossValCalculate <- function(row, MCMCOut, dataDimensions, saveData){
##   simDataArray <- array(MCMCOut[row,], dim = c(dataDimensions))
##   evalCode0 <- paste0("simData <- simDataArray[",rep(",", leaveOutIndex - 1), i, paste0(rep(",", length(dataDimensions) - leaveOutIndex)),"]")
##   eval(parse(text = evalCode0))
##   discrep <- simData - saveData
##   return(discrep)
## }



crossValCalculate <- function(row, MCMCOut, dataDimensions, saveData){
  simDataArray <- array(MCMCOut[row,], dim = c(dataDimensions))
  discrep <- sum((simDataArray - saveData)^2)
  return(discrep)
}


crossValidateOne <- function(model, dataNames,
                             MCMCIter,
                             burnIn, thin,
                             leaveOutIndex){
  ## fill in each element of data along leaveOutIndex with NAs.  then,
  ## data na values will be filled in as each mcmc runs.  These
  ## estimated data values can be compared to known data values, and
  ## the average loss (0/1) can be taken over all MCMC runs.  then
  ## take the average of these for all data points? woo!
  newModel <- model$newModel()
  compileNimble(newModel)
  modelMCMCConf <- configureMCMC(newModel, monitors = dataNames, thin = thin)
  modelMCMC <- buildMCMC(modelMCMCConf)
  C.modelMCMC <- compileNimble(modelMCMC,
                               project = newModel)
                               ## resetFunctions = (i != 1))
  
  numBlocks <- newModel$getVarInfo(dataNames)[['maxs']][leaveOutIndex]
  dataDimensions <- newModel$getVarInfo(dataNames)[['maxs']]
  saveData <- array(values(newModel, dataNames), dim = c(dataDimensions))

  calcCrossVal <- function(i){
    tempData <- saveData
    newModel$resetData()
    evalCode1 <- paste0("tempData[",rep(",", leaveOutIndex - 1), i,
                        paste0(rep(",", length(dataDimensions) -
                                   leaveOutIndex)),"] <- NA")
    eval(parse(text = evalCode1))
    evalCode2 <- paste0("modelDataList <- list(", dataNames, "= tempData)")
    eval(parse(text=evalCode2))
    newModel$setData(modelDataList)
    ## modelMCMCConf <- configureMCMC(newModel, monitors = dataNames, thin = thin)
    ## modelMCMC <- buildMCMC(modelMCMCConf)
    ## C.modelMCMC <- compileNimble(modelMCMC,
    ##                              project = newModel,
    ##                              resetFunctions = (i != 1))
    C.modelMCMC$run(MCMCIter)
    MCMCout <- as.matrix(C.modelMCMC$mvSamples)
    sampNum <- dim(MCMCout)[1] 
    crossValValue <- sapply(ceiling(burnIn/thin):sampNum,
                            crossValCalculate, MCMCout,
                            dataDimensions, saveData)
    crossValAverage <- mean(crossValValue)
    return(crossValAverage)
  }
  crossVal <- mean(sapply(1:numBlocks, calcCrossVal))
  ## model$setData(saveData)
  return(crossValAverage)
}
output <- crossValidateOne(dyesModel, "y", 1000, 300, 2, 2)

