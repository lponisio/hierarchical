library(parallel)
## draws a random number from the posterior, reassigned it as the
## paramter value, simulates data from the model


## calculate proportion of deviances that are greater or equal to
## the observed value
calcCPPP <- nimbleFunction(
  setup = function(model, dataNames, paramNames, MCMCIter, thin){
    mcmcSetup <- configureMCMC(model, thin = thin, monitors = paramNames)
    mcmcBuild <- buildMCMC(mcmcSetup)
    mcmcMV <- mcmcBuild$mvSamples
    paramDependencies <- model$getDependencies(paramNames)
  },
  run = function(NSamp = integer(0)){
     ppInd = numeric(NSamp)
    # 
    mcmcBuild$run(MCMCIter)
    observedDisc <- calculate(model)
     for(i in 1:NSamp){
      randNum <- ceiling(runif(1, 0, MCMCIter/thin -1 ))
      nimCopy(mcmcMV, model, paramNames, row = randNum)
      calculate(model, paramDependencies)
      simulate(model, dataNames, includeData = TRUE)
      otherDisc <- calculate(model)
      if(otherDisc >= observedDisc)
        ppInd[i] <- 1
      else
         ppInd[i] <- 0
     }
    prePP <- mean(ppInd)
    returnType(double(0))
    return(prePP)    
  })


nimGenerateCPPP <-  nimbleFunction(
  setup = function(model, dataNames,
                          paramNames, ## vector of parameters to monitor
                          MCMCIter, ## number of samples
                          NSamp, thin){
    my_initializeModel <- initializeModel(model)
    
    cpppFunc <- calcCPPP(model, dataNames, paramNames, MCMCIter, thin)
  },
  
  run = function(NPDist = integer(0)){
    cpppInds = numeric(NPDist)
   
    my_initializeModel$run()
    obsCppp <- cpppFunc$run(NSamp)
    for(i in 1:NPDist){
      simulate(model,  includeData =  TRUE)
      calculate(model)
       sampleCppp <- cpppFunc$run(NSamp)
       if(sampleCppp <= obsCppp)
         cpppInds[i] <- 1
       else
        cpppInds[i] <- 0
    }
    cpppOut <- mean(cpppInds)
    returnType(double(0))
    return(cpppOut)
  }
)                        


