library(parallel)
## draws a random number from the posterior, reassigned it as the
## paramter value, simulates data from the model
pppFunc <- nimbleFunction(
  setup = function(model, dataNames, paramNames,
    mcmcMV, MCMCIter, thin, burnIn){
    paramDependencies <- model$getDependencies(paramNames)
  },
  run = function(N = integer(0)){
    deviances <- numeric(N)
    for(i in 1:N){
      randNum <- ceiling(runif(1, 0, (MCMCIter)/thin - burnIn - 1 ))
      nimCopy(mcmcMV, model, paramNames, row = burnIn + randNum)
      calculate(model, paramDependencies)
      simulate(model, dataNames, includeData = TRUE)
      deviances[i] <- calculate(model)
    }
    return(deviances)
    returnType(double(1))
  })



discAverageFunc <- nimbleFunction(
  setup = function(model, paramNames, mcmcMV, MCMCIter, thin, burnIn){
  },
  run = function(){
    mcmcSamps <-floor(MCMCIter/thin - burnIn)
    discMean <- 0
    for(i in 1:mcmcSamps){
      copy(mcmcMV, model, paramNames, paramNames, row = i + burnIn)
      discMean <- discMean + calculate(model)
    }
    discMean <- discMean / mcmcSamps
    returnType(double(0))
    return(discMean)
  }
  )

## calculate proportion of deviances that are greater or equal to the
## observed value
calcCPPP <- function(MCMCIter,
                     NSamp,
                     C.pppFunc,
                     C.averageFunc,
                     cppp.C.mcmc){
  cppp.C.mcmc$run(MCMCIter)
  observedDisc <- C.averageFunc$run()
  if(!is.finite(observedDisc)) return(NA)
  otherDiscs <- C.pppFunc$run(NSamp)
  pre.pp <- mean(otherDiscs >= observedDisc)
  if(!is.finite(pre.pp))    pre.pp <- NA
  return(pre.pp)    
}


generateCPPP <-  function(R.model,
                          orig.C.model,
                          orig.C.mcmc,
                          orig.mcmc,
                          dataNames, ## names of the data column
                          paramNames, ## vector of parameters to monitor
                          MCMCIter, ## number of mcmc iterations
                          NSamp,## number of samples from posterior
                          NPDist, ## number of simulated PPP values
                          burnInProportion, ## proportion of mcmc to drop
                          thin,...){ ## thinning used in original mcmc run

  burnIn <- ceiling(burnInProportion*(MCMCIter/thin))
  
  origData <- nimble:::values(orig.C.model, dataNames)
  
  ## sample posterior, simulate data from sample 
  paramDependencies <- orig.C.model$getDependencies(paramNames)
  mcmcMV <- orig.mcmc$mvSamples
  modelpppFunc <- pppFunc(R.model,
                          dataNames,
                          paramNames,
                          mcmcMV,
                          MCMCIter,
                          thin,
                          burnIn)

  modelAverageFunc <- discAverageFunc(R.model,
                                      paramNames,
                                      mcmcMV,
                                      MCMCIter,
                                      thin,
                                      burnIn)

  C.pppFunc <- compileNimble(modelpppFunc,
                             project = R.model)
  C.averageFunc <- compileNimble(modelAverageFunc,
                                 project = R.model)

  ## calculate deviances
  obs.cppp <- calcCPPP(MCMCIter,
                       NSamp,
                       C.pppFunc,
                       C.averageFunc,
                       orig.C.mcmc)
  print(obs.cppp)

  ## refits model with sampled data, reruns, enter inner loop,
  ## calculates distbution of PPPs
  simPppDist <- function(iteration){
    print(iteration)
    simulate(orig.C.model,  includeData =  TRUE)
    out <- calcCPPP(MCMCIter,
                    NSamp,
                    C.pppFunc,
                    C.averageFunc,
                    orig.C.mcmc)
    return(out)
  }
  sim.cppp <- unlist(lapply(1:NPDist, simPppDist))

  ## calculates the number of simulated ppp that fall below the obs
  out.cppp <- mean(obs.cppp <= sim.cppp, na.rm = TRUE)  
  
  nimble:::values(orig.C.model, dataNames) <- origData 
  
  return(list(cppp=out.cppp,
              obs=obs.cppp,
              sim=sim.cppp))
}

