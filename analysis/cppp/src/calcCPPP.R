library(parallel)
## draws a random number from the posterior, reassigned it as the
## paramter value, simulates data from the model

discFuncGenerator <- nimbleFunction(
  setup = function(model){},
  run = function(){
    output <- calculate(model)
    returnType(double(0))
    return(output)
  }
)

discFuncGenerator <- nimbleFunction(
  setup = function(model, ...){
    params <- list(...)
    dataNames <- params[[1]]
  },
  run = function(){
    output <- max(values(model,dataNames))
    returnType(double(0))
    return(output)
  }
)


pppFunc <- nimbleFunction(
  setup = function(model, dataNames, paramNames,
                   mcmcMV, MCMCIter, thin, burnIn, averageParams, ...){
    paramDependencies <- model$getDependencies(paramNames)
    discFunc <- discFuncGenerator(model, ...)
  },
  run = function(N = integer(0)){
    output <- numeric(N)
    if(averageParams == 0){
      discMean <- discFunc()
    }
    else{
      mcmcSamps <-floor(MCMCIter/thin - burnIn)
      discMean <- 0
      for(i in 1:mcmcSamps){
        copy(mcmcMV, model, paramNames, paramNames, row = i + burnIn)
        discMean <- discMean + discFunc()
      }
      discMean <- discMean / mcmcSamps
    }
    if(is.nan(discMean)) return(NA)
    for(i in 1:N){
      randNum <- ceiling(runif(1, 0, (MCMCIter)/thin - burnIn - 1 ))
      nimCopy(mcmcMV, model, paramNames, row = burnIn + randNum)
      calculate(model, paramDependencies)
      simulate(model, dataNames, includeData = TRUE)
      deviance <- discFunc()
      if(deviance >= discMean) output[i] <- 1
      else output[i] <- 0
    }
    out <- mean(output)
    returnType(double(0))
    return(out)
  })




## calculate proportion of deviances that are greater or equal to the
## observed value
calcCPPP <- function(MCMCIter,
                     NSamp,
                     C.pppFunc,
                     cppp.C.mcmc,
                     firstRun){
  if(firstRun  == 0)
    cppp.C.mcmc$run(MCMCIter)
  
  pre.pp <- C.pppFunc$run(NSamp)
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
                          thin,
                          discFuncGenerator, 
                          averageParams, 
                          ...){ ## thinning used in original mcmc run

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
                          burnIn,
                          averageParams,
                          ...)


  # modelAverageFunc <- discAverageFunc(R.model,
  #                                     paramNames,
  #                                     mcmcMV,
  #                                     MCMCIter,
  #                                     thin,
  #                                     burnIn)

  C.pppFunc <- compileNimble(modelpppFunc,
                             project = R.model)
  # C.averageFunc <- compileNimble(modelAverageFunc,
  #                                project = R.model)

  ## calculate deviances
  obs.cppp <- calcCPPP(MCMCIter,
                       NSamp,
                       C.pppFunc,
                       orig.C.mcmc,
                       firstRun = 1)
  print(obs.cppp)

  ## refits model with sampled data, reruns, enter inner loop,
  ## calculates distbution of PPPs
  simPppDist <- function(iteration){
    print(iteration)
    simulate(orig.C.model,  includeData =  TRUE)
    out <- calcCPPP(MCMCIter,
                    NSamp,
                    C.pppFunc,
                    orig.C.mcmc,
                    firstRun = 0)
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

