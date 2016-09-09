library(parallel)
## draws a random number from the posterior, reassigned it as the
## paramter value, simulates data from the model
pppFunc <- nimbleFunction(
  setup = function(model, dataNames, paramNames, mcmcMV, MCMCIter, thin, burnIn){
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




generateCPPP <-  function(R.model,
                          C.model,
                          C.mcmc,
                          mcmc,
                          dataNames,
                          paramNames, ## vector of parameters to monitor
                          MCMCIter, ## number of samples
                          NSamp,
                          NPDist,
                          burnInProportion,
                          thin,...){
  
  burnIn <- ceiling(burnInProportion*(MCMCIter/thin))
  
  ## calculate proportion of deviances that are greater or equal to
  ## the observed value
  calcCPPP <- function(MCMCIter, C.model, NSamp,
                       C.pppFunc, C.averageFunc,
                       C.mcmc, paramNames){
    
    C.mcmc$run(MCMCIter)
    observedDisc <- C.averageFunc$run()
    if(!is.finite(observedDisc)) return(NA)

    otherDiscs <- C.pppFunc$run(NSamp)
    pre.pp <- mean(otherDiscs >= observedDisc)
    if(!is.finite(pre.pp))    pre.pp <- NA
    return(pre.pp)    
  }

  ## simulate the distribution of posterior predictive pvalues
  simPppDist <- function(interation){
    simulate(C.model,  includeData =  TRUE)
    out <- calcCPPP(MCMCIter, C.model, NSamp, C.pppFunc, C.averageFunc, C.mcmc, paramNames)
    return(out)
  }
  
  
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

  ## sample posterior, simulate data from sample 
  paramDependencies <- C.model$getDependencies(paramNames)
  mcmcMV <- mcmc$mvSamples
  modelpppFunc <- pppFunc(R.model,
                          dataNames,
                          paramNames,
                          mcmcMV,
                          MCMCIter,
                          thin,
                          burnIn)
  
  modelAverageFunc <- discAverageFunc(R.model, paramNames, mcmcMV, MCMCIter, thin, burnIn)

  C.pppFunc <- compileNimble(modelpppFunc, project = R.model)
  C.averageFunc <- compileNimble(modelAverageFunc, project = R.model)

  ## calculate deviances
  obs.cppp <- calcCPPP(MCMCIter, C.model, NSamp, C.pppFunc, C.averageFunc, C.mcmc, paramNames)

  ## refits model with sampled data, reruns, enter inner loop,
  ## calculates distbution of PPPs
  sim.cppp <- unlist(mclapply(1:NPDist, simPppDist))

  ## calculates the number of simulated ppp that fall below the obs
  out.cppp <- mean(obs.cppp <= sim.cppp, na.rm = TRUE)  
  
  return(list(cppp=out.cppp,
              obs=obs.cppp,
              sim=sim.cppp))
}

