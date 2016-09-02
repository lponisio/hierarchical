library(parallel)
## draws a random number from the posterior, reassigned it as the
## paramter value, simulates data from the model
pppFunc <- nimbleFunction(
  setup = function(model, dataName, paramNames, MCMCmv, MCMCIter, thin){
    paramDependencies <- model$getDependencies(paramNames)
  },
  run = function(N = integer(0)){
    deviances <- numeric(N)
    for(i in 1:N){
      randNum <- ceiling(runif(1, 0, MCMCIter/thin -1 ))
      nimCopy(MCMCmv, model, paramNames, row = randNum)
      calculate(model, paramDependencies)
      simulate(model, dataName)
      deviances[i] <- calculate(model)
    }
    return(deviances)
    returnType(double(1))
  }
  )



generateCPPP <-  function(R.model,
                          C.model,
                          C.mcmc,
                          mcmc,
                          dataName,
                          paramNames, ## vector of parameters to monitor
                          MCMCIter, ## number of samples
                          NSamp,
                          NPDist,
                          thin,...){
  ## calculate proportion of deviances that are greater or equal to
  ## the observed value
  calcCPPP <- function(MCMCIter, C.model, NSamp,
                       C.pppFunc, C.mcmc){

    C.mcmc$run(MCMCIter)  
    observedDisc <- calculate(C.model)  
    otherDiscs <- C.pppFunc$run(NSamp)
    pre.pp <- mean(otherDiscs >= observedDisc)
    return(pre.pp)    
  }

  ## simulate the distribution of posterior predictive pvalues
  simPppDist <- function(interation,
                         MCMCIter,
                         C.mcmc,
                         C.model,
                         paramNames,
                         paramDependencies,
                         NSamp, thin){
    randNum <- ceiling(runif(1, 0, MCMCIter/thin - 1))
    mcmc.samples <- as.matrix(C.mcmc$mvSamples)
    values(C.model, paramNames) <- mcmc.samples[randNum,]
    calculate(C.model, paramDependencies)
    simulate(C.model, dataName)
    out <- calcCPPP(MCMCIter, C.model, NSamp, C.pppFunc, C.mcmc)
    return(out)
  }

  ## sample posterior, simulate data from sample 
  paramDependencies <- C.model$getDependencies(paramNames)
  mcmcMV <- mcmc$mvSamples
  modelpppFunc <- pppFunc(R.model,
                          dataName,
                          paramNames,
                          mcmcMV,
                          MCMCIter,
                          thin)

  C.pppFunc <- compileNimble(modelpppFunc, project = R.model)

  ## calculate deviances
  obs.cppp <- calcCPPP(MCMCIter, C.model, NSamp, C.pppFunc, C.mcmc)

  ## refits model with sampled data, reruns, enter inner loop,
  ## calculates distbution of PPPs
  sim.cppp <- unlist(mclapply(1:NPDist, simPppDist,
                              MCMCIter,
                              C.mcmc,
                              C.model,
                              paramNames,
                              paramDependencies,
                              NSamp,
                              thin))
  
  ## calculates the number of simulated ppp that fall below the obs
  out.cppp <- mean(obs.cppp <= sim.cppp)  
  
  return(list(cppp=out.cppp,
              obs=obs.cppp,
              sim=sim.cppp))
}
