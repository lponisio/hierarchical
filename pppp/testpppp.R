library(nimble)
library(parallel)
options(mc.cores=2)

pppFunc <- nimbleFunction(
  setup = function(model, dataName, paramNames, MCMCmv, MCMCIter){
    postLength <- MCMCIter
    paramDependencies <- model$getDependencies(paramNames)
  },
  run = function(N = integer(0)){
    deviances <- numeric(N)
    for(i in 1:N){
      randNum <- ceiling(runif(1, 0, postLength))
      nimCopy(MCMCmv, model, paramNames, row = randNum)
      calculate(model, paramDependencies)
      simulate(model, dataName)
      deviances[i] <- calculate(model)
    }
    return(deviances)
    returnType(double(1))
  }
  )



pumpCode <- nimbleCode({
  for (i in 1:N){
    theta[i] ~ dgamma(alpha,beta)
    lambda[i] <- theta[i]*t[i]
    x[i] ~ dpois(lambda[i])
  }
  alpha ~ dexp(1.0)
  beta ~ dgamma(0.1,1.0)
})

pumpConsts <- list(N = 10,
                   t = c(94.3, 15.7, 62.9, 126, 5.24,
                     31.4, 1.05, 1.05, 2.1, 10.5))
pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22))
pumpInits <- list(alpha = 1, beta = 1,
                  theta = rep(0.1, pumpConsts$N))

generatecPPP <-  function(code, ## model created by nimbleCode()
                          constants,
                          data,
                          dataName,
                          inits, ## inital values for parameters to estimate
                          paramNames, ## vector of parameters to monitor
                          MCMCIter, ## number of samples
                          NSamp,
                          NPDist,
                          thin,...){ ## thinning rate
  calcCPPP <- function(MCMCIter, C.model, NSamp){
    C.mcmc$run(MCMCIter)  
    observedDisc <- calculate(C.model)  
    otherDiscs <- C.pppFunc$run(NSamp)
    pre.pp <- mean(otherDiscs >= observedDisc)
    return(pre.pp)    
  }

  ## build model
  R.model <- nimbleModel(code=code,
                         constants=constants,
                         data=data,
                         inits=inits,
                         check=FALSE,...)
  message('R model created')
  
  ## configure and build mcmc
  mcmc.spec <- configureMCMC(R.model,
                             print=FALSE,
                             monitors = paramNames,
                             thin=thin)
  mcmc <- buildMCMC(mcmc.spec)
  message('MCMC built')
  
  mcmcMV <- mcmc$mvSamples
  modelpppFunc <- pppFunc(R.model, dataName, paramNames, mcmcMV, MCMCIter)
  
  ## compile model in C++
  C.model <- compileNimble(R.model)
  C.mcmc <- compileNimble(mcmc, project = R.model)
  C.pppFunc <- compileNimble(modelpppFunc, project = R.model)
  paramDependencies <- C.model$getDependencies(paramNames)
  message('NIMBLE model compiled')

  ## run model
  print('running model')
  
  obs.cppp <- calcCPPP(MCMCIter, C.model, NSamp)

  simPppDist <- function(interation,
                         MCMCIter,
                         C.mcmc,
                         C.model,
                         paramNames,
                         paramDependencies,
                         NSamp){
    randNum <- ceiling(runif(1, 0, MCMCIter))
    mcmc.samples <- as.matrix(C.mcmc$mvSamples)
    values(C.model, paramNames) <- mcmc.samples[randNum,]
    calculate(C.model, paramDependencies)
    simulate(C.model, dataName)
    out <- calcCPPP(MCMCIter, C.model, NSamp)
    return(out)
  }

  sim.cppp <- unlist(mclapply(1:NPDist, simPppDist,
                     MCMCIter,
                     C.mcmc,
                     C.model,
                     paramNames,
                     paramDependencies,
                     NSamp))
  
  out.cppp <- mean(obs.cppp <= sim.cppp)  
  
  return(list(cppp=out.cppp,
              obs=obs.cppp,
              sim=sim.cppp))
}

generatecPPP(code=pumpCode,
             constants =pumpConsts,
             data = pumpData,
             dataName = 'x',
             inits = pumpInits, 
             paramNames = c('alpha','beta'), 
             MCMCIter = 1000, 
             NSamp = 1000,
             NPDist = 100,
             thin = 2)
