library(parallel)
library(coda)
## draws a random number from the posterior, reassigned it as the
## paramter value, simulates data from the model

virtualDiscFunction <- nimbleFunctionVirtual(run = function()
                                             returnType(double(0))
                                             )

likeDiscFuncGenerator <- nimbleFunction(
  setup = function(model, ...){},
  run = function(){
    output <- calculate(model)
    returnType(double(0))
    return(output)
  },
  contains = virtualDiscFunction
  )

maxDiscFuncGenerator <- nimbleFunction(
  setup = function(model, ...){
    params <- list(...)
    dataNames <- params[[1]]
  },
  run = function(){
    output <- max(values(model, dataNames))
    returnType(double(0))
    return(output)
  },
  contains = virtualDiscFunction
  )


pppFunc <- nimbleFunction(
  setup = function(model, dataNames, paramNames,
    mcmcMV,
    MCMCIter,
    thin,
    burnIn,
    averageParams,
    discFuncGenerator,
    ...){
    paramDependencies <- model$getDependencies(paramNames)
    
    discFunction <- nimbleFunctionList(virtualDiscFunction)
    discFunction[[1]] <- discFuncGenerator(model, ...) 
  },
  run = function(N = integer(0)){
    output <- numeric(N)
    if(averageParams == 0){
      discMean <- discFunction[[1]]$run()
    }
    else{
      mcmcSamps <-floor(MCMCIter/thin - burnIn)
      discMean <- 0
      for(i in 1:mcmcSamps){
        copy(mcmcMV, model, paramNames, paramNames, row = i + burnIn)
        discMean <- discMean + discFunction[[1]]$run()
      }
      discMean <- discMean / mcmcSamps
    }
    if(is.nan(discMean)) return(NA)
    for(i in 1:N){
      randNum <- ceiling(runif(1, 0, (MCMCIter)/thin - burnIn - 1 ))
      nimCopy(mcmcMV, model, paramNames, row = burnIn + randNum)
      calculate(model, paramDependencies)
      simulate(model, dataNames, includeData = TRUE)
      deviance <- discFunction[[1]]$run()
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
  if(firstRun  == 0){
    cppp.C.mcmc$run(MCMCIter)
    samples <- mcmc(as.matrix(cppp.C.mcmc$mvSamples))
  } else{
    samples <- NA
  }
  
  pre.pp <- C.pppFunc$run(NSamp)
  if(!is.finite(pre.pp))    pre.pp <- NA
  return(list(pre.pp = pre.pp,
              samples = samples))    
}


generateCPPP <-  function(R.model,
                          orig.C.model,
                          orig.C.mcmc,
                          orig.mcmc,
                          dataNames, ## names of the data column
                          paramNames, ## vector of parameters to monitor
                          NSamp,## number of samples from posterior
                          NPDist, ## number of simulated PPP values
                          burnInProp, ## proportion of mcmc to drop
                          discFuncGenerator, 
                          averageParams,
                          returnChains = TRUE,
                          ...){
  if(!inherits(R.model, "RmodelBaseClass")){
    stop("R.model is not an Rmodel")
  }
  if(!inherits(orig.C.model, "CmodelBaseClass")){
    stop("orig.C.model is not an Cmodel")
  }
  if(burnInProp >= 1 | burnInProp < 0){
    stop("burnInProp needs to be between 0 and 1")
  }
  
  thin <- orig.C.mcmc$thin
  MCMCIter <- nrow(as.matrix(orig.C.mcmc$mvSamples))*thin
  
  if(MCMCIter <= 1){
    stop("MCMC must be run on C model before assessment")
  }
  
  burnIn <- ceiling(burnInProp*(MCMCIter/thin))

  if(NSamp > MCMCIter){
    stop("number of samples from posterior must be < number of MCMC iterations")
  }
  
  testDataNames <- try(R.model[[dataNames]], silent=TRUE)
  if(inherits(testDataNames, "try-error")){
    stop(paste("dataNames", dataNames,
               "is not the name of the data in model"))
  } else{
    test2DataNames <- all(R.model$expandNodeNames(dataNames) %in%
                          R.model$getNodeNames(dataOnly=TRUE))
    if(test2DataNames == FALSE){
      stop(paste("dataNames", dataNames,
                 "is not the name of the data in model"))
    }
  }
  testParamNames <- lapply(paramNames, function(x){
    test.this.param <- try(R.model[[x]], silent=TRUE)
    if(inherits(test.this.param, "try-error")){
      stop(paste("paramNames", x,
                 "are not parameters in model"))
    }
  })
  test2ParamNames <- all(R.model$expandNodeNames(paramNames) %in%
                         R.model$getNodeNames(includeData=FALSE,
                                              stochOnly=TRUE))
  if(test2ParamNames == FALSE){
    stop(paste("paramNames", paramNames,
               "are not parameters in model"))
  }


  
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
                          discFuncGenerator = discFuncGenerator,
                          ...)

  C.pppFunc <- compileNimble(modelpppFunc,
                             project = R.model)

  ## calculate deviances
  obs.cppp <- calcCPPP(MCMCIter,
                       NSamp,
                       C.pppFunc,
                       orig.C.mcmc,
                       firstRun = 1)$pre.pp

  ## refits model with sampled data, reruns, enter inner loop,
  ## calculates distbution of PPPs
  simPppDist <- function(iteration){
    message(paste("refitting data iteration", iteration))
    simulate(orig.C.model,  includeData =  TRUE)
    out <- calcCPPP(MCMCIter,
                    NSamp,
                    C.pppFunc,
                    orig.C.mcmc,
                    firstRun = 0)
    return(out)
  }

  sim.cppp <- mclapply(1:NPDist, simPppDist)
  sim.ppp <- sapply(sim.cppp, function(x) x$pre.pp)
  sim.samples <- lapply(sim.cppp, function(x) x$samples)

  ## calculates the number of simulated ppp that fall below the obs
  out.cppp <- mean(obs.cppp <= sim.ppp,
                   na.rm = TRUE)
  chain.diag <- do.call(rbind, lapply(
    lapply(sim.samples, geweke.diag), function(x) x$z))
  
  nimble:::values(orig.C.model, dataNames) <- origData 

  out <- list(cppp=out.cppp,
              obs.ppp=obs.cppp,
              sim.cpp.dist=sim.ppp,
              chain.diagnostics= chain.diag,
              samples=sim.samples)
  if(!returnChains){
    out$samples <- NULL
  }
  return(out)
}

