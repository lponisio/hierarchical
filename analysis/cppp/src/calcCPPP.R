library(parallel)
library(coda)

pppFuncVirtual <- nimbleFunctionVirtual(
  run = function(N = integer(0)) returnType(double(0))
  )

virtualDiscFunction <- nimbleFunctionVirtual(
  run = function() returnType(double(0))
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
    setup = function(model,
                     dataNames,
                     paramNames,
                     mcmcMV,
                     MCMCIter,
                     thin,
                     burnIn,
                     averageParams,
                     discFunc,
                     NbootReps,
                     ...){
        blockMV <- modelValues(R.model)
        paramDependencies <- model$getDependencies(paramNames)
        discFunction <- nimbleFunctionList(virtualDiscFunction)
        discFunction[[1]] <- discFunc(model, ...)
    },
    run = function(NpostSamp = integer(0), useBlockMV = logical(0)){
        output <- numeric(NpostSamp)
        if(averageParams == 0){
            discMean <- discFunction[[1]]$run()
        }
        else{
            mcmcSamps <-floor(MCMCIter/thin - burnIn)
            discMean <- 0
            for(i in 1:mcmcSamps){
              if(useBlockMV == TRUE){
                copy(blockMV, model, paramNames, paramNames, row = i)
              }
              else{
                copy(mcmcMV, model, paramNames, paramNames, row = i + burnIn)
              }
              calculate(model, paramDependencies)
              discMean <- discMean + discFunction[[1]]$run()
            }
            discMean <- discMean/mcmcSamps
        }
        if(is.nan(discMean)) return(NA)
        for(i in 1:NpostSamp){
            randNum <- ceiling(runif(1, 0, (MCMCIter)/thin - burnIn - 1 ))
            if(useBlockMV == TRUE){
              nimCopy(blockMV, model, paramNames, row = randNum)
            }
            else{
              nimCopy(mcmcMV, model, paramNames, row = burnIn + randNum)
            }
            calculate(model, paramDependencies)
            simulate(model, dataNames, includeData = TRUE)
            deviance <- discFunction[[1]]$run()
            if(deviance >= discMean) output[i] <- 1
            else output[i] <- 0
        }
        out <- mean(output)
        returnType(double(0))
        return(out)
    },
    methods = list(
        getSD = function(NpostSamp=double(0)){
            blockpppValues <- numeric(NbootReps) ## block estimates of ppp
            l <- ceiling(min(1000, (MCMCIter/thin - burnIn)/20)) ##length of each
            ##block, ensures
            ##it's not too big
            ## total number of blocks available
            q <- (MCMCIter/thin - burnIn) - l + 1
            ##to sample from
            ## number of blocks to use for ppp
            h <- ceiling((MCMCIter/thin - burnIn)/l)
            ##function calculation
            resize(blockMV, h*l) ## size our model value object to be
            ## approximately of size m (number of mc
            ## samples)

            for(r in 1:NbootReps){
                for(i in 1:h){
                    randNum <- runif(1,0,1)
                    ##random starting index for blocks (post burn-in)
                    randIndex <- ceiling(randNum*q)
                    for(j in 1:l){
                        ## fill in blockMV with chosen blocks
                        nimCopy(mcmcMV, blockMV, paramNames, paramNames,
                                burnIn + randIndex - 1 + j,  (i - 1)*l + j)
                    }
                }
                ##as per Caffo, calculate both Q functions using the same
                ##samples from the latent variables
                blockpppValues[r]  <- run(NpostSamp, TRUE)
            }

            blockpppValuesSD <- sd(blockpppValues)
            returnType(double())
            return(blockpppValuesSD)
        }
    ))




## calculate proportion of deviances that are greater or equal to the
## observed value
calcCPPP <- function(MCMCIter,
                     burnIn,
                     NpostSamp,
                     C.pppFunc,
                     cppp.C.mcmc,
                     firstRun,
                     runUntilConverged = NULL,
                     maxIter = NULL,
                     convStep= NULL){
    convergeTest <- NA
    samples <- NA
    if(firstRun  == FALSE){
        cppp.C.mcmc$run(MCMCIter)
        samples <- mcmc(as.matrix(cppp.C.mcmc$mvSamples)[-c(1:burnIn),])
        convergeTest <- geweke.diag(samples)$z
        if(runUntilConverged == TRUE){
            ## use Geweke diagnostic to see if MCMC has converged
            ## any na values get set to a very large z value so the MCMC
            ## will continue to be run. note that nan values correspond to
            ## posterior samples that are constant, indicating a lack of
            ## convergence
            convergeTest[!is.finite(convergeTest)] <- 10
            zVal <- 1.96
            while(any(abs(convergeTest) > zVal)){
                cppp.C.mcmc$run(MCMCIter*convStep, reset = FALSE)
                this.niter <- MCMCIter + convStep*MCMCIter
                samples <- mcmc(as.matrix(cppp.C.mcmc$mvSamples)[-c(1:burnIn),])
                convergeTest <- geweke.diag(samples)$z
                convergeTest[!is.finite(convergeTest)] <- 10
                if(this.niter > maxIter)
                    convergeTest <- 0
            }
        }
    }
    pre.pp <- C.pppFunc$run(NpostSamp, FALSE)
    bootSD <- C.pppFunc$getSD(NpostSamp)
    if(!is.finite(pre.pp))    pre.pp <- NA
    return(list(pre.pp = pre.pp,
                bootSD = bootSD,
                samples = samples,
                converge.stat = convergeTest))
}


generateCPPP <-  function(R.model,
                          orig.C.model,
                          orig.C.mcmc,
                          orig.mcmc,
                          dataNames, ## names of the data column
                          paramNames, ## vector of parameters to monitor
                          NpostSamp,## number of samples from posterior
                          NPDist, ## number of simulated PPP values
                          burnInProp, ## proportion of mcmc to drop
                          discFuncGenerator,
                          averageParams,
                          returnChains = TRUE,
                          runUntilConverged = FALSE,
                          maxIter = 1*10^4,
                          convStep = 0.5,
                          nRepBoot, ## number of bootstrap samples
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

  if(NpostSamp > MCMCIter){
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


  ## keep track of the real data
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
                          discFunc = discFuncGenerator,
                          NbootReps=nRepBoot)
  C.pppFunc <- compileNimble(modelpppFunc,
                             project = R.model)


  ## calculate deviances
  obs.cppp <- calcCPPP(MCMCIter,
                       burnIn,
                       NpostSamp=NpostSamp,
                       C.pppFunc,
                       orig.C.mcmc,
                       firstRun = TRUE)

  ## refits model with sampled data, reruns, enter inner loop,
  ## calculates distbution of PPPs
  simPppDist <- function(iteration){
    message(paste("refitting data iteration", iteration))
    simulate(orig.C.model,  includeData =  TRUE)
    out <- calcCPPP(MCMCIter,
                    burnIn,
                    NpostSamp,
                    C.pppFunc,
                    orig.C.mcmc,
                    firstRun = FALSE,
                    runUntilConverged=runUntilConverged,
                    maxIter = maxIter,
                    convStep=convStep)
    return(out)
  }
  
  ## simulate the ppp values
  sim.ppp.output <- mclapply(1:NPDist, simPppDist)

  ## extract simulated ppp and boot SDs
  sim.ppp <- sapply(sim.ppp.output,
                    function(x) x$pre.pp)
  sim.vars <- (sapply(sim.ppp.output,
                      function(x) x$bootSD))^2 

  ## approximate the distbution of observed ppp and simulated ppp
  diff.ppp <- obs.cppp$pre.pp - sim.ppp 
  diff.vars.ppp <- obs.cppp$bootSD^2 + sim.vars

  ## simulate the ppp values
  sim.ppp.output <- mclapply(1:NPDist, simPppDist)

  ## extract simulated ppp and boot SDs
  sim.ppp <- sapply(sim.ppp.output,
                    function(x) x$pre.pp)
  sim.vars <- (sapply(sim.ppp.output,
                      function(x) x$bootSD))^2

  ## approximate the distbution of observed ppp and simulated ppp
  diff.ppp <- obs.cppp$pre.pp - sim.ppp
  diff.vars.ppp <- obs.cppp$bootSD^2 + sim.vars

  ## calculate the average prob that simulated values are less than
  ## the observed
  cppp <- mean(pnorm(0, diff.ppp, diff.vars.ppp),
               na.rm=TRUE)

  ## extract chains and diagnostics
  sim.samples <- lapply(sim.ppp.output,
                        function(x) x$samples)
  chain.diag <- sapply(sim.ppp.output,
                       function(x) x$converge.stat)

  ## output real data and model
  nimble:::values(orig.C.model, dataNames) <- origData

  out <- list(cppp=cppp,
              obs.ppp=c(estimate=obs.cppp$pre.pp,
                bootVar=(obs.cppp$bootSD)^2),
              sim.cpp.dist=cbind(esimate=sim.ppp,
                bootVar=sim.vars),
              chain.diagnostics= chain.diag,
              samples=sim.samples)
  if(!returnChains){
    out$samples <- NULL
  }
  return(out)
}

